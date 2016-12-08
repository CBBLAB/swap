#include "sequence.h"
#include "mympi.h"
#include "kmerGraph.h"

unsigned long long kmerGraph::reverseComplement(unsigned long long a, parameter *parameters)
{
        unsigned long long rev = 0;
        for(int i=0; i<parameters->hashLength; i++){
                rev <<= 2;
                rev |= (3-a&3);
                a >>= 2;
        }
        return rev;
}

unsigned long long kmerGraph::stringToLongLong(const char *buf, int start, int end, parameter *parameters)
{
        unsigned long long ret = 0;
        for(int i=start; i<end; i++){
                ret <<= 2;
                ret |= (unsigned long long)parameters->nucleotideValue[buf[i]]&3;
        }
        return ret;
}

string kmerGraph::longLongToString(unsigned long long a, parameter *parameters)
{
	string descriptor;
        descriptor.clear();
        for(int i=0;i<parameters->hashLength;i++)
        {
                descriptor += parameters->nucleotideArray[a%4];
                a = a>>2;
        }
        reverse(descriptor.begin(), descriptor.end());
	return descriptor;
}

int kmerGraph::arcPos(unsigned long long &A, unsigned long long &B, char directA, char directB, int hashLength)
{
	int lastB;
	if(directB=='+')	lastB = B%4;
	else			lastB = 3-(B>>((hashLength-1)*2));
	
	if(directA=='-')	lastB += 4;
	return lastB;
}

/*
unsigned long long kmerGraph::getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{
	unsigned long long ret = 0, tmpID = kmerID;
	for(int i=0;i<hashLength;i++)
	{
		ret = ret * 4;
//		ret |= ((unsigned long long)3 - (tmpID&(unsigned long long)3));
		ret += (3 - tmpID%4);
		tmpID = tmpID / 4;	
	}
//	ret = ret ^ kmerID;
//	return ( (ret % (unsigned long long) MersPrime) % (unsigned long long) MPIcontrol->nprocs);	
//	return ( (ret ) % (unsigned long long) MPIcontrol->nprocs);	
	if(ret<kmerID)	ret = kmerID;
	
	double tmp = ((sqrt(5.0)-1)/2) * ret ;
	double rr  = tmp - floor(tmp);
	return floor(rr * MPIcontrol->nprocs);
}
*/

unsigned long long kmerGraph::getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{
        unsigned long long revKmer = 0, tmpID = kmerID;
        for(int i=0;i<hashLength;i++)
        {
                revKmer = revKmer * 4;
                revKmer += ( 3 - tmpID%4 );
                tmpID = tmpID / 4;
        }

        if(revKmer>kmerID) revKmer = kmerID;

         unsigned int factor = 19;
         unsigned int numBytes = (hashLength + 3) / 4;

         unsigned int sum = 0;
         for(int i = 0; i < numBytes; i++){
                 sum = sum * factor + (revKmer & 0xFF);
                 revKmer >>= 8;
         }
         return sum % MPIcontrol->nprocs;
}


int kmerGraph::constructKmerGraph(sequence *sequences, parameter *parameters, MPIEnviroment *MPIcontrol)
{
	unsigned long long tot_kmer_num = 0, totSize=0;

	unsigned long long readStep = MPIcontrol->nprocs/4;
//	unsigned long long readStep = REA_PROCESS_STEP/MPIcontrol->nprocs;
   	for(int i=0; i<readStep&& (read_pos+i) < sequences->readCount; i++)
	{
		int len = strlen(sequences->reads[read_pos+i]);
		if(len<parameters->hashLength) continue;

		len = len - parameters->hashLength;
		tot_kmer_num += len;	
	}	
	
	size = 2*tot_kmer_num;
	kmers = new unsigned long long 	[size];
	arcs  = new unsigned char 	[size];

	unsigned long long kmer_index=0;
	for(int i=0; i<readStep && (read_pos+i)<sequences->readCount; i++)
        {
         	if((read_pos+i)%100000==0) 
		{	
			char localstr[100];
			sprintf(localstr, "construct graph: %d (readCount=%d)", read_pos+i, sequences->readCount);
			MPIcontrol->print(localstr);
		}

                char *buf = sequences->reads[read_pos+i];
                unsigned long long curNodeID=0,     primeNodeID=0,   	twinNodeID=0;
		unsigned long long preNodeID=0,     prePrimeNodeID=0; 

                int len = strlen(buf);
                if(len < parameters->hashLength) continue;
                curNodeID = this->stringToLongLong(buf, 0, parameters->hashLength-1, parameters);
                for(int j=0; j<=len-parameters->hashLength; j++) {
                	curNodeID <<= 2;
                	curNodeID |= (parameters->nucleotideValue[buf[j+parameters->hashLength-1]]&(3ull));
                	curNodeID &= parameters->MASK;

                	twinNodeID  = this->reverseComplement(curNodeID, parameters);
                	primeNodeID = (curNodeID>twinNodeID) ? curNodeID : twinNodeID;

			if(j==0)	{
				preNodeID = curNodeID;
				prePrimeNodeID = primeNodeID; 
				continue;
			}

			kmers[kmer_index]     = prePrimeNodeID;	arcs[kmer_index] = 0;
			kmers[kmer_index+1]   = primeNodeID;	arcs[kmer_index+1]   = 0;
			
			int pos_ret;
			//preNode -> curNode  add(A, B, +, +, last(B+))
			if( preNodeID == prePrimeNodeID && curNodeID == primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '+', '+', parameters->hashLength);	
				arcs[kmer_index]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '-', '-', parameters->hashLength); 
				arcs[kmer_index+1] |=  ((unsigned char) 1<<pos_ret);
			}
			//preNode -> ~curNode
			else if(preNodeID == prePrimeNodeID && curNodeID != primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '+', '-', parameters->hashLength);	
				arcs[kmer_index]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '+', '-', parameters->hashLength); 
				arcs[kmer_index+1] |=  ((unsigned char) 1<<pos_ret);
			}
			//~preNode -> curNode
			else if(preNodeID != prePrimeNodeID && curNodeID == primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '-', '+', parameters->hashLength);	
				arcs[kmer_index]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '-', '+', parameters->hashLength); 
				arcs[kmer_index+1] |=  ((unsigned char) 1<<pos_ret);
			}			
			//~preNode -> ~curNode
			else if(preNodeID != prePrimeNodeID && curNodeID != primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '-', '-', parameters->hashLength);	
				arcs[kmer_index]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '+', '+', parameters->hashLength); 
				arcs[kmer_index+1] |=  ((unsigned char) 1<<pos_ret);
			}
			preNodeID = curNodeID;
			prePrimeNodeID = primeNodeID;
			kmer_index += 2;
		}
	//	printf("--------+++++++++++++++++++++++++++''''''''''''''''\n\n");

                delete sequences->reads[read_pos+i];
		sequences->reads[read_pos+i] = NULL;
	}
	assert(kmer_index == size);

	read_pos = read_pos+readStep; 
   	if(read_pos > sequences->readCount)	read_pos = sequences->readCount;

	unsigned long long leftReads = sequences->readCount - read_pos;
//	unsigned long long totLeftReads; 
	unsigned long long totkmers=0;	

	MPI_Allreduce(&kmer_index, &totkmers, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD); 

//      	printf("Proc %d: kmer_index=%llu, leftReads=%llu\n", MPIcontrol->rank, kmer_index, leftReads);
	return totkmers;
}

void kmerGraph::printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	/*hash_map<unsigned long long, kmer>::iterator it;	
	for(it=kmers.begin();it!=kmers.end();it++)
	{
	
		kmer tmp = it->second;
		string descriptor = this->longLongToString(tmp.kmerID, parameters);		
		printf("proc:%d ID=%lld, Descriptor=%s\n",MPIcontrol->rank, it->first, descriptor.c_str());
		for(int i=0;i<4;i++)
		{
			if(tmp.arcs&(1<<i))
			printf("proc:%d -(%c|%d)-",MPIcontrol->rank, parameters->nucleotideArray[i], tmp.multiplicity[i]); 
		} 
		printf("\n");

		for(int i=4;i<8;i++)
		{
			if(tmp.arcs&(1<<i))
			printf("proc:%d -(%c|%d)-",MPIcontrol->rank, parameters->nucleotideArray[i-4], tmp.multiplicity[i]); 
		} 
		printf("\nproc:%d--------------------------------------\n", MPIcontrol->rank);
	}*/
}

void kmerGraph::distributeKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	unsigned long long *kmers_send = new unsigned long long [size];
	assert(kmers_send != NULL);
	unsigned char      *arcs_send  = new unsigned char      [size];
	assert(arcs_send != NULL);
	for(unsigned long long i=0; i < size; i++)	kmers_send[i] = arcs_send[i] = 0;

	int nproc = MPIcontrol->nprocs;
	
//****	
	int  *send_size = new int [nproc];
	assert(send_size!=NULL);
	int  *recv_size = new int [nproc];
	assert(recv_size!=NULL);
	int  *send_pos  = new int [nproc];
	assert(send_pos!=NULL);
	int  *recv_pos  = new int [nproc];
	assert(recv_pos!=NULL);

	for(int i=0; i<nproc; i++)	send_size[i] = recv_size[i] = 0;

	unsigned long long send_sum=0,  recv_sum=0;
	unsigned long long ProcsID;	

	for(unsigned long long i=0; i<size; i++)
	{
		ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
		send_size[ProcsID]++;
	}

//****	
	MPI_Alltoall(send_size,1, MPI_INT, recv_size, 1, MPI_INT, MPI_COMM_WORLD);

	for(int i=0; i<nproc; i++)	send_sum += send_size[i];
	for(int i=0; i<nproc; i++)	recv_sum += recv_size[i];
	assert(send_sum == size);
	
//****	
//	printf("Proc%d: send_sum=%llu recv_sum=%llu\n", MPIcontrol->rank, send_sum, recv_sum);

	send_pos[0] = 0;
	for(int i=1; i<nproc; i++)	send_pos[i] = send_pos[i-1] + send_size[i-1];
	recv_pos[0] = 0;
	for(int i=1; i<nproc; i++)	recv_pos[i] = recv_pos[i-1] + recv_size[i-1];

/*
        printf("proc%d: ", MPIcontrol->rank);
        for(int i=0;i<size;i++) printf("send_size[%d]=%d ", i, send_size[i]);
        printf("\n");
	printf("proc%d: ", MPIcontrol->rank);
        for(int i=0;i<size;i++) printf("recv_size[%d]=%d ", i, recv_size[i]);
        printf("\n");
*/
	for(int i=0; i<size; i++)
	{
		ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
		kmers_send[send_pos[ProcsID]]  = kmers[i]; 
		arcs_send[send_pos[ProcsID]++] = arcs[i];	
	}

	send_pos[0] = 0;
	for(int i=1; i<nproc; i++)	send_pos[i] = send_pos[i-1] + send_size[i-1];

	delete kmers;
	delete arcs;
	kmers = new unsigned long long [recv_sum];
	assert(kmers!=NULL);
	arcs  = new unsigned char      [recv_sum];
	assert(arcs!=NULL);


	clock_t t1 = clock();
	MPI_Alltoallv(kmers_send,send_size,send_pos, MPI_LONG_LONG_INT, kmers,recv_size,recv_pos, MPI_LONG_LONG_INT, MPI_COMM_WORLD);	
	MPI_Alltoallv(arcs_send, send_size,send_pos, MPI_CHAR,          arcs, recv_size,recv_pos, MPI_CHAR,          MPI_COMM_WORLD);
	clock_t t2 = clock();
	commtime = commtime + t2-t1;
/*
	MPI_Status *status     = new MPI_Status  [nproc];
	MPI_Request *send_reqs = new MPI_Request [nproc];
	MPI_Request *recv_reqs = new MPI_Request [nproc];

	
//	for(int i=0;i<send_sum;i++)	printf("proc:%d send_kmer[%d]=%ld\n", MPIcontrol->rank, i, send_buf[i].kmerID);

	for(int i=0;i<nproc;i++)
		MPI_Irecv(kmers+recv_pos[i], recv_size[i], MPI_LONG_LONG_INT, i, MPIcontrol->rank, MPI_COMM_WORLD, &recv_reqs[i]);

	for(int i=0;i<nproc;i++)
	{
		MPI_Send(kmers_send+send_pos[i], send_size[i], MPI_LONG_LONG_INT, i,i,MPI_COMM_WORLD);
	}
	MPI_Waitall(nproc, recv_reqs, status);
//	MPI_Waitall(nproc, send_reqs, status);
//	MPI_Waitall(nproc, send_reqs, status);

	for(int i=0;i<nproc;i++)
		MPI_Irecv(arcs+recv_pos[i], recv_size[i], MPI_CHAR, i, MPIcontrol->rank, MPI_COMM_WORLD, &recv_reqs[i]);

	for(int i=0;i<nproc;i++)
		MPI_Send(arcs_send+send_pos[i], send_size[i], MPI_CHAR, i,i,MPI_COMM_WORLD);

	MPI_Waitall(nproc, recv_reqs, status);
//	MPI_Waitall(nproc, send_reqs, status);
	
//	for(int i=0;i<recv_sum;i++)	printf("proc:%d recv_kmer[%d]=%ld\n", MPIcontrol->rank, i, recv_buf[i].kmerID);

	delete status;
	delete send_reqs;
	delete recv_reqs;
*/
	delete kmers_send;
	delete arcs_send;
	kmers_send = NULL;
	arcs_send = NULL;
	size  = recv_sum;
	
	delete send_size;
	delete recv_size;
	delete send_pos;
	delete recv_pos;
	send_size = recv_size = send_pos = recv_pos = NULL;
}
/*
int main (int argc, char *argv[])
{
   	int i, size, namelen, lgsize;
	MPIEnviroment MPIcontrol;
	MPIcontrol.init(argc, argv);

	parameter parameters;
	sequence sequences;
	parameters.getParameters(argc, argv);	
	sequences.getSequences(&parameters, &MPIcontrol);

	kmerGraph mygraph;

	while(mygraph.constructKmerGraph(&sequences, &parameters, &MPIcontrol)!=0)
	{
		//mygraph.printKmerGraph(&parameters, &MPIcontrol);

	//	char message[100];
		//sprintf(message,"kmer number is %d\n", mygraph.kmers.size());	
		//MPIcontrol.print(message);

		mygraph.distributeKmerGraph(&parameters, &MPIcontrol);	

		delete mygraph.kmers;
		delete mygraph.arcs;
	
	//	sprintf(message, "distribute kmers finished\n");
	//	MPIcontrol.print(message);
	}

	MPIcontrol.finalize();
    	return (0);
}*/
