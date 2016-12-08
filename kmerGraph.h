#ifndef _KMERGRAPH_H_
#define _KMERGRAPH_H_

#include <map>
#include <string>
#include <ext/hash_map>
#include <unistd.h>
#include "sequence.h"
#include "mympi.h"

using namespace std;

class kmerGraph
{
public:
	unsigned long long *kmers;
	unsigned char *arcs;
	unsigned long long size;
	unsigned long long read_pos;
	double commtime;
	MPI_Datatype commType;
	
	unsigned long long reverseComplement(unsigned long long kmerDescriptor, parameter *parameters);
    	static unsigned long long stringToLongLong(const char *buf, int start, int end, parameter *parameters);
	static string longLongToString(unsigned long long a, parameter *parameters);
	unsigned long long getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol);
	int arcPos(unsigned long long &A, unsigned long long &B, char directA, char directB, int hashLength);

	kmerGraph()
	{
		kmers    = NULL;
		arcs     = NULL;
		size     = 0;
		read_pos = 0;
		commtime = 0;

		MPI_Type_contiguous(KMER_COMM_TYPE_LEN, MPI_UNSIGNED_LONG_LONG, &commType);
		MPI_Type_commit(&commType);
	}
	
	~kmerGraph()
	{
		delete(kmers);
		delete(arcs);
		kmers = NULL;
		arcs = NULL;
	}

public:
	int  constructKmerGraph(sequence *sequences, parameter *parameters, MPIEnviroment *MPIcontrol);
	void printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void distributeKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
};

#endif
