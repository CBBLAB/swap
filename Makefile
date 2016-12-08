CC=mpic++

all	:swap stats
swap    : graph.o kmerGraph.o sequence.o mympi.o
	$(CC) -O2 graph.o kmerGraph.o sequence.o mympi.o  -o swap
stats   :stats.cpp
	$(CC) stats.cpp -o stats	

graph.o : graph.cpp kmerGraph.h mympi.h sequence.h
	$(CC) -Wno-deprecated -c -O2 graph.cpp -o graph.o
kmerGraph.o : kmerGraph.cpp kmerGraph.h mympi.h sequence.h
	$(CC) -Wno-deprecated -c -O2 kmerGraph.cpp -o kmerGraph.o
sequence.o : sequence.cpp sequence.h mympi.h
	$(CC) -Wno-deprecated -c -O2 sequence.cpp -o sequence.o
mympi.o : mympi.cpp mympi.h
	$(CC) -Wno-deprecated -c -O2 mympi.cpp -o mympi.o

clean :
	rm -f stats swap graph.o kmerGraph.o sequence.o mympi.o

#	$(CC) -c -Wno-deprecated -g graph.cpp -o graph.o
