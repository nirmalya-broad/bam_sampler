# compiler to use
CC = c++
CFLAGS=-O3 -std=c++14 -g

LOCALPATH=/broad/IDP-Dx_work/nirmalya/local
INC = -I${LOCALPATH}/include
LIBDIR=${LOCALPATH}/lib/
#LIBS=${LIBDIR}/libhts.so  -lz -lboost_program_options -lboost_regex
LIBS=-lhts  -lz -lboost_program_options -lboost_regex


all: tools
    
tools:
	$(CC) $(CFLAGS) $(INC) sample_bam.cpp -o sample_bam  $(LIBS) 
	$(CC) $(CFLAGS) $(INC) remove_bu.cpp -o remove_bu  $(LIBS) 
	$(CC) $(CFLAGS) $(INC) seed_gen.cpp -o seed_gen  $(LIBS) 
	$(CC) $(CFLAGS) $(INC) frag_count.cpp -o frag_count  $(LIBS) 

