# Makefile for ACO++

CPP= g++
CFLAGS= -O3 -w -std=c++11
EXE=acothop

all: thop.o utilities.o ants.o inout.o unix_timer.o ls.o parse.o acothop.o 
	$(CPP) $(CFLAGS) -o $(EXE) thop.o utilities.o ants.o inout.o unix_timer.o ls.o parse.o acothop.o 

thop.o: thop.c thop.h
	$(CPP) $(CFLAGS) -c thop.c

utilities.o: utilities.c utilities.h
	$(CPP) $(CFLAGS) -c utilities.c

ants.o: ants.c ants.h
	$(CPP) $(CFLAGS) -c ants.c

inout.o: inout.c inout.h
	$(CPP) $(CFLAGS) -c inout.c

unix_timer.o: unix_timer.c timer.h
	$(CPP) $(CFLAGS) -c unix_timer.c

ls.o: ls.c ls.h
	$(CPP) $(CFLAGS) -c ls.c

parse.o: parse.c parse.h
	$(CPP) $(CFLAGS) -c parse.c

acothop.o: ants.h utilities.h inout.h thop.h timer.h ls.h acothop.c
	$(CPP) $(CFLAGS) -c acothop.c

clean:
	rm -f $(EXE) *.o
