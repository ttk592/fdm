CPP=g++
CC=gcc
AR=ar


CFLAGS=-O3 -DNDEBUG
#CFLAGS=-g 
#CFLAGS=-g -pg
CFLAGS += -Wall -Wextra  
CFLAGS += -I../../src
CFLAGS += -DUSE_FENV
CFLAGS += -std=c++11
#CFLAGS += -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0

INCLUDES=pde1d.h

LFLAGS=

OUT=exec


PROGS=pde1d

all:		${PROGS}

pde1d:		pde1d.o 
		${CPP} pde1d.o $(LFLAGS) -o pde1d

%.o: %.cpp	${INCLUDES}
		${CPP} ${CFLAGS} -c $< -o $@

%.o: %.c
		${CC} ${CFLAGS} -c $< -o $@
clean:
		rm -f *.o ${PROGS}
