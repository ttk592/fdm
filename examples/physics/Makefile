CPP=g++
CC=gcc
AR=ar

CFLAGS=-O3 -DNDEBUG
#CFLAGS=-g 
#CFLAGS=-g -pg
CFLAGS += -Wall -Wextra -I../../src
CFLAGS += -DUSE_FENV
CFLAGS += -std=c++11
#CFLAGS += -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0

INCLUDES=

LFLAGS=

OUT=exec


PROGS=heat

all:		${PROGS}


heat:		heat.o png_io.o
		${CPP} heat.o png_io.o $(LFLAGS) -lpng -o heat


%.o: %.cpp	${INCLUDES}
		${CPP} ${CFLAGS} -c $< -o $@

%.o: %.c
		${CC} ${CFLAGS} -c $< -o $@
clean:
		rm -f *.o ${PROGS}
