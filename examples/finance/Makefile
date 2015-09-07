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

INCLUDES=heston.hpp heston_pde.hpp option.hpp market.hpp result.hpp

LFLAGS=


PROGS=example

all:		${PROGS}

example:	example.o heston.o heston_analytic.o
		${CPP} example.o heston.o heston_analytic.o $(LFLAGS) -o example

%.o: %.cpp	${INCLUDES}
		${CPP} ${CFLAGS} -c $< -o $@

%.o: %.c
		${CC} ${CFLAGS} -c $< -o $@
clean:
		rm -f *.o ${PROGS}
