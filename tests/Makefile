CPP=g++
CC=gcc
AR=ar


CFLAGS=-O2
CFLAGS +=-g
#CFLAGS=-g -pg
CFLAGS += -Wall -Wextra
CFLAGS += -I../src
CFLAGS += -DUSE_FENV
#CFLAGS += -DDEBUG
#CFLAGS += -DNDEBUG 
CFLAGS += -std=c++11

INCLUDES=

LFLAGS=-lboost_unit_test_framework

OUT=exec


PROGS=interpolation svector grid grid_interpol algorithm

all:		${PROGS}



# cancelling default rules as they would short cut assembly and object stages
% : %.cpp


%: %.o
		${CPP} $< $(LFLAGS) -o $@


%.o: %.cpp	${INCLUDES}
		${CPP} ${CFLAGS} -c $< -o $@


# mark .o files as precious,
# otherwise would be deleted as considered intermediates
.PRECIOUS: %.o 


clean:
		rm -f *.o ${PROGS} interpol.csv
