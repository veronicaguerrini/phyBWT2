#CCLIB=-fopenmp -lsdsl -ldivsufsort -ldivsufsort64
VLIB= -g -O0

#LIB_DIR = ${HOME}/lib
#INC_DIR = ${HOME}/include
MY_CXX_FLAGS= -std=c++11 -Wall -fomit-frame-pointer -Wno-comment
#-D_FILE_OFFSET_BITS=64

MY_CXX_OPT_FLAGS= -O3 -m64
#MY_CXX_OPT_FLAGS= $(VLIB)
MY_CXX=g++


##

LIBOBJ = \
	external/malloc_count/malloc_count.o

##

M64 = 0
DEBUG = 0
SHORT = 1
##

LFLAGS = -lm -ldl

DEFINES = -DDEBUG=$(DEBUG) -DM64=$(M64) -DSHORT=$(SHORT) -DminLCP=$(minLCP)

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(LFLAGS) $(DEFINES) -I$(INC_DIR) -L$(LIB_DIR)

##

all: compile

clean:
	\rm -f *.o  external/*.o external/malloc_count/*.o phyBWT2 create_cda

##

compile: phyBWT2 create_cda

phyBWT2: src/phyBWT2.cpp ${LIBOBJ} 
	$(MY_CXX) src/phyBWT2.cpp $(CCLIB) -o phyBWT2 ${LIBOBJ} $(CXX_FLAGS) 

create_cda: src/create_cda.cpp ${LIBOBJ} 
	$(MY_CXX) src/create_cda.cpp $(CCLIB) -o create_cda ${LIBOBJ} $(CXX_FLAGS) 
