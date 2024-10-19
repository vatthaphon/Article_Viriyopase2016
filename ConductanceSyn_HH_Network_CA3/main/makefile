GSL_HOME 		= /home/viriyopa/utils/gsl
MY_HOME			= /home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3

GSL_INCLUDE 	= $(GSL_HOME)/include
GSL_LIB 		= $(GSL_HOME)/lib

MY_INCLUDE 		= $(MY_HOME)/include
MY_MAIN			= $(MY_HOME)/main
MY_MISC			= $(MY_HOME)/misc

OUTPUT = hh_network_PING_ING_WB
#OUTPUT = hh_network_PING_ING_WB_incdec_I
#OUTPUT = hh_network_PING_ING_WB_PERIODIC_EXTERNAL_INPUTS
#OUTPUT = hh_network_PING_ING_BW
#OUTPUT = hh_network_PING_ING_BW_PulseSyn
#OUTPUT = hh_network_PING_ING_BW_incdec_I
#OUTPUT = hh_network_PING_ING_WHITE_WB

INPUT = main

all: 
	make -f $(MY_MISC)/makefile
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MAIN)/$(INPUT).c -o $(MY_MAIN)/$(INPUT).o
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -L$(MY_MISC) -L$(GSL_LIB) -Werror -Wall -o $(MY_MAIN)/$(OUTPUT) $(MY_MAIN)/$(INPUT).o -lcsln -lgsl -lgslcblas  -lm 
