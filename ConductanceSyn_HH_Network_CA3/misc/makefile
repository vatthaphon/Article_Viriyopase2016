GSL_HOME 		= /home/viriyopa/utils/gsl
MY_HOME			= /home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3

MY_INCLUDE 		= $(MY_HOME)/include
MY_MAIN			= $(MY_HOME)/main
MY_MISC 		= $(MY_HOME)/misc
GSL_INCLUDE 	= $(GSL_HOME)/include
GSL_LIB 		= $(GSL_HOME)/lib

all: interpolate_White_V_c interpolate_White_T0_c interpolate_White_n_c interpolate_White_h_c interpolate_White_t2spk_c interpolate_CA1_h_CaH_c interpolate_CA1_m_KDR_c interpolate_CA1_m_KM_c interpolate_CA1_h_CaT_c interpolate_CA1_h_KDR_c interpolate_CA1_h_NaT_c interpolate_CA1_m_CaH_c interpolate_CA1_m_CaT_c interpolate_CA1_V_c interpolate_CA1_T0_c interpolate_CA3_h_CaH_c interpolate_CA3_m_KDR_c interpolate_CA3_m_KM_c interpolate_CA3_h_CaT_c interpolate_CA3_h_KDR_c interpolate_CA3_h_NaT_c interpolate_CA3_m_CaH_c interpolate_CA3_m_CaT_c interpolate_CA3_V_c interpolate_CA3_T0_c interpolate_BW_T0_c interpolate_BW_n_c interpolate_BW_h_c interpolate_BW_V_c interpolate_WB_V_c interpolate_WB_T0_c interpolate_WB_n_c interpolate_WB_h_c interpolate_WB_t2spk_c interpolate_CA3_t2spk_c interpolate_CA1_t2spk_c m_time_c add_Input_c cmp_c matrix_c brunel_wang_2003_network_c create_static_lib

create_static_lib:
	ar -cr $(MY_MISC)/libcsln.a $(MY_MISC)/interpolate_White_V.o $(MY_MISC)/interpolate_White_T0.o $(MY_MISC)/interpolate_White_n.o $(MY_MISC)/interpolate_White_h.o $(MY_MISC)/interpolate_White_t2spk.o $(MY_MISC)/interpolate_CA1_h_CaH.o $(MY_MISC)/interpolate_CA1_m_KDR.o $(MY_MISC)/interpolate_CA1_m_KM.o $(MY_MISC)/interpolate_CA1_h_CaT.o $(MY_MISC)/interpolate_CA1_h_KDR.o $(MY_MISC)/interpolate_CA1_h_NaT.o $(MY_MISC)/interpolate_CA1_m_CaH.o $(MY_MISC)/interpolate_CA1_m_CaT.o $(MY_MISC)/interpolate_CA1_V.o $(MY_MISC)/interpolate_CA1_T0.o $(MY_MISC)/interpolate_CA3_h_CaH.o $(MY_MISC)/interpolate_CA3_m_KDR.o $(MY_MISC)/interpolate_CA3_m_KM.o $(MY_MISC)/interpolate_CA3_h_CaT.o $(MY_MISC)/interpolate_CA3_h_KDR.o $(MY_MISC)/interpolate_CA3_h_NaT.o $(MY_MISC)/interpolate_CA3_m_CaH.o $(MY_MISC)/interpolate_CA3_m_CaT.o $(MY_MISC)/interpolate_CA3_V.o $(MY_MISC)/interpolate_CA3_T0.o $(MY_MISC)/interpolate_BW_T0.o $(MY_MISC)/interpolate_BW_n.o $(MY_MISC)/interpolate_BW_h.o $(MY_MISC)/interpolate_BW_V.o $(MY_MISC)/interpolate_WB_V.o $(MY_MISC)/interpolate_WB_T0.o $(MY_MISC)/interpolate_WB_n.o $(MY_MISC)/interpolate_WB_h.o $(MY_MISC)/interpolate_CA1_t2spk.o $(MY_MISC)/interpolate_CA3_t2spk.o $(MY_MISC)/interpolate_WB_t2spk.o $(MY_MISC)/brunel_wang_2003_network.o $(MY_MISC)/m_time.o $(MY_MISC)/add_Input.o $(MY_MISC)/cmp.o $(MY_MISC)/matrix.o

interpolate_CA1_h_KDR_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_h_KDR.c -o $(MY_MISC)/interpolate_CA1_h_KDR.o
interpolate_CA1_h_CaH_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_h_CaH.c -o $(MY_MISC)/interpolate_CA1_h_CaH.o
interpolate_CA1_m_KDR_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_m_KDR.c -o $(MY_MISC)/interpolate_CA1_m_KDR.o
interpolate_CA1_m_KM_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_m_KM.c -o $(MY_MISC)/interpolate_CA1_m_KM.o
interpolate_CA1_h_CaT_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_h_CaT.c -o $(MY_MISC)/interpolate_CA1_h_CaT.o
interpolate_CA1_h_NaT_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_h_NaT.c -o $(MY_MISC)/interpolate_CA1_h_NaT.o
interpolate_CA1_m_CaH_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_m_CaH.c -o $(MY_MISC)/interpolate_CA1_m_CaH.o
interpolate_CA1_m_CaT_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_m_CaT.c -o $(MY_MISC)/interpolate_CA1_m_CaT.o
interpolate_CA1_V_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_V.c -o $(MY_MISC)/interpolate_CA1_V.o
interpolate_CA1_T0_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_T0.c -o $(MY_MISC)/interpolate_CA1_T0.o

interpolate_CA3_h_KDR_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_h_KDR.c -o $(MY_MISC)/interpolate_CA3_h_KDR.o
interpolate_CA3_h_CaH_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_h_CaH.c -o $(MY_MISC)/interpolate_CA3_h_CaH.o
interpolate_CA3_m_KDR_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_m_KDR.c -o $(MY_MISC)/interpolate_CA3_m_KDR.o
interpolate_CA3_m_KM_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_m_KM.c -o $(MY_MISC)/interpolate_CA3_m_KM.o
interpolate_CA3_h_CaT_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_h_CaT.c -o $(MY_MISC)/interpolate_CA3_h_CaT.o
interpolate_CA3_h_NaT_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_h_NaT.c -o $(MY_MISC)/interpolate_CA3_h_NaT.o
interpolate_CA3_m_CaH_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_m_CaH.c -o $(MY_MISC)/interpolate_CA3_m_CaH.o
interpolate_CA3_m_CaT_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_m_CaT.c -o $(MY_MISC)/interpolate_CA3_m_CaT.o
interpolate_CA3_V_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_V.c -o $(MY_MISC)/interpolate_CA3_V.o
interpolate_CA3_T0_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_T0.c -o $(MY_MISC)/interpolate_CA3_T0.o

interpolate_WB_V_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_WB_V.c -o $(MY_MISC)/interpolate_WB_V.o
interpolate_WB_T0_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_WB_T0.c -o $(MY_MISC)/interpolate_WB_T0.o
interpolate_WB_n_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_WB_n.c -o $(MY_MISC)/interpolate_WB_n.o
interpolate_WB_h_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_WB_h.c -o $(MY_MISC)/interpolate_WB_h.o
interpolate_WB_t2spk_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_WB_t2spk.c -o $(MY_MISC)/interpolate_WB_t2spk.o

interpolate_White_V_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_White_V.c -o $(MY_MISC)/interpolate_White_V.o
interpolate_White_T0_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_White_T0.c -o $(MY_MISC)/interpolate_White_T0.o
interpolate_White_n_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_White_n.c -o $(MY_MISC)/interpolate_White_n.o
interpolate_White_h_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_White_h.c -o $(MY_MISC)/interpolate_White_h.o
interpolate_White_t2spk_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_White_t2spk.c -o $(MY_MISC)/interpolate_White_t2spk.o
	
interpolate_BW_T0_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_BW_T0.c -o $(MY_MISC)/interpolate_BW_T0.o	
interpolate_BW_V_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_BW_V.c -o $(MY_MISC)/interpolate_BW_V.o	
interpolate_BW_h_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_BW_h.c -o $(MY_MISC)/interpolate_BW_h.o	
interpolate_BW_n_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_BW_n.c -o $(MY_MISC)/interpolate_BW_n.o	

interpolate_CA3_t2spk_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA3_t2spk.c -o $(MY_MISC)/interpolate_CA3_t2spk.o

interpolate_CA1_t2spk_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/interpolate_CA1_t2spk.c -o $(MY_MISC)/interpolate_CA1_t2spk.o

brunel_wang_2003_network_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -c $(MY_MISC)/brunel_wang_2003_network.c -o $(MY_MISC)/brunel_wang_2003_network.o

m_time_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/m_time.c -o $(MY_MISC)/m_time.o

add_Input_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/add_Input.c -o $(MY_MISC)/add_Input.o
	
cmp_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/cmp.c -o $(MY_MISC)/cmp.o

matrix_c:
	gcc -I$(GSL_INCLUDE) -I$(MY_INCLUDE) -Werror -Wall -c $(MY_MISC)/matrix.c -o $(MY_MISC)/matrix.o
