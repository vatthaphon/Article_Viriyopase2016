/*
 * csln.h
 *
 *  Created on: May 31, 2012
 *      Author: viriyopa
 */

#ifndef CSLN_H_
#define CSLN_H_

/**************************** Directive declaration ****************************/
/* Parameters related to GSL */
#define EPSABS 		(1e-8)
#define EPSREL 		(1e-8)
#define GSL_TYPE  	gsl_odeiv2_step_rkf45
//#define GSL_TYPE  	gsl_odeiv2_step_rk2

/* Parameters related to Lyapunov exponents */
#define IS_DETERMINE_FULL_LEs 	(0)
//#define N_FIRST_LEs 			(2)
#define N_FIRST_LEs 			(1)
//#define PERT_SIZE_LEs 			(1e-6)
#define PERT_SIZE_LEs 			(0.5*1e-3)

#define IS_CONN_WITH_0E_0I		0
#define IS_CONN_WITH_0E_1I		1
#define IS_CONN_WITH_1E_0I		2
#define IS_CONN_WITH_1E_1I		3

/* Parameters related to special functions */
#define	M_EXP	exp
#define M_LOG	log

//#define	M_EXP	gsl_sf_exp
//#define M_LOG	gsl_sf_log

/* Type of interneurons */
#define WANG_BUZSAKI		(0)
#define BORGER_WALKER		(1)
//#define INTN_TYPE			(BORGER_WALKER)
#define INTN_TYPE			(WANG_BUZSAKI)

//#define WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI // To use this, we have also to define INTN_TYPE == (WANG_BUZSAKI)

/* Dump how the neurons are connected */
//#define IS_DUMP_NETWORK_EXITS	(1)
#define IS_DUMP_NETWORK_EXITS	(0)

/**************************** Constant declaration ****************************/
#define BIEXP_SYNAPSE		(0)
#define PULSE_SYNAPSE		(1)

#define IS_SELF_II_ALLOW 	1
#define IS_SELF_EE_ALLOW 	1
#define IS_SELF_II_NOTALLOW 0
#define IS_SELF_EE_NOTALLOW 0

#define	CELL_NOGEN_SPK_STATE 	(0)
#define	CELL_GEN_SPK_STATE 		(1)

#define PYR_CELL_CA1	0
#define PYR_CELL_CA3	1

#define TOO_BIG_NUMBER		(1e+14)
#define TOO_SMALL_NUMBER	(1e-14)

#define SELF_SPK 0
#define RECV_SPK 1

#define MINUS_INFINITY -1000000000.0

#define LIMIT_ITER 1
#define UNLIMIT_ITER 0

#define PI 3.14159265

#define NO_I_NEURONS -1
#define NO_E_NEURONS -2

/* Characteristic of poisson external inputs */
#define POISSON_EXTERNAL_INPUTS 			0
#define PERIODIC_EXTERNAL_INPUTS 			1
#define PERIODIC_POISSON_EXTERNAL_INPUTS 	2

/**************************** Struct declaration ****************************/
struct spiking
{
	struct spiking *pre;

	double	time;

	struct spiking *next;
};

struct post_syn_osc
{
	struct post_syn_osc 	*pre;

	int 				id;

	struct post_syn_osc 	*next;
};

struct syn_props
{
	double g_syn_AMPA_on_E;
	double g_syn_AMPA_on_I;
	double g_syn_GABA_on_E;
	double g_syn_GABA_on_I;

	double v_syn_AMPA;
	double v_syn_GABA;

	double tau_m_on_E;
	double tau_m_on_I;

	double tau_l_AMPA_on_E;
	double tau_l_AMPA_on_I;
	double tau_l_GABA_on_E;
	double tau_l_GABA_on_I;

	double tau_r_AMPA_on_E;
	double tau_r_AMPA_on_I;
	double tau_r_GABA_on_E;
	double tau_r_GABA_on_I;

	double tau_d_AMPA_on_E;
	double tau_d_AMPA_on_I;
	double tau_d_GABA_on_E;
	double tau_d_GABA_on_I;
};

struct E_neuron_props
{
	double V_rest;     	// Resting membrane potential [mV]
	double V_th;       	// The spike threshold [mV]
	double V_reset;    	// The reset potential [mV]
	double V_refrac;   	// Absolute refractory period [ms]
	double tau_m;      	// Membrane time constants [ms]
	double C_m;       	// Capacitance [pF]
	double R_m;
};

struct I_neuron_props
{
	/* Typical interneuron properties */
	double V_rest;      // Resting membrane potential [mV]
	double V_th;       	// The spike threshold [mV]
	double V_reset;     // The reset potential [mV]
	double V_refrac;    // Absolute refractory period [ms]
	double tau_m;       // Membrane time constants [ms]
	double C_m;       	// Capacitance [pF]
	double R_m;
};

struct network_props
{
	int	N_i;       	// >0,The number of interneurons
	int	N_e;        // >0,The number of pyramidal neurons
	int	N_ext;  	// >0,1, The number of external oscillators that generate the external inputs

	double p_EE;   	// 0.2,from E to E
	double p_EI;   	// 0.2,from E to I
	double p_IE;   	// 0.2,from I to E
	double p_II;   	// 0.2,from I to I

	double p_EE_gap_junction;
	double p_EI_gap_junction;
	double p_IE_gap_junction;
	double p_II_gap_junction;
};

struct factor_refrac_st
{
	int	is_not_in_refrac;
	struct factor_refrac_st *next;
};

struct ext_input_props
{
//	double lambda_on_E;     // firing rate. Lambda spikes per ms.
//	double lambda_on_I;     // firing rate. Lambda spikes per ms.
	double lambda_AMPA_on_E;    // firing rate. Lambda spikes per ms.
	double lambda_AMPA_on_I;    // firing rate. Lambda spikes per ms.
	double lambda_GABA_on_E;    // firing rate. Lambda spikes per ms.
	double lambda_GABA_on_I;    // firing rate. Lambda spikes per ms.


	double v_syn_AMPA;
	double v_syn_GABA;

	double g_ext_AMPA_on_E;	// [nS]
	double g_ext_AMPA_on_I; // [nS]

	double g_ext_GABA_on_E;	// [nS]
	double g_ext_GABA_on_I; // [nS]

	double tau_m_on_I;
	double tau_m_on_E;

	double tau_l_AMPA_on_E;
	double tau_l_AMPA_on_I;
	double tau_l_GABA_on_E;
	double tau_l_GABA_on_I;

	double tau_r_AMPA_on_E;
	double tau_r_AMPA_on_I;
	double tau_r_GABA_on_E;
	double tau_r_GABA_on_I;

	double tau_d_AMPA_on_E;
	double tau_d_AMPA_on_I;
	double tau_d_GABA_on_E;
	double tau_d_GABA_on_I;


	double RmIe_E, RmIe_I;
};

struct sim_props
{
	double 	dt;	// Step size. [dt]=mS
	int		Nt; // The number of rounds of simulation, i.e. tEnd=dt*Nt
	long	seed;
	int		isUsingSeed;
};

struct X_elem
{
	double 	x_E, x_I, x_E_ext, x_I_ext;

	struct X_elem *next;
};

struct S_elem
{
	double 	s_E, s_I, s_E_ext, s_I_ext;

	struct S_elem *next;
};

typedef struct post_syn_osc 	post_syn_osc_struct;

struct ode_solver_params
{
	struct E_neuron_props 	*e_neuron_props;
	struct I_neuron_props 	*i_neuron_props;

	struct syn_props		*syn_props;
	struct ext_input_props	*ext_input_props;

	struct network_props	*network_props;

	struct factor_refrac_st	*factor_refrac_ptr;
	struct factor_refrac_st	*factor_refrac_ptr_ptraj;

	struct S_elem			*head_S;
	struct X_elem			*head_X;

	struct S_elem			*head_S_ptraj;
	struct X_elem			*head_X_ptraj;

	int						LEs_i;
	int						N_LEs;
	double 					*v_perts_ptr;
	double 					*s_E_perts_ptr;
	double 					*x_E_perts_ptr;
	double 					*s_I_perts_ptr;
	double 					*x_I_perts_ptr;
	int						*is_connected_with_E;
	int						*is_connected_with_I;

	double 					*RmIe_arrays_ptr;

	int						is_initV_mode;

	post_syn_osc_struct		**pre_syn_Mat_ptr;
};

struct v_pert_elem
{
	struct v_pert_elem		*prev;

	double dv;
	int is_from_E;

	struct v_pert_elem		*next;
};

struct input
{
	struct input *pre;

	double 	arr_time;	// [ms]
	double 	N_input;	// The number of inputs that are activated by the external synapses
	int		from;		// Oscillator ID that generates this event

	struct v_pert_elem *v_pert;
	double 	dvdt;

	struct input *next;
};


typedef struct input 			input_struct;
typedef struct spiking 			spiking_struct;
typedef struct syn_props 		syn_props_struct;
typedef struct E_neuron_props 	E_neuron_props_struct;
typedef struct I_neuron_props 	I_neuron_props_struct;
typedef struct network_props 	network_props_struct;
typedef struct ext_input_props 	ext_input_props_struct;
typedef struct sim_props 		sim_props_struct;
typedef struct ode_solver_params ode_solver_params_struct;
typedef struct factor_refrac_st	factor_refrac_st_struct;
typedef struct X_elem			X_elem_struct;
typedef struct S_elem			S_elem_struct;
typedef struct v_pert_elem		v_pert_elem_struct;

/**************************** Function declaration ****************************/
double 	interpolate_CA1_t2spk(double I);
double 	interpolate_CA1_h_KDR(double I);
double 	interpolate_CA1_h_CaH(double I);
double 	interpolate_CA1_m_KDR(double I);
double 	interpolate_CA1_m_KM(double I);
double 	interpolate_CA1_h_CaT(double I);
double 	interpolate_CA1_h_NaT(double I);
double 	interpolate_CA1_m_CaH(double I);
double 	interpolate_CA1_m_CaT(double I);
double 	interpolate_CA1_V(double I);
double 	interpolate_CA1_T0(double I);

double 	interpolate_CA3_t2spk(double I);
double 	interpolate_CA3_h_KDR(double I);
double 	interpolate_CA3_h_CaH(double I);
double 	interpolate_CA3_m_KDR(double I);
double 	interpolate_CA3_m_KM(double I);
double 	interpolate_CA3_h_CaT(double I);
double 	interpolate_CA3_h_NaT(double I);
double 	interpolate_CA3_m_CaH(double I);
double 	interpolate_CA3_m_CaT(double I);
double 	interpolate_CA3_V(double I);
double 	interpolate_CA3_T0(double I);

double 	interpolate_WB_t2spk(double I);
double 	interpolate_WB_h(double I);
double 	interpolate_WB_n(double I);
double 	interpolate_WB_T0(double I);
double 	interpolate_WB_V(double I);

double 	interpolate_White_t2spk(double I);
double 	interpolate_White_h(double I);
double 	interpolate_White_n(double I);
double 	interpolate_White_T0(double I);
double 	interpolate_White_V(double I);

double 	interpolate_BW_h(double I);
double 	interpolate_BW_n(double I);
double 	interpolate_BW_T0(double I);
double 	interpolate_BW_V(double I);

int 	brunel_wang_2003_network(int argc, char **argv, int isTellDone);

void 	add_Input(v_pert_elem_struct *v_pert, double dvdt, double t, double N_input, input_struct *head[], input_struct *tail[], double latency, int from, int input_i, int is_all_delays_the_same);
void 	add_Input_v2(int N_inputs_X_bef_or_eq[], int N_inputs_X_aft[], int N_inputs_EI_bef_or_eq[], int N_inputs_EI_aft[], v_pert_elem_struct *v_pert, double dvdt, double t, double N_input, input_struct *head[], input_struct *tail[], double latency, int from, int i_LE);

double 	getAEle(double array[], int x, int y, int x_l, int y_l);
double 	getIntPtrEle(int *int_ptr, int x, int y, int x_l, int y_l);
double 	getDoublePtrEle(double *double_ptr, int x, int y, int x_l, int y_l);

void 	setAEle(double array[], double ele, int x, int y, int x_l, int y_l);
void 	setIntPtrEle(int *int_ptr, double ele, int x, int y, int x_l, int y_l);
void 	setDoublePtrEle(double *double_ptr, double ele, int x, int y, int x_l, int y_l);

int 	cmp(double x, double y, double tol_eq);

char 	*getCurrentTime();

#endif /* CSLN_H_ */
