/*
 * brunel_wang_2003_network.c
 *
 *  Created on: Feb 6, 2012
 *      Author: viriyopa
 *
 *	Based on the paper by Brunel and Wang 2003
 *	LIF: tau_m*dV/dt=V_rest-V-R_m*I_syn, derived from Eq. 5.13 in
 * Peter & Abbot book, where I_syn=g*s(t)*(V-V_syn).
 * Units: [tau_m]=ms,[V]=mV,[t]=ms,[V_rest]=mV,[R_m]=GOhm,[g]=nS. R_m is
 * obtained from R_m=tau_m/C_m where C_m is the capacitance with [C_m]=pF.
 * Then, we need to know C_m but it isn't given in the paper. Hence, we take the
 * values from Geisler et al.,2005: C_m=200 pF (interneuron) and C_m=250 pF
 * (pyramidal). However, we can directly calculate C_m from Brunel and Wang's paper given in my matlab 'Determine_resistance.m':
 * C_m=200 pF (interneuron) and C_m=500 pF (pyramidal).
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "csln.h"

/**/
double t_obs_beg = -1;
double t_obs_end = -1;
//double t_obs_beg = 35.0;
//double t_obs_end = 41.5;
//int osc_obs = 498;
int osc_obs = 19;
int is_first = 1;
/**/

/****************************** Patameters ******************************/
/******* ODE solver ********/
static int is_use_GSL_2_solve_ODE = 1;		// Value is determined by the calling program.
static int is_use_explicit_Euler = 1;		// Active when is_use_GSL_2_solve_ODE == 0

/******* Random numbers ********/
static int is_use_rnd_init_perts_seed = 1;	// Value is determined by the calling program.
static long rnd_init_perts_seed = 1;		// Value is determined by the calling program.

static int is_rand_V_init_X = 1;			// Value is determined by the calling program. If is_rand_V_init_X = 0, the cells start with the steady states.
static int is_use_rnd_init_V = 1;			// Activated when is_rand_V_init_X == 1.
static long rnd_init_V = 1;					// Value is determined by the calling program.

static int is_use_rnd_init_RmIe_seed = 1;	// Value is determined by the calling program.
static long rnd_init_RmIe_seed = 1;			// Value is determined by the calling program.

static int is_use_rnd_determineConnections = 1;	// Value is determined by the calling program.
static long rnd_determineConnections = 1;		// Value is determined by the calling program.

static int is_use_rnd_duringExecution = 1;	// This random number is used during execution.
static long rnd_duringExecution = 1;		// Value is determined by the calling program. This effects the White noise and the Poisson noise.

static int is_use_rnd_white_noise_E = 1;
static long rnd_white_noise_E = 1;
static double sigma_white_noise_E = 0.0;
static gsl_rng *r_white_noise_E = NULL;

static int is_use_rnd_white_noise_I = 1;
static long rnd_white_noise_I = 1;
static double sigma_white_noise_I = 0.0;
static gsl_rng *r_white_noise_I = NULL;

/******** Constant external currents ********/
static int is_hetero_Ie_E = 0;		// Value is determined by the calling program.
static int is_hetero_Ie_I = 0;		// Value is determined by the calling program.

static double	sigma_Ie_E = 0.0;	// the standard deviation of the normal distribution [pA], Value is determined by the calling program.
static double	sigma_Ie_I = 0.0;	// the standard deviation of the normal distribution [pA], Value is determined by the calling program.

/******** Mode of the external input ********/
//static int	type_of_external_inputs = POISSON_EXTERNAL_INPUTS;	// Will be override by the passing parameters.
static int	type_of_external_inputs = PERIODIC_EXTERNAL_INPUTS;	// Will be override by the passing parameters.

/* type_of_external_inputs = PERIODIC_EXTERNAL_INPUTS */
//static int 	time_shift							= 0;	// [ms]
static int 	time_shift							= 88.04 + 4.0;	// [ms]

/* type_of_external_inputs = PERIODIC_POISSON_EXTERNAL_INPUTS */
static double periodic_poisson_A				= 1;
static double periodic_poisson_f				= 50;	// [Hz]
static double periodic_poisson_time_shift 		= 0; 	// [ms]

/* Continuous cosine external inputs */
static int		is_Turn_On_Cont_Sin_Ext_I		= 0;
static int 		is_All_phi_the_same				= 1;
static double 	A_cos_ext_input 				= 0.05;
static double 	f_cos_ext_input 				= 0.04;

/******** Test slowly increase and decrease the external input ********/
static int 		is_slowly_inc_dec_ext_input_to_I_cells	= 0;
static int 		is_slowly_inc_dec_ext_input_to_E_cells	= 0;

static int		I0_I_mode						= 1;	// 1, increasing current mode. -1, decreasing current mode.
static int		I0_I_running_I					= 1;

static double	I0_I_begin						= 180.6;	// microA.
static double	I0_I_end						= 210.8;	// microA.
static int		I0_I_N							= 2;	// I0_I_N is like linspace(I0_I_begin, I0_I_end, I0_I_N) and linspace(I0_I_end, I0_I_begin, I0_I_N) in Matlab
static int		I0_I_at_which_I					= 1;	// '1' means 'I0_I_begin', 'I0_I_N' means 'I0_I_end'.
static double 	I0_I_t_begin					= 100;	// when t>I0_I_t_begin, slowly-inc. and dec. are turned on.
//static double	I0_I_t_int						= 100;	// An interval after which we increase or decrease the external current.
static double	I0_I_t_int						= 200;	// An interval after which we increase or decrease the external current.
//static double	I0_I_t_int						= 500;	// An interval after which we increase or decrease the external current.
//static double	I0_I_t_int						= 1000;	// An interval after which we increase or decrease the external current.
//static double	I0_I_t_int						= 2000;	// An interval after which we increase or decrease the external current.
//static double	I0_I_t_int						= 3000;	// An interval after which we increase or decrease the external current.

static int		I0_E_mode						= 1;	// 1, increasing current mode. -1, decreasing current mode.
static int		I0_E_running_I					= 1;
static double	I0_E_begin						= 1174.0;	// microA.
static double	I0_E_end						= 2168.0;	// microA.
static int		I0_E_N							= 20;	// I0_I_N is like linspace(I0_I_begin, I0_I_end, I0_I_N) and linspace(I0_I_end, I0_I_begin, I0_I_N) in Matlab
static int		I0_E_at_which_I					= 1;	// '1' means 'I0_I_begin', 'I0_I_N' means 'I0_I_end'.
static double 	I0_E_t_begin					= 1000;	// when t>I0_I_t_begin, slowly-inc. and dec. are turned on.
static double	I0_E_t_int						= 100;	// An interval after which we increase or decrease the external current.

/******** Type of synpase ********/
static int 	is_Cond_synapse 			= 1;	// If is_Cond_synapse = 1, we use (V-E) instead of (V_th-E).

static int	syn_type_E2E				= BIEXP_SYNAPSE;	// Type declaration is given in the header file.

static int	syn_type_E2I				= BIEXP_SYNAPSE;
//static int	syn_type_E2I				= PULSE_SYNAPSE;
static int	is_E2I_suprathreshold		= 1;				// Active when syn_type_E2I == PULSE_SYNAPSE;
static double	E2I_suprathreshold_new_V= 60.0;				// Active when syn_type_E2I == PULSE_SYNAPSE;
static int 	is_I_spike					= 0;

static int	syn_type_I2E				= BIEXP_SYNAPSE;
static int	syn_type_I2I				= BIEXP_SYNAPSE;
/******** Distribution of delays ********/
static int 	is_dist_EE_latency = 0;	// 0 by default
static int 	is_dist_EI_latency = 0;	// 0 by default
static int 	is_dist_IE_latency = 0;	// 0 by default
static int 	is_dist_II_latency = 0;	// 0 by default

static double	sigma_EE_latency = -1;
static double	sigma_EI_latency = -1;
static double	sigma_IE_latency = -1;
static double	sigma_II_latency = -1;

static double	mean_standard_EE_latency = 2.50;
static double	mean_standard_EI_latency = 1.30;
static double	mean_standard_IE_latency = 0.95;
static double	mean_standard_II_latency = 0.60;

#define	sigma_standard_latency_factor (100.0)
static double	sigma_standard_EE_latency = 0.3 + (0.3/100.0)*sigma_standard_latency_factor;
static double	sigma_standard_EI_latency = 0.3 + (0.3/100.0)*sigma_standard_latency_factor;
static double	sigma_standard_IE_latency = 0.1 + (0.1/100.0)*sigma_standard_latency_factor;
static double	sigma_standard_II_latency = 0.1 + (0.1/100.0)*sigma_standard_latency_factor;

static int	allow_EE_connections	= IS_SELF_EE_NOTALLOW;
static int 	allow_II_connections	= IS_SELF_II_NOTALLOW;

/******** Lyapunov exponents ********/
static int 	is_rescaling_pert					= 1;
static double beg_rescaling_pert_after			= -1.0;	// [ms]
static int 	is_using_linearized_perturbations 	= 0;

static int 	is_apply_the_actual_perturbation_on_pert_vectors = 0; // This option is active when is_turn_on_pert == 1.When it's 1, the perturbation is on the perturbation vectors.

static int 	is_show_end_LEs						= 0;
static int 	is_show_end_pert 					= 0;		// Show the LE exponents and their associate direction.

static int 	is_change_pert_directions 			= 0;
static double t_to_apply_the_actual_perturbation = 25; 	// [ms]
static int 	show_Every_LEs_Pert_each_second 	= 10; 		// [ms], record the perturbations for around each show_Every_LEs_Pert_each_second ms.
static double t_Begin_Rec_LEs_percent 			= 0.0;		// Percent of the simulation time that we used for calculating LE.
static double t_Begin_cal_avg_avg_LEs_percent 	= 20.0;		// Percent of the simulation time that we used for calculating avarage of LE, !!! t_Begin_Rec_LEs <= t_Begin_cal_avg_avg_LEs
static int 	is_exit_execution					= 0;

/******** Parameters of E-cells (Nowacki et al., 2011) ********/
#define 		nw_A_e 		(21590*1e-8)			// [cm^2], Because Nowacki gave the parameters in term of /cm^2, we need to multiply A_i to get rid of /cm^2.
#define			nw_N_states	(8) 					// The number of states to describe a neuron EXCLUDING V.

static double	nw_E_Na 		= 60.0;					// [mV]
static double	nw_E_Ca 		= 90.0;					// [mV]
static double	nw_E_K 			= -85.0;				// [mV]
static double	nw_E_L 			= -65.0;				// [mV]

static double	nw_C_m			= (1.00*1e+6)*nw_A_e;	// [pF]

static double	nw_g_NaT 		= (65.0*1e+6)*nw_A_e;	// [nS]
static double	nw_g_NaP 		= (0.1*1e+6)*nw_A_e;	// [nS]
static double	nw_g_CaH 		= (2.6*1e+6)*nw_A_e;	// [nS]
static double	nw_g_L 			= (0.02*1e+6)*nw_A_e;	// [nS]

// CA3
#define 		nw_CA3_steady_V  (-76.678601773367888);	// [mV]
#define 		nw_CA3_g_KDR 	 ((10.0*1e+6)*nw_A_e);	// [nS]
#define 		nw_CA3_g_KM 	 ((1.65*1e+6)*nw_A_e);	// [nS]
#define 		nw_CA3_g_CaT 	 ((0.74*1e+6)*nw_A_e);	// [nS]
#define 		nw_CA3_Iapp_th	  (175.0)				// [pA]. The minimum current which the cell oscillates.
#define 		nw_CA3_Iapp_Maxth (1657.0)				// [pA]. The maximum current which the cell oscillates.

// CA1
#define 		nw_CA1_steady_V  (-75.328738986887444);	// [mV]
#define 		nw_CA1_g_KDR 	 ((9.5*1e+6)*nw_A_e);	// [nS]
#define 		nw_CA1_g_KM 	 ((0.8*1e+6)*nw_A_e);	// [nS]
#define 		nw_CA1_g_CaT 	 ((0.6*1e+6)*nw_A_e);	// [nS]
#define 		nw_CA1_Iapp_th	  (95.0)				// [pA]. The minimum current which the cell oscillates.
#define 		nw_CA1_Iapp_Maxth (275.0)				// [pA]. The maximum current which the cell oscillates.

// For using in the equations
static int		type_Pyr_cells	= PYR_CELL_CA1;
static double	nw_steady_V		= nw_CA1_steady_V;		// [mV]
static double	nw_g_KDR 		= nw_CA1_g_KDR;			// [nS]
static double	nw_g_KM 		= nw_CA1_g_KM;			// [nS]
static double	nw_g_CaT 		= nw_CA1_g_CaT;			// [nS]

/******** Parameters of Interneuron (Borger and Walker, 2013) ********/
#define 		bw_A_i 		(18069*1e-8)			// [cm^2], Because Borger gave the parameters in term of /cm^2, we need to multiply A_i to get rid of /cm^2.
#define			bw_N_states	(2)						// The number of states to describe a neuron EXCLUDING V.
#define 		bw_Iapp_th	(1174.485)			// [pA]. The current above which the cell generates an action potential starting from the resting state.

static double	bw_E_Na		= 60.0;					// [mV]
static double	bw_E_K 		= -90.0;				// [mV]
static double	bw_E_L 		= -70.0;				// [mV]

static double 	bw_steady_V = -69.83;				// [mV]

static double	bw_C_m		= (1.0*1e+6)*bw_A_i;	// [pF]

static double 	bw_g_L 		= (0.5*1e+6)*bw_A_i; 	// [nS]
static double 	bw_g_Na 	= (112.0*1e+6)*bw_A_i;	// [nS]
static double 	bw_g_K 		= (224.0*1e+6)*bw_A_i;	// [nS]

/******** Parameters of Interneuron (Wang and Buzsaki, 1996) ********/
#define 		wb_A_i 		(18069*1e-8)			// [cm^2], Because Wang gave the parameters in term of /cm^2, we need to multiply A_i to get rid of /cm^2.
#define			wb_N_states	(2)						// The number of states to describe a neuron EXCLUDING V.
#define 		wb_Iapp_th	(29.0)					// [pA]. The current above which the cell generates an action potential starting from the resting state.

static double 	wb_E_L 		= -65.0; 				// [mV]
static double 	wb_E_Na 	= 55.0; 				// [mV]
static double 	wb_E_K 		= -90.0; 				// [mV]
static double	wb_Phi		= 5.0;

static double 	wb_steady_V = -64.017564909636803;	// [mV]

static double	wb_C_m		= (1.0*1e+6)*wb_A_i;	// [pF]

static double 	wb_g_L 		= (0.1*1e+6)*wb_A_i; 	// [nS]
static double 	wb_g_Na 	= (35.0*1e+6)*wb_A_i;	// [nS]
static double 	wb_g_K 		= (9.0*1e+6)*wb_A_i;	// [nS]
static double 	wb_g_gap_junction = 0.1;			// [nS], Check Traub et al., J of Neuros 2001.

/******** Parameters of Interneuron (White et al, 1998) ********/
#define 		w_Iapp_th	(-115.0)					// Need to look [pA]. The current above which the cell generates an action potential starting from the resting state.

static double 	w_E_L 		= -60.0; 				// [mV]
static double 	w_E_Na 		= 45.0; 				// [mV]
static double 	w_E_K 		= -80.0; 				// [mV]

static double 	w_steady_V 	= -63.63;	// Need to look [mV]

static double	w_C_m		= (1.0*1e+6)*wb_A_i;	// [pF]

static double 	w_g_L 		= (0.1*1e+6)*wb_A_i; 	// [nS]
static double 	w_g_Na 		= (30.0*1e+6)*wb_A_i;	// [nS]
static double 	w_g_K 		= (20.0*1e+6)*wb_A_i;	// [nS]
static double 	w_g_gap_junction = 0.1;				// [nS], Check Traub et al., J of Neuros 2001.

/****************************** Static global variables ******************************/
/* Initial observation time */
//static double OBS_TIMES_ICELLS	= 10.0;
//static double OBS_TIMES_ECELLS	= 9.15 - 1.3 - 2.5 - 4;

static double OBS_TIMES_ICELLS	= 4.0 + 7.0;
static double OBS_TIMES_ECELLS	= 4.0;

/* Typical pyramidal properties */
static double V_e_rest;      	// Resting membrane potential [mV]
static double V_e_th;       	// The spike threshold [mV]
static double V_e_reset;       	// The reset potential [mV]
static double t_e_refrac;      	// Absolute refractory period [ms]
static double tau_m_e;        	// Membrane time constants [ms]
static double C_m_e;       		// Capacitance [pF]
static double R_m_e;			// Resistance

/* Typical interneuron properties */
static double V_i_rest; 		// Resting membrane potential [mV]
static double V_i_th;       	// The spike threshold [mV]
static double V_i_reset;       	// The reset potential [mV]
static double t_i_refrac;      	// Absolute refractory period [ms]
static double tau_m_i;        	// Membrane time constants [ms]
static double C_m_i;       		// Capacitance [pF]
static double R_m_i;			// Resistance

/* Global vars */
static double 	global_t 			= 0.0;
static double 	t_before_call_gsl 	= 0.0;
static int		osc_E_ext_id, osc_I_ext_id;

/* Characteristic of poisson external inputs */
static int	n_periodic_AMPA_external_inputs_to_E = 0;
static int	n_periodic_AMPA_external_inputs_to_I = 0;
static int 	n_periodic_GABA_external_inputs_to_E = 0;
static int 	n_periodic_GABA_external_inputs_to_I = 0;

static int 	is_all_delays_the_same = 0;

/****************************** Static function declarations ******************************/
static void adapt_I0_Icells(double t, int N_i, int N_e, double RmIe_arrays[], I_neuron_props_struct	*I_neuron_p);
static void adapt_I0_Ecells(double t, int N_i, int N_e, double RmIe_arrays[], E_neuron_props_struct *E_neuron_p);
static double bw_f_h(double V, double h);
static double bw_f_n(double V, double n);
static double bw_get_I_L(double g_L, double v);
static double bw_get_I_K(double g_K, double v, double n);
static double bw_get_I_Na(double g_Na, double v, double h);
static double bw_get_I_L(double g_L, double v);
static double bw_get_I_K(double g_K, double v, double n);
static double bw_get_I_Na(double g_Na, double v, double h);
static double bw_m_inf_V(double V);
static double bw_h_inf_V(double V);
static double bw_n_inf_V(double V);
static double bw_alpha_m_V(double V);
static double bw_beta_m_V(double V);
static double bw_alpha_h_V(double V);
static double bw_beta_h_V(double V);
static double bw_alpha_n_V(double V);
static double bw_beta_n_V(double V);
static double cal_obs_times(int i_osc, int N_i, double RmIe_arrays[]);
static double nw_rhs_gating_variable(double V, double V_x, double k_x, double x, double tau_x);
static double nw_cal_x_inf(double V, double V_x, double k_x);
static int is_GJ_connected(int post, int pre, post_syn_osc_struct *pre_syn_Mat[]);
static int connect_Neurons_via_GJ(post_syn_osc_struct *pre_syn_Mat[], int X, int Y);
static double wb_get_I_gap(double V_post, const double V_dV[], post_syn_osc_struct *pre_syn_list);
post_syn_osc_struct	*getPre_syn_MatPtrEle(post_syn_osc_struct **pre_syn_Mat_ptr, int x, int y, int x_l, int y_l);
static void init_Cells(double v[], double m_wb_HH_th, double wb_Icells[], int state_wb_Icells[], double m_mp_HH_th, double mp_Ecells[], int state_mp_Ecells[], int N_i, int N_e, ode_solver_params_struct *ode_solver_params, double RmIe_arrays[]);
static int isHHSpk(int *state_HH, double pre_x, double x, double x_th, double x_go_down_th);
static double wb_get_n_inf(double v);
static double wb_get_h_inf(double v);
static double wb_get_m_inf(double v);
static double wb_get_f_h(double v, double h);
static double wb_get_f_n(double v, double n);
static double wb_get_I_Na(double g_Na, double v, double h);
static double wb_get_I_K(double g_K, double v, double n);
static double wb_get_I_L(double g_L, double v);
static void my_ode_solver(double t, double dt, double perts_gsl[], int N_i, int N_osc, ode_solver_params_struct *ode_solver_params);
static void explicitEuler(double t, double dt, double perts_gsl[], int N_i, int N_osc, ode_solver_params_struct *ode_solver_params);
static double get_Rm(int i, int N_i, double Rm_E, double Rm_I);
static void init_RmIe(int N_i, int N_e, double Rm_E, double Rm_I, double RmIe_E, double RmIe_I, double RmIe_arrays[]);
static double get_effective_latency(double mean_tau_l, double sigma_tau_l);
static void cal_Cv_Ce_Ci(int i_Osc, int N_i, double V, double s_E, double s_I, double s_E_ext, double s_I_ext,
		syn_props_struct *syn_props, ext_input_props_struct *ext_input_props,
		double *Cv, double *Ce, double *Ci);
static double update_dV(double dt, double Cv, double Ce, double Ci, double dV, double dS_E, double dX_E, double dS_I, double dX_I, double tau_r_E, double tau_d_E, double tau_r_I, double tau_d_I, int is_connected_with_E, int is_connected_with_I);
static void update_V_perts(
		int is_connected_with_E, int is_connected_with_I,
		double Cv, double Ce, double Ci,
		int is_not_in_refrac,
		double tau_m_E, double tau_r_E, double tau_d_E, double tau_m_I, double tau_r_I, double tau_d_I,
		double t, int N_LE, int N_osc, int i_Osc,
		int N_inputs_EI, int N_inputs_E, int N_inputs_I,
		double A_E, double A_I,
		double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[],
		input_struct *head_inputArrays_from_EI[],
		int N_inputs_EI_bef_or_eq[], int N_inputs_EI_aft[],
		int N_inputs_E_bef_or_eq[], int N_inputs_E_aft[],
		int N_inputs_I_bef_or_eq[], int N_inputs_I_aft[]);
static double update_dS(double dS_prev, double dX_prev, double dt, double tau_r, double tau_d, double A);
static double update_dX(double dX_prev, double dt, double tau_r, double A);
static double cal_S_kernel(double t, double tau_r, double tau_d);
static void openFiles(
		char *abs_file_name_Str, 	int is_rec_spkp, 		FILE **pFileSpkProfile,
		char *abs_file_name_Str1, 	int is_rec_traj, 		FILE **pFileSpkProfile1,
		char *abs_file_name_Str2, 	int is_rec_LEs, 		FILE **pFileSpkProfile2,
		char *abs_file_name_Str3, 	int is_rec_LEs_End_Pert,FILE **pFileSpkProfile3,
		char *abs_file_name_Str4, 	int is_rec_LEs_Pert, 	FILE **pFileSpkProfile4,
		char *abs_file_name_Str5, 	int is_rec_LEs_End, 	FILE **pFileSpkProfile5);

static void print_meta_files(char *self_FN);
static void init_SX_link_list(
		int N_osc,
		factor_refrac_st_struct **head_factor_refrac_ptr, factor_refrac_st_struct **tmp_factor_refrac_ptr,
		S_elem_struct **head_S, S_elem_struct **tmp_S,
		X_elem_struct **head_X, X_elem_struct **tmp_X);
static void inc_X_and_update_S_X_pert_after_spike_arrival(double V[], S_elem_struct *head_S, syn_props_struct *syn_props, ext_input_props_struct *ext_input_props, int is_connected_with_E[], int is_connected_with_I[], factor_refrac_st_struct *head_factor_refrac_ptr, int N_LE, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], double t, int N_i, int N_e, X_elem_struct *head_X, input_struct *head_inputArrays[], input_struct *tail_inputArrays[], double *RmIe_arrays);
static void update_S_X(double dt, int N_i, int N_osc, S_elem_struct *head_S, X_elem_struct *head_X, syn_props_struct *syn_props, ext_input_props_struct *ext_input_props);
static double cal_s_syn(double t, double t0, double tau_r, double tau_d, double s0, double x0);
static void get_synapse_kinetics(
		int osc_id, int N_i,
		double *g_syn_E, double *g_syn_I, double *g_syn_E_ext, double *g_syn_I_ext,
		double *V_syn_E, double *V_syn_I, double *V_syn_E_ext, double *V_syn_I_ext,
		double *tau_m_E, double *tau_m_I, double *tau_m_E_ext, double *tau_m_I_ext,
		double *tau_r_E, double *tau_r_I, double *tau_r_E_ext, double *tau_r_I_ext,
		double *tau_d_E, double *tau_d_I, double *tau_d_E_ext, double *tau_d_I_ext,
		syn_props_struct *syn_props, ext_input_props_struct *ext_input_props);
static int RHS_v_v_perts(double t, const double V_dV[], double f_V_dV[], void *params);
static double get_R_m(int id_osc, int N_i);
static double get_tau_m(int id_osc, int N_i);
static int chk_is_all_delays_the_same(double *tau_l_AMPA_on_E, double *tau_l_AMPA_on_I, double *tau_l_GABA_on_E, double *tau_l_GABA_on_I, double *tau_l_ext_AMPA_on_E, double *tau_l_ext_AMPA_on_I, double *tau_l_ext_GABA_on_E, double *tau_l_ext_GABA_on_I);
static double get_RmIe(int i, int N_i, double RmIe_E, double RmIe_I);
static double update_V_and_PertVectors(ode_solver_params_struct *ode_solver_params, double epsabs, double epsrel, int N_i, int N_osc, int N_LE, double dt, double V[], double wb_Icells[], double mp_Ecells[], double v_perts[]);
static int display(int show_Every_i, int show_Every, double cur_t, double end_t, int isShowTime);
static void get_cand_pert(double pert_size, int N_osc, double cand[]);
static int isLinearlyDependent(double x1[], double x2[], int N, double tol_eq);
static double cal_norm(double x[], int N);
static void norm_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int is_connected_with_E[], int is_connected_with_I[], double norm[]);
static void norm_perts(double rescaling_factor[], int N_LE, int N_osc, double pert_size, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int is_connected_with_E[], int is_connected_with_I[], double norm[]);
//static void dot_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], double is_connected_with_E[], double is_connected_with_I[], double dot[]);
//static void sum_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double dot[], double norm[], double sum_RHS[]);
//static void subtract_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double sum_RHS[]);
static int rec_LEs(double t, double t_Begin_Rec_LEs, int N_LE, int cnt_LEs, double LEs[], FILE *pFileSpkProfile2);
static int cal_LEs(int N_of_Samples_of_LEs, int N_LE, int N_osc, double dt, double t, double t_Begin_Rec_LEs, double pert_size, double LEs[], double avg_LEs[], double norm[], double avg_avg_LEs[], int *N_of_Samples_of_avg_avg_LEs, double t_Begin_cal_avg_avg_LEs);
static void GSR(int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int is_connected_with_E[], int is_connected_with_I[], double norm[]);
static void init_perts(
		/* Input */
		int N_LE,
		int N_osc,
		double pert_size,
		int is_connected_with_E[], int is_connected_with_I[],
		/* Output */
		double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int *N_states);
static void init_arrays(double A[],double x,int N);
static void init_arrays_int(int A[], int x, int N);
static void update_factor_refrac(factor_refrac_st_struct *factor_refrac_ptr, int is_just_after_refrac[], double factor_counter[], sim_props_struct* sim_p, network_props_struct *network_p);
static void init_struct_arrays(void *A[],int N);
static void addSpkArrivalArrays_Ext_inputs(double t, input_struct *head_inputArrays[], input_struct *tail_inputArrays[], network_props_struct *network_p, ext_input_props_struct *ext_input_props, sim_props_struct *sim_p);
static void addSpkArrivalArrays(int N_LE, double v_perts[], double dvdt, int pre, double t, input_struct *head_inputArrays[], input_struct *tail_inputArrays[], post_syn_osc_struct *post_syn_Mat[], syn_props_struct *syn_p, network_props_struct *network_p);
static void addSpikingInfo(int i, double t, spiking_struct *head_spikingArrays[], spiking_struct *tail_spikingArrays[], int N);
static void free_struct_arrays(void *head_Arrays[],void *tail_Arrays[],int N,int type);
static int	recResults(double data,FILE *pFile,int cnt,int isAppendSeparator);
static void determineConnections(int is_connected_with_E[], int is_connected_with_I[],
		post_syn_osc_struct *post_syn_Mat[], post_syn_osc_struct *pre_syn_Mat[],
		double p_EE,double p_EI,double p_IE,double p_II,
		double p_EE_gap_junction, double p_EI_gap_junction, double p_IE_gap_junction, double p_II_gap_junction,
		int N_i, int N_e, int isGJExactNumberOfConns);
static int rec_Params(double wb_V_pk, double wb_Icells_pk[], double mp_V_pk, double mp_Ecells_pk[], int N_e, int N_i, double mp_Ecells[], double wb_Icells[], double v[], FILE *pFileSpkProfile,double RmIe_E, double RmIe_I, syn_props_struct *syn_p,E_neuron_props_struct *E_neuron_p,I_neuron_props_struct *I_neuron_p,network_props_struct *network_p,ext_input_props_struct *ext_input_props,sim_props_struct *sim_p,int cnt);
static void max_min_connected_Osc(
		S_elem_struct *head_S,
		S_elem_struct **S_max_connected_I, double *g_syn_E_max_CI, double *g_syn_I_max_CI, double *g_syn_E_ext_max_CI, double *g_syn_I_ext_max_CI, double *V_syn_E_max_CI, double *V_syn_I_max_CI, double *V_syn_E_ext_max_CI, double *V_syn_I_ext_max_CI,
		S_elem_struct **S_min_connected_I, double *g_syn_E_min_CI, double *g_syn_I_min_CI, double *g_syn_E_ext_min_CI, double *g_syn_I_ext_min_CI, double *V_syn_E_min_CI, double *V_syn_I_min_CI, double *V_syn_E_ext_min_CI, double *V_syn_I_ext_min_CI,
		S_elem_struct **S_max_connected_E, double *g_syn_E_max_CE, double *g_syn_I_max_CE, double *g_syn_E_ext_max_CE, double *g_syn_I_ext_max_CE, double *V_syn_E_max_CE, double *V_syn_I_max_CE, double *V_syn_E_ext_max_CE, double *V_syn_I_ext_max_CE,
		S_elem_struct **S_min_connected_E, double *g_syn_E_min_CE, double *g_syn_I_min_CE, double *g_syn_E_ext_min_CE, double *g_syn_I_ext_min_CE, double *V_syn_E_min_CE, double *V_syn_I_min_CE, double *V_syn_E_ext_min_CE, double *V_syn_I_ext_min_CE,
		int *max_connected_I, int *min_connected_I, int *max_connected_E, int *min_connected_E,
		network_props_struct *network_p, post_syn_osc_struct *post_syn_Mat[],
		syn_props_struct *syn_props, ext_input_props_struct *ext_input_props);
static int rec_Traj(int N_osc, double V[],
		S_elem_struct *S_max_connected_I, double g_syn_E_max_CI, double g_syn_I_max_CI, double g_syn_E_ext_max_CI, double g_syn_I_ext_max_CI, double V_syn_E_max_CI, double V_syn_I_max_CI, double V_syn_E_ext_max_CI, double V_syn_I_ext_max_CI,
		S_elem_struct *S_min_connected_I, double g_syn_E_min_CI, double g_syn_I_min_CI, double g_syn_E_ext_min_CI, double g_syn_I_ext_min_CI, double V_syn_E_min_CI, double V_syn_I_min_CI, double V_syn_E_ext_min_CI, double V_syn_I_ext_min_CI,
		S_elem_struct *S_max_connected_E, double g_syn_E_max_CE, double g_syn_I_max_CE, double g_syn_E_ext_max_CE, double g_syn_I_ext_max_CE, double V_syn_E_max_CE, double V_syn_I_max_CE, double V_syn_E_ext_max_CE, double V_syn_I_ext_max_CE,
		S_elem_struct *S_min_connected_E, double g_syn_E_min_CE, double g_syn_I_min_CE, double g_syn_E_ext_min_CE, double g_syn_I_ext_min_CE, double V_syn_E_min_CE, double V_syn_I_min_CE, double V_syn_E_ext_min_CE, double V_syn_I_ext_min_CE,
		FILE *pFileSpkProfile1, int max_connected_I, int min_connected_I, int max_connected_E, int min_connected_E, int cnt_traj);
static void init_global_variables();

#if (IS_DUMP_NETWORK_EXITS == 1)
static void dump_network(FILE *pFile, int N_e, int N_i, double V[], post_syn_osc_struct *gj_Mat[], post_syn_osc_struct *post_syn_Mat[]);
#endif

int brunel_wang_2003_network(int argc, char **argv, int isTellDone)
{
	init_global_variables();
	/******************************* Parameters *******************************/
	char 	*end;
	char 	*self_FN=*argv++;

	char 	*job_msg=*argv++;							//printf("%s\n",job_msg);fflush(stdout);
	char 	*abs_file_name_Str	=strtok(job_msg, " ");	printf("[%s] %s\n", getCurrentTime(), abs_file_name_Str);fflush(stdout);
	char 	*abs_file_name_Str1	=strtok(NULL, " ");		//printf("%s\n",abs_file_name_Str1);fflush(stdout);
	char 	*abs_file_name_Str2	=strtok(NULL, " ");		//printf("%s\n",abs_file_name_Str2);fflush(stdout);
	char 	*abs_file_name_Str3	=strtok(NULL, " ");		//printf("%s\n",abs_file_name_Str3);fflush(stdout);
	char 	*abs_file_name_Str4	=strtok(NULL, " ");		//printf("%s\n",abs_file_name_Str4);fflush(stdout);
	char 	*abs_file_name_Str5	=strtok(NULL, " ");		//printf("%s\n",abs_file_name_Str5);fflush(stdout);

	/* Deterministically random, obsolete */
	long	seed		=atol(strtok(NULL, " "));		//printf("%d\n",seed);
	int		isUsingSeed	=atoi(strtok(NULL, " "));

	/* Network properties, 0,...,N_i-1=Inh.; N_i,...,N_i+N_e-1=Exc.; N_i+N_e=Ext. */
	int	N_i		=atoi(strtok(NULL, " "));       		//printf("%d\n",N_i); // >0,The number of interneurons
	int	N_e		=atoi(strtok(NULL, " "));       		//printf("%d\n",N_e);// >0,The number of pyramidal neurons
	int	N_ext	=atoi(strtok(NULL, " "));  				//printf("%d\n",N_ext);// >0,1, The number of external oscillators that generate the external inputs

	/* Probability for each pair groups of neurons */
	double p_EE	=strtod(strtok(NULL, " "), &end);   	// E to E
	double p_EI	=strtod(strtok(NULL, " "), &end);   	// E to I
	double p_IE	=strtod(strtok(NULL, " "), &end);   	// I to E
	double p_II	=strtod(strtok(NULL, " "), &end);   	// I to I

	/* Typical pyramidal properties */
	V_e_rest    =strtod(strtok(NULL, " "), &end);       // Resting membrane potential [mV]
	V_e_th      =strtod(strtok(NULL, " "), &end);       // The spike threshold [mV]
	V_e_reset   =strtod(strtok(NULL, " "), &end);       // The reset potential [mV]
	t_e_refrac  =strtod(strtok(NULL, " "), &end);       // Absolute refractory period [ms]
	tau_m_e     =strtod(strtok(NULL, " "), &end);       // Membrane time constants [ms]
	C_m_e       =strtod(strtok(NULL, " "), &end);       // Capacitance [pF]
	R_m_e		=tau_m_e/C_m_e;

	/* Typical interneuron properties */
	V_i_rest    =strtod(strtok(NULL, " "), &end); 		// Resting membrane potential [mV]
	V_i_th      =strtod(strtok(NULL, " "), &end);       // The spike threshold [mV]
	V_i_reset   =strtod(strtok(NULL, " "), &end);       // The reset potential [mV]
	t_i_refrac  =strtod(strtok(NULL, " "), &end);       // Absolute refractory period [ms]
	tau_m_i     =strtod(strtok(NULL, " "), &end);       // Membrane time constants [ms]
	C_m_i       =strtod(strtok(NULL, " "), &end);       // Capacitance [pF]
	R_m_i		=tau_m_i/C_m_i;

	/* Typical synaptic properties */
	double g_syn_AMPA_on_E=strtod(strtok(NULL, " "), &end);		// [nS]
	double g_syn_AMPA_on_I=strtod(strtok(NULL, " "), &end);    	// [nS]
	double g_syn_GABA_on_E=strtod(strtok(NULL, " "), &end);    	// [nS]
	double g_syn_GABA_on_I=strtod(strtok(NULL, " "), &end);    	// [nS]

	double v_syn_AMPA=strtod(strtok(NULL, " "), &end);          // [mV]
	double v_syn_GABA=strtod(strtok(NULL, " "), &end);         	// [mV]

	double tau_m_on_E=strtod(strtok(NULL, " "), &end);     		// [ms]
	double tau_m_on_I=strtod(strtok(NULL, " "), &end);     		// [ms]

	double tau_l_AMPA_on_E=strtod(strtok(NULL, " "), &end);     // [ms],Synaptic latency
	double tau_l_AMPA_on_I=strtod(strtok(NULL, " "), &end);     // [ms],Synaptic latency
	double tau_l_GABA_on_E=strtod(strtok(NULL, " "), &end);     // [ms],Synaptic latency
	double tau_l_GABA_on_I=strtod(strtok(NULL, " "), &end);     // [ms],Synaptic latency

	double tau_r_AMPA_on_E=strtod(strtok(NULL, " "), &end);     // [ms]
	double tau_r_AMPA_on_I=strtod(strtok(NULL, " "), &end);     // [ms]
	double tau_r_GABA_on_E=strtod(strtok(NULL, " "), &end);     // [ms]
	double tau_r_GABA_on_I=strtod(strtok(NULL, " "), &end);     // [ms]

	double tau_d_AMPA_on_E=strtod(strtok(NULL, " "), &end);     // [ms]
	double tau_d_AMPA_on_I=strtod(strtok(NULL, " "), &end);     // [ms]
	double tau_d_GABA_on_E=strtod(strtok(NULL, " "), &end);     // [ms]
	double tau_d_GABA_on_I=strtod(strtok(NULL, " "), &end);     // [ms]

	/* Typical external input properties */
	double lambda_AMPA_on_E=strtod(strtok(NULL, " "), &end);    // ++firing rate. Lambda spikes per ms.
	double lambda_AMPA_on_I=strtod(strtok(NULL, " "), &end);    // ++firing rate. Lambda spikes per ms.
	double lambda_GABA_on_E=strtod(strtok(NULL, " "), &end);    // ++firing rate. Lambda spikes per ms.
	double lambda_GABA_on_I=strtod(strtok(NULL, " "), &end);    // ++firing rate. Lambda spikes per ms.

	double v_ext_syn_AMPA=strtod(strtok(NULL, " "), &end);
	double v_ext_syn_GABA=strtod(strtok(NULL, " "), &end);

	double g_ext_AMPA_on_E=strtod(strtok(NULL, " "), &end);   	// [nS]
	double g_ext_AMPA_on_I=strtod(strtok(NULL, " "), &end);   	// [nS]
	double g_ext_GABA_on_E=strtod(strtok(NULL, " "), &end);   	// [nS]
	double g_ext_GABA_on_I=strtod(strtok(NULL, " "), &end);   	// [nS]

	double tau_m_ext_on_E=strtod(strtok(NULL, " "), &end);     	// [ms]
	double tau_m_ext_on_I=strtod(strtok(NULL, " "), &end);     	// [ms]

	double tau_l_ext_AMPA_on_E=strtod(strtok(NULL, " "), &end);
	double tau_l_ext_AMPA_on_I=strtod(strtok(NULL, " "), &end);
	double tau_l_ext_GABA_on_E=strtod(strtok(NULL, " "), &end);
	double tau_l_ext_GABA_on_I=strtod(strtok(NULL, " "), &end);

	double tau_r_ext_AMPA_on_E=strtod(strtok(NULL, " "), &end);
	double tau_r_ext_AMPA_on_I=strtod(strtok(NULL, " "), &end);
	double tau_r_ext_GABA_on_E=strtod(strtok(NULL, " "), &end);
	double tau_r_ext_GABA_on_I=strtod(strtok(NULL, " "), &end);

	double tau_d_ext_AMPA_on_E=strtod(strtok(NULL, " "), &end);
	double tau_d_ext_AMPA_on_I=strtod(strtok(NULL, " "), &end);
	double tau_d_ext_GABA_on_E=strtod(strtok(NULL, " "), &end);
	double tau_d_ext_GABA_on_I=strtod(strtok(NULL, " "), &end);

	double RmIe_E=strtod(strtok(NULL, " "), &end);						// [mV]
	double RmIe_I=strtod(strtok(NULL, " "), &end);						// printf("%1.16lf\n", RmIe_I);fflush(stdout); //[mV]

	/* Simulation parameters */
	double 	dt=strtod(strtok(NULL, " "), &end);     					// Step size. [mS]
	int		Nt=atoi(strtok(NULL, " "));    								// The number of rounds of simulation, i.e. tEnd=dt*Nt

	char *DONE=strtok(NULL, " ");

	int is_turn_on_pert = atoi(strtok(NULL, " "));
	double tmp_pert_size = strtod(strtok(NULL, " "), &end);

	double 	show_each_second = strtod(strtok(NULL, " "), &end);	// Display output every show_each_second second.
	int 	is_rec_spkp = atoi(strtok(NULL, " "));
	int 	is_rec_traj = atoi(strtok(NULL, " "));				// Trajectories of the neurons are recorded when is_rec_traj == 1.
	int 	is_rec_LEs = atoi(strtok(NULL, " "));
	int 	is_rec_LEs_End_Pert = atoi(strtok(NULL, " "));
	int		is_rec_LEs_Pert = atoi(strtok(NULL, " "));
	int		is_rec_LEs_End = atoi(strtok(NULL, " "));

	is_hetero_Ie_E = atoi(strtok(NULL, " "));
	is_hetero_Ie_I = atoi(strtok(NULL, " "));

	sigma_Ie_E = strtod(strtok(NULL, " "), &end);
	sigma_Ie_I = strtod(strtok(NULL, " "), &end);

	rnd_init_perts_seed = atol(strtok(NULL, " "));
	rnd_init_V = atol(strtok(NULL, " "));
	rnd_init_RmIe_seed = atol(strtok(NULL, " "));
	rnd_determineConnections = atol(strtok(NULL, " "));
	rnd_duringExecution = atol(strtok(NULL, " "));

	/* White noise */
	sigma_white_noise_E = strtod(strtok(NULL, " "), &end);
	sigma_white_noise_I = strtod(strtok(NULL, " "), &end);
	rnd_white_noise_E = atol(strtok(NULL, " "));
	rnd_white_noise_I = atol(strtok(NULL, " "));

	is_use_GSL_2_solve_ODE = atoi(strtok(NULL, " "));

	type_of_external_inputs = atoi(strtok(NULL, " "));

	double p_EE_gap_junction = strtod(strtok(NULL, " "), &end);
	double p_EI_gap_junction = strtod(strtok(NULL, " "), &end);
	double p_IE_gap_junction = strtod(strtok(NULL, " "), &end);
	double p_II_gap_junction = strtod(strtok(NULL, " "), &end);

	wb_g_gap_junction = strtod(strtok(NULL, " "), &end);

	is_rand_V_init_X = atoi(strtok(NULL, " "));

	wb_g_L = (strtod(strtok(NULL, " "), &end)*1e+6)*wb_A_i;

	type_Pyr_cells = atoi(strtok(NULL, " "));

	allow_EE_connections = atoi(strtok(NULL, " "));
	allow_II_connections = atoi(strtok(NULL, " "));

	int isGJExactNumberOfConns = atoi(strtok(NULL, " "));
	/******************************* Variables *******************************/
	// CA1 or CA3 pyramidal cells?
	switch(type_Pyr_cells)
	{
	    case PYR_CELL_CA1:
	    	nw_g_KDR 		= nw_CA1_g_KDR;
	    	nw_g_KM 		= nw_CA1_g_KM;
	    	nw_g_CaT 		= nw_CA1_g_CaT;
	    	nw_steady_V		= nw_CA1_steady_V;

	        break;
	    case PYR_CELL_CA3:
	    	nw_g_KDR 		= nw_CA3_g_KDR;
	    	nw_g_KM 		= nw_CA3_g_KM;
	    	nw_g_CaT 		= nw_CA3_g_CaT;
	    	nw_steady_V		= nw_CA3_steady_V;

	        break;
	}

	printf("[%s] is_turn_on_pert = %d, pert_size = %1.16lf\n", getCurrentTime(), is_turn_on_pert, tmp_pert_size);fflush(stdout);

	is_all_delays_the_same = chk_is_all_delays_the_same(&tau_l_AMPA_on_E, &tau_l_AMPA_on_I, &tau_l_GABA_on_E, &tau_l_GABA_on_I, &tau_l_ext_AMPA_on_E, &tau_l_ext_AMPA_on_I, &tau_l_ext_GABA_on_E, &tau_l_ext_GABA_on_I);

	int i;

	// Voltages of both pyramidal cells and interneurons
	double V[N_i + N_e];

	// gating variables of I-cells
#if (INTN_TYPE == WANG_BUZSAKI)
	double wb_Icells[wb_N_states*N_i];
#endif
#if (INTN_TYPE == BORGER_WALKER)
	double wb_Icells[bw_N_states*N_i];
#endif
	int	   state_wb_Icells[N_i];		// State of wb_Icells. CELL_NOGEN_SPK_STATE=not excited, CELL_GEN_SPK_STATE=excited
	double pre_m_inf[N_i];
	double m_wb_HH_th = 0.9, m_wb_go_down_HH_th = 0.4; // Changing these values has to match with interpolate_WB_X

	// gating variables of E-cells
	double mp_Ecells[nw_N_states*N_e];
	int	   state_mp_Ecells[N_e];		// State of mp_Ecells. CELL_NOGEN_SPK_STATE=not excited, CELL_GEN_SPK_STATE=excited
	double pre_m[N_e];
	double m_mp_HH_th = 0.9, m_mp_go_down_HH_th = 0.4; // Changing these values has to match with interpolate_CA1{3}_X

	double wb_V_pk;
#if (INTN_TYPE == WANG_BUZSAKI)
	double wb_Icells_pk[wb_N_states];
#endif
#if (INTN_TYPE == BORGER_WALKER)
	double wb_Icells_pk[bw_N_states];
#endif
	double mp_V_pk;
	double mp_Ecells_pk[nw_N_states];

	// Synapses
	S_elem_struct *head_S = NULL, *tmp_S = NULL;
	X_elem_struct *head_X = NULL, *tmp_X = NULL;

	post_syn_osc_struct	*post_syn_Mat[N_i + N_e];		// post_syn_Mat[i] contains a linked list of postsynaptic neurons affected by presynaptic neuron i.
	post_syn_osc_struct	*pre_syn_Mat[N_i + N_e];

	double 					factor_counter 			[N_i + N_e];
	factor_refrac_st_struct	*head_factor_refrac_ptr = NULL, *tmp_factor_refrac_ptr = NULL;
	int 					is_just_after_refrac 	[N_i + N_e];

	input_struct 	*head_inputArrays	[N_i + N_e];	// This is used for calculating syn. strength more quickly.
	input_struct 	*tail_inputArrays	[N_i + N_e];	// This is used for adding syn. strength more quickly.

	spiking_struct	*head_spikingArrays	[N_i + N_e];	// Keep the spiking information for each neuron.
	spiking_struct	*tail_spikingArrays	[N_i + N_e];	// Keep the spiking information for each neuron.

	syn_props_struct		*syn_p		=NULL;
	E_neuron_props_struct	*E_neuron_p	=NULL;
	I_neuron_props_struct	*I_neuron_p	=NULL;
	network_props_struct	*network_p	=NULL;
	ext_input_props_struct	*ext_input_props=NULL;
	sim_props_struct		*sim_p		=NULL;

	/* Record trajectory variables */
	double g_syn_E_max_CI, g_syn_I_max_CI, g_syn_E_ext_max_CI, g_syn_I_ext_max_CI;
	double g_syn_E_min_CI, g_syn_I_min_CI, g_syn_E_ext_min_CI, g_syn_I_ext_min_CI;
	double g_syn_E_max_CE, g_syn_I_max_CE, g_syn_E_ext_max_CE, g_syn_I_ext_max_CE;
	double g_syn_E_min_CE, g_syn_I_min_CE, g_syn_E_ext_min_CE, g_syn_I_ext_min_CE;
	double V_syn_E_max_CI, V_syn_I_max_CI, V_syn_E_ext_max_CI, V_syn_I_ext_max_CI;
	double V_syn_E_min_CI, V_syn_I_min_CI, V_syn_E_ext_min_CI, V_syn_I_ext_min_CI;
	double V_syn_E_max_CE, V_syn_I_max_CE, V_syn_E_ext_max_CE, V_syn_I_ext_max_CE;
	double V_syn_E_min_CE, V_syn_I_min_CE, V_syn_E_ext_min_CE, V_syn_I_ext_min_CE;
	int max_connected_I, min_connected_I, max_connected_E, min_connected_E;
	S_elem_struct *S_max_connected_I = NULL, *S_min_connected_I = NULL, *S_max_connected_E = NULL, *S_min_connected_E = NULL;

	FILE *pFileSpkProfile = NULL, *pFileSpkProfile1 = NULL, *pFileSpkProfile2 = NULL, *pFileSpkProfile3 = NULL, *pFileSpkProfile4 = NULL, *pFileSpkProfile5 = NULL;

	/* Perturbation variables */
	int N_LE;		// How many L.E. do we want to see?

	/* If the number of LEs is more than the number of states, we set the number of LEs to be the number of states. */
	if (IS_DETERMINE_FULL_LEs == 1)
	{
		N_LE = N_i + N_e;
	}
	else
	{
		N_LE = N_FIRST_LEs;

		if ((N_i + N_e) < N_LE)
		{
			N_LE = N_i + N_e;
		}
	}

	double 	pert_size = PERT_SIZE_LEs;

	double 	v_perts				[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of voltage.
	double 	s_E_perts			[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of S from E neurons.
	double 	x_E_perts			[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of X from E neurons.
	double 	s_I_perts			[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of S from I neurons.
	double 	x_I_perts			[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of X from I neurons.

	double 	v_perts_beg			[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of voltage.
	double 	s_E_perts_beg		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of S from E neurons.
	double 	x_E_perts_beg		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of X from E neurons.
	double 	s_I_perts_beg		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of S from I neurons.
	double 	x_I_perts_beg		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of X from I neurons.

	double 	v_perts_end			[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of voltage.
	double 	s_E_perts_end		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of S from E neurons.
	double 	x_E_perts_end		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of X from E neurons.
	double 	s_I_perts_end		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of S from I neurons.
	double 	x_I_perts_end		[N_LE * (N_i + N_e)];	// L.E. x a perturbation vector of X from I neurons.

	int		is_connected_with_E [N_i + N_e];			// Only from the E neurons. Not include the external inputs
	int		is_connected_with_I [N_i + N_e];			// Only from the I neurons. Not include the external inputs
	int		N_states_of_pert_vector;									// The number of states (V, s_E, x_E, s_I, x_I). The number is calculated only with the neccessary connections.

	double 	LEs[N_LE], avg_LEs[N_LE], avg_avg_LEs[N_LE];
	int 	N_of_Samples_of_LEs = 0;
	double 	t_Begin_Rec_LEs 			= Nt*dt*t_Begin_Rec_LEs_percent/100.0;		// [ms]
	double 	t_Begin_cal_avg_avg_LEs 	= Nt*dt*t_Begin_cal_avg_avg_LEs_percent/100.0;		// [ms], !!! t_Begin_Rec_LEs <= t_Begin_cal_avg_avg_LEs

	int 	N_of_Samples_of_avg_avg_LEs = 0;

	/* GSL variables */
	double epsabs = EPSABS, epsrel = EPSREL;

	/* External constant currents */
	double RmIe_arrays[N_i + N_e];

	/* Keep the first spike time */
	int 	is_first_spike[N_i + N_e];
	double	first_spike_time[N_i + N_e];

	/********************************************* Pre-execution processes *********************************************/
	/******************************* Initialization *******************************/
	syn_p		=(syn_props_struct *)			malloc(sizeof(syn_props_struct));
	E_neuron_p	=(E_neuron_props_struct *)		malloc(sizeof(E_neuron_props_struct));
	I_neuron_p	=(I_neuron_props_struct *)		malloc(sizeof(I_neuron_props_struct));
	network_p	=(network_props_struct *)		malloc(sizeof(network_props_struct));
	ext_input_props	=(ext_input_props_struct *)	malloc(sizeof(ext_input_props_struct));
	sim_p		=(sim_props_struct *)			malloc(sizeof(sim_props_struct));

	syn_p->g_syn_AMPA_on_E	=g_syn_AMPA_on_E;
	syn_p->g_syn_AMPA_on_I	=g_syn_AMPA_on_I;
	syn_p->g_syn_GABA_on_E	=g_syn_GABA_on_E;
	syn_p->g_syn_GABA_on_I	=g_syn_GABA_on_I;

	syn_p->v_syn_AMPA		=v_syn_AMPA;
	syn_p->v_syn_GABA		=v_syn_GABA;

	syn_p->tau_m_on_E		=tau_m_on_E;
	syn_p->tau_m_on_I		=tau_m_on_I;

	syn_p->tau_l_AMPA_on_E	=tau_l_AMPA_on_E;
	syn_p->tau_l_AMPA_on_I	=tau_l_AMPA_on_I;
	syn_p->tau_l_GABA_on_E	=tau_l_GABA_on_E;
	syn_p->tau_l_GABA_on_I	=tau_l_GABA_on_I;

	syn_p->tau_r_AMPA_on_E	=tau_r_AMPA_on_E;
	syn_p->tau_r_AMPA_on_I	=tau_r_AMPA_on_I;
	syn_p->tau_r_GABA_on_E	=tau_r_GABA_on_E;
	syn_p->tau_r_GABA_on_I	=tau_r_GABA_on_I;

	syn_p->tau_d_AMPA_on_E	=tau_d_AMPA_on_E;
	syn_p->tau_d_AMPA_on_I	=tau_d_AMPA_on_I;
	syn_p->tau_d_GABA_on_E	=tau_d_GABA_on_E;
	syn_p->tau_d_GABA_on_I	=tau_d_GABA_on_I;

	E_neuron_p->V_rest		=V_e_rest;
	E_neuron_p->V_th		=V_e_th;
	E_neuron_p->V_reset		=V_e_reset;
	E_neuron_p->V_refrac	=t_e_refrac;
	E_neuron_p->tau_m		=tau_m_e;
	E_neuron_p->C_m			=C_m_e;
	E_neuron_p->R_m			=R_m_e;

	I_neuron_p->V_rest		=V_i_rest;
	I_neuron_p->V_th		=V_i_th;
	I_neuron_p->V_reset		=V_i_reset;
	I_neuron_p->V_refrac	=t_i_refrac;
	I_neuron_p->tau_m		=tau_m_i;
	I_neuron_p->C_m			=C_m_i;
	I_neuron_p->R_m			=R_m_i;

	network_p->N_i			=N_i;
	network_p->N_e			=N_e;
	network_p->N_ext		=N_ext;

	network_p->p_EE			=p_EE;
	network_p->p_EI			=p_EI;
	network_p->p_IE			=p_IE;
	network_p->p_II			=p_II;

	network_p->p_EE_gap_junction =p_EE_gap_junction;
	network_p->p_EI_gap_junction =p_EI_gap_junction;
	network_p->p_IE_gap_junction =p_IE_gap_junction;
	network_p->p_II_gap_junction =p_II_gap_junction;

	ext_input_props->lambda_AMPA_on_E	=lambda_AMPA_on_E;
	ext_input_props->lambda_AMPA_on_I	=lambda_AMPA_on_I;
	ext_input_props->lambda_GABA_on_E	=lambda_GABA_on_E;
	ext_input_props->lambda_GABA_on_I	=lambda_GABA_on_I;

	ext_input_props->v_syn_AMPA			=v_ext_syn_AMPA;
	ext_input_props->v_syn_GABA			=v_ext_syn_GABA;

	ext_input_props->g_ext_AMPA_on_E	=g_ext_AMPA_on_E;
	ext_input_props->g_ext_AMPA_on_I	=g_ext_AMPA_on_I;
	ext_input_props->g_ext_GABA_on_E	=g_ext_GABA_on_E;
	ext_input_props->g_ext_GABA_on_I	=g_ext_GABA_on_I;

	ext_input_props->tau_m_on_E 		=tau_m_ext_on_E;
	ext_input_props->tau_m_on_I 		=tau_m_ext_on_I;

	ext_input_props->tau_l_AMPA_on_E	=tau_l_ext_AMPA_on_E;
	ext_input_props->tau_l_AMPA_on_I	=tau_l_ext_AMPA_on_I;
	ext_input_props->tau_l_GABA_on_E	=tau_l_ext_GABA_on_E;
	ext_input_props->tau_l_GABA_on_I	=tau_l_ext_GABA_on_I;

	ext_input_props->tau_r_AMPA_on_E	=tau_r_ext_AMPA_on_E;
	ext_input_props->tau_r_AMPA_on_I	=tau_r_ext_AMPA_on_I;
	ext_input_props->tau_r_GABA_on_E	=tau_r_ext_GABA_on_E;
	ext_input_props->tau_r_GABA_on_I	=tau_r_ext_GABA_on_I;

	ext_input_props->tau_d_AMPA_on_E	=tau_d_ext_AMPA_on_E;
	ext_input_props->tau_d_AMPA_on_I	=tau_d_ext_AMPA_on_I;
	ext_input_props->tau_d_GABA_on_E	=tau_d_ext_GABA_on_E;
	ext_input_props->tau_d_GABA_on_I	=tau_d_ext_GABA_on_I;

	ext_input_props->RmIe_E				=RmIe_E;
	ext_input_props->RmIe_I				=RmIe_I;

	sim_p->dt 							=dt;
	sim_p->Nt							=Nt;
	sim_p->seed							=seed;
	sim_p->isUsingSeed					=isUsingSeed;

	init_struct_arrays((void *)post_syn_Mat, N_i+N_e);
	init_struct_arrays((void *)pre_syn_Mat, N_i+N_e);

	/* Initialize storage to keep the first spike times */
	init_arrays(first_spike_time, -1.0, N_i + N_e);
	init_arrays_int(is_first_spike, 1, N_i + N_e);

	/* Determine connected nodes. Topology of the network can change if we move the function to somewhere else because of the random number. */
	init_arrays_int(is_connected_with_E, 0, N_i + N_e);
	init_arrays_int(is_connected_with_I, 0, N_i + N_e);

	/* Initialize the connections */
	printf("[%s] Initialize connections\n", getCurrentTime());fflush(stdout);
	determineConnections(is_connected_with_E, is_connected_with_I, post_syn_Mat, pre_syn_Mat, p_EE, p_EI, p_IE, p_II, p_EE_gap_junction, p_EI_gap_junction, p_IE_gap_junction, p_II_gap_junction, N_i, N_e, isGJExactNumberOfConns);

	init_struct_arrays((void *)head_inputArrays		, N_i + N_e);
	init_struct_arrays((void *)tail_inputArrays		, N_i + N_e);

	init_struct_arrays((void *)head_spikingArrays	, N_i + N_e);
	init_struct_arrays((void *)tail_spikingArrays	, N_i + N_e);

	init_RmIe(N_i, N_e, R_m_e, R_m_i, RmIe_E, RmIe_I, RmIe_arrays);

	/* Initialize state to detect spiking */
#if (INTN_TYPE == WANG_BUZSAKI)
	init_arrays(wb_Icells, 0.0, wb_N_states*N_i);
#endif
#if (INTN_TYPE == BORGER_WALKER)
	init_arrays(wb_Icells, 0.0, bw_N_states*N_i);
#endif
	init_arrays_int(state_wb_Icells, CELL_NOGEN_SPK_STATE, N_i);
	init_arrays(mp_Ecells, 0.0, nw_N_states*N_e);
	init_arrays_int(state_mp_Ecells, CELL_NOGEN_SPK_STATE, N_e);

	init_arrays(factor_counter, 0.0, N_i + N_e);

	init_SX_link_list(N_i + N_e, &head_factor_refrac_ptr, &tmp_factor_refrac_ptr, &head_S, &tmp_S, &head_X, &tmp_X);

	init_arrays_int(is_just_after_refrac, 0, N_i + N_e);

	init_arrays(v_perts_beg, 0.0, N_LE * (N_i + N_e));
	init_arrays(s_E_perts_beg, 0.0, N_LE * (N_i + N_e));
	init_arrays(x_E_perts_beg, 0.0, N_LE * (N_i + N_e));
	init_arrays(s_I_perts_beg, 0.0, N_LE * (N_i + N_e));
	init_arrays(x_I_perts_beg, 0.0, N_LE * (N_i + N_e));

	init_arrays(v_perts_end, 0.0, N_LE * (N_i + N_e));
	init_arrays(s_E_perts_end, 0.0, N_LE * (N_i + N_e));
	init_arrays(x_E_perts_end, 0.0, N_LE * (N_i + N_e));
	init_arrays(s_I_perts_end, 0.0, N_LE * (N_i + N_e));
	init_arrays(x_I_perts_end, 0.0, N_LE * (N_i + N_e));


	/* Perturbation variables */
	if ((is_apply_the_actual_perturbation_on_pert_vectors == 1) && (is_turn_on_pert == 1))
	{
		init_arrays(v_perts, 0.0, N_LE * (N_i + N_e));
		init_arrays(s_E_perts, 0.0, N_LE * (N_i + N_e));
		init_arrays(x_E_perts, 0.0, N_LE * (N_i + N_e));
		init_arrays(s_I_perts, 0.0, N_LE * (N_i + N_e));
		init_arrays(x_I_perts, 0.0, N_LE * (N_i + N_e));
	}
	else
	{
		init_perts(N_LE, N_i + N_e, pert_size, is_connected_with_E, is_connected_with_I, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, &N_states_of_pert_vector);
	}

	init_arrays(LEs, 0.0, N_LE);
	init_arrays(avg_LEs, 0.0, N_LE);
	init_arrays(avg_avg_LEs, 0.0, N_LE);

	sigma_EE_latency = (sigma_standard_EE_latency/mean_standard_EE_latency)*syn_p->tau_l_AMPA_on_E;
	sigma_EI_latency = (sigma_standard_EI_latency/mean_standard_EI_latency)*syn_p->tau_l_AMPA_on_I;
	sigma_IE_latency = (sigma_standard_IE_latency/mean_standard_IE_latency)*syn_p->tau_l_GABA_on_E;
	sigma_II_latency = (sigma_standard_II_latency/mean_standard_II_latency)*syn_p->tau_l_GABA_on_I;

	// White noise
	r_white_noise_E = gsl_rng_alloc(gsl_rng_taus);
	if (is_use_rnd_white_noise_E == 1)
	{
		gsl_rng_set(r_white_noise_E, rnd_white_noise_E);
	}

	r_white_noise_I = gsl_rng_alloc (gsl_rng_taus);
	if (is_use_rnd_white_noise_I == 1)
	{
		gsl_rng_set(r_white_noise_I, rnd_white_noise_I);
	}

	/******************************* Execution *******************************/
	/* Copy beginning perturbations */
	for (i = 0; i < N_i + N_e; i++)
	{
		int j;
		for (j = 0; j < N_LE; j++)
		{
			setAEle(v_perts_beg, getAEle(v_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
			setAEle(s_E_perts_beg, getAEle(s_E_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
			setAEle(x_E_perts_beg, getAEle(x_E_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
			setAEle(s_I_perts_beg, getAEle(s_I_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
			setAEle(x_I_perts_beg, getAEle(x_I_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
		}
	}

	/* Print versions of files*/
	print_meta_files(self_FN);

	/* Assign ID for the external inputs */
	osc_E_ext_id = N_i + N_e;
	osc_I_ext_id = N_i + N_e + 1;

	/* Initialize GSL's parameters */
	ode_solver_params_struct *ode_solver_params = (ode_solver_params_struct *)malloc(sizeof(ode_solver_params_struct));
	ode_solver_params->e_neuron_props 	= E_neuron_p;
	ode_solver_params->i_neuron_props 	= I_neuron_p;
	ode_solver_params->syn_props 		= syn_p;
	ode_solver_params->ext_input_props 	= ext_input_props;
	ode_solver_params->factor_refrac_ptr= head_factor_refrac_ptr;
	ode_solver_params->network_props 	= network_p;
	ode_solver_params->head_S			= head_S;
	ode_solver_params->head_X			= head_X;
	ode_solver_params->v_perts_ptr		= &v_perts[0];
	ode_solver_params->s_E_perts_ptr	= &s_E_perts[0];
	ode_solver_params->x_E_perts_ptr	= &x_E_perts[0];
	ode_solver_params->s_I_perts_ptr	= &s_I_perts[0];
	ode_solver_params->x_I_perts_ptr	= &x_I_perts[0];
	ode_solver_params->is_connected_with_E = &is_connected_with_E[0];
	ode_solver_params->is_connected_with_I = &is_connected_with_I[0];
	ode_solver_params->N_LEs			= N_LE;
	ode_solver_params->RmIe_arrays_ptr	= &RmIe_arrays[0];
	ode_solver_params->is_initV_mode 	= 0;
	ode_solver_params->pre_syn_Mat_ptr 	= &pre_syn_Mat[0];

	/* Initialize Cells */
	printf("[%s] Initialize V\n", getCurrentTime());fflush(stdout);

	init_Cells(V, m_wb_HH_th, wb_Icells, state_wb_Icells, m_mp_HH_th, mp_Ecells, state_mp_Ecells, N_i, N_e, ode_solver_params, RmIe_arrays);

	printf("Init V:%1.16lf,%1.16lf \n", V[0], V[1]);

	if (is_use_rnd_duringExecution == 1)
	{
		srand48(rnd_duringExecution);
	}

	/* Prepare recoding files */
	openFiles(abs_file_name_Str, 	is_rec_spkp, 		&pFileSpkProfile,
			abs_file_name_Str1, is_rec_traj, 		&pFileSpkProfile1,
			abs_file_name_Str2, is_rec_LEs, 		&pFileSpkProfile2,
			abs_file_name_Str3, is_rec_LEs_End_Pert,&pFileSpkProfile3,
			abs_file_name_Str4, is_rec_LEs_Pert, 	&pFileSpkProfile4,
			abs_file_name_Str5, is_rec_LEs_End,		&pFileSpkProfile5);

	/* Determine the maximum and minimum connected oscillators */
	max_min_connected_Osc(
			head_S,
			&S_max_connected_I, &g_syn_E_max_CI, &g_syn_I_max_CI, &g_syn_E_ext_max_CI, &g_syn_I_ext_max_CI, &V_syn_E_max_CI, &V_syn_I_max_CI, &V_syn_E_ext_max_CI, &V_syn_I_ext_max_CI,
			&S_min_connected_I, &g_syn_E_min_CI, &g_syn_I_min_CI, &g_syn_E_ext_min_CI, &g_syn_I_ext_min_CI, &V_syn_E_min_CI, &V_syn_I_min_CI, &V_syn_E_ext_min_CI, &V_syn_I_ext_min_CI,
			&S_max_connected_E, &g_syn_E_max_CE, &g_syn_I_max_CE, &g_syn_E_ext_max_CE, &g_syn_I_ext_max_CE, &V_syn_E_max_CE, &V_syn_I_max_CE, &V_syn_E_ext_max_CE, &V_syn_I_ext_max_CE,
			&S_min_connected_E, &g_syn_E_min_CE, &g_syn_I_min_CE, &g_syn_E_ext_min_CE, &g_syn_I_ext_min_CE, &V_syn_E_min_CE, &V_syn_I_min_CE, &V_syn_E_ext_min_CE, &V_syn_I_ext_min_CE,
			&max_connected_I, &min_connected_I, &max_connected_E, &min_connected_E,
			network_p, post_syn_Mat, syn_p, ext_input_props);

	int 	show_Every_i = 0, show_Every = floor(show_each_second/dt);	// floor(msec./dt)
	int		show_Every_LEs_Pert_i = 0, show_Every_LEs_Pert = floor(show_Every_LEs_Pert_each_second/dt);	// floor(msec./dt)
	int 	t_i;
	int 	cnt_traj = 0, cnt_LEs = 0, cnt_LEs_End = 0,cnt_LEs_End_Pert = 0, cnt_LEs_Pert = 0;

	/* Write meta data */
	if (is_rec_LEs == 1)
	{
		cnt_LEs = recResults((double)N_LE, pFileSpkProfile2, cnt_LEs, 0);
	}

	if (is_rec_LEs_End == 1)
	{
		cnt_LEs_End = recResults((double)N_LE, pFileSpkProfile5, cnt_LEs_End, 0);
	}

	if (is_rec_LEs_End_Pert == 1)
	{
		cnt_LEs_End_Pert = recResults((double)N_LE, pFileSpkProfile3, cnt_LEs_End_Pert, 0);
		cnt_LEs_End_Pert = recResults((double)(N_i + N_e), pFileSpkProfile3, cnt_LEs_End_Pert, 0);

		int j;
		for (j = 0; j < N_i + N_e; j++)
		{
			cnt_LEs_End_Pert = recResults((double)(is_connected_with_E[j]), pFileSpkProfile3, cnt_LEs_End_Pert, 0);
			cnt_LEs_End_Pert = recResults((double)(is_connected_with_I[j]), pFileSpkProfile3, cnt_LEs_End_Pert, 0);
		}
	}

	if (is_rec_LEs_Pert == 1)
	{
		cnt_LEs_Pert = recResults((double)N_LE, pFileSpkProfile4, cnt_LEs_Pert, 0);
		cnt_LEs_Pert = recResults((double)(N_i + N_e), pFileSpkProfile4, cnt_LEs_Pert, 0);

		int j;
		for (j = 0; j < N_i + N_e; j++)
		{
			cnt_LEs_Pert = recResults((double)(is_connected_with_E[j]), pFileSpkProfile4, cnt_LEs_Pert, 0);
			cnt_LEs_Pert = recResults((double)(is_connected_with_I[j]), pFileSpkProfile4, cnt_LEs_Pert, 0);
		}

	}

#if (IS_DUMP_NETWORK_EXITS == 1)
	printf("[%s] Dump network\n", getCurrentTime());fflush(stdout);
	dump_network(pFileSpkProfile, N_e, N_i, V, pre_syn_Mat, post_syn_Mat);

	return EXIT_SUCCESS;
#endif

	/********************************************* Execution processes *********************************************/
	printf("[%s] Begin\n", getCurrentTime());fflush(stdout);

	int is_Already_Pert = 0;

	global_t = 0;
	for (t_i = 0; t_i < Nt; t_i++)
	{
		/* Store the current time */
		t_before_call_gsl = global_t;

		/* Record perturbation vectors */
		if (is_rec_LEs_Pert == 1)
		{
			show_Every_LEs_Pert_i = display(show_Every_LEs_Pert_i, show_Every_LEs_Pert, global_t, Nt*dt, 0);

			if (show_Every_LEs_Pert_i == 0)
			{
				for (i = 0; i < N_LE; i++)
				{
					int j;
					for (j = 0; j < N_i + N_e; j++)
					{
						cnt_LEs_Pert = recResults(getAEle(v_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile4, cnt_LEs_Pert, 0);

						if (is_connected_with_E[j] == 1)
						{
							cnt_LEs_Pert = recResults(getAEle(s_E_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile4, cnt_LEs_Pert, 0);
							cnt_LEs_Pert = recResults(getAEle(x_E_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile4, cnt_LEs_Pert, 0);
						}

						if (is_connected_with_I[j] == 1)
						{
							cnt_LEs_Pert = recResults(getAEle(s_I_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile4, cnt_LEs_Pert, 0);
							cnt_LEs_Pert = recResults(getAEle(x_I_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile4, cnt_LEs_Pert, 0);
						}
					}
				}
			}
		}

		/* Display */
		show_Every_i = display(show_Every_i, show_Every, global_t, Nt*dt, 1);

		/* Shall we dynamically change the external input to the I-cells? */
		if (is_slowly_inc_dec_ext_input_to_I_cells == 1)
		{
			adapt_I0_Icells(global_t, N_i, N_e, RmIe_arrays, I_neuron_p);
		}

		/* Shall we dynamically change the external input to the E-cells? */
		if (is_slowly_inc_dec_ext_input_to_E_cells == 1)
		{
			adapt_I0_Ecells(global_t, N_i, N_e, RmIe_arrays, E_neuron_p);
		}

		/* Store m_inf of interneurons */
		for (i = 0; i < N_i; i++)
		{
			pre_m_inf[i] = wb_get_m_inf(V[i]);
		}

		/* Store m of E-cells */
		for (i = 0; i < N_e; i++)
		{
			pre_m[i] = nw_cal_x_inf(V[N_i + i], -37, 5);
		}

		/* Update the states V, S, and X as well as return the next time step */
		global_t = update_V_and_PertVectors(ode_solver_params, epsabs, epsrel, N_i, N_i + N_e, N_LE, dt, V, wb_Icells, mp_Ecells, v_perts);

		/* Absolute refractory. Because we want the out of refrac event happens at t+dt, we need to play this function after update_V_and_PertVectors. */
		update_factor_refrac(head_factor_refrac_ptr, is_just_after_refrac, factor_counter, sim_p, network_p);

		/* Update S and X of synpases */
		update_S_X(dt, N_i, N_i + N_e, head_S, head_X, syn_p, ext_input_props);

		/* Check if neurons fire or just pass the refractory period */
		tmp_factor_refrac_ptr = head_factor_refrac_ptr;
		tmp_S = head_S;
		tmp_X = head_X;
		int i_osc;
		for (i_osc = 0; i_osc < N_i + N_e; i_osc++)
		{
			if (i_osc < N_i)
			{
				int tmp_state = state_wb_Icells[i_osc];

				if (isHHSpk(&tmp_state, pre_m_inf[i_osc], wb_get_m_inf(V[i_osc]), m_wb_HH_th, m_wb_go_down_HH_th) == 1)
				{
					/* Record the first spiking time */
					if (is_first_spike[i_osc] == 1)
					{
						first_spike_time[i_osc] = global_t;
						is_first_spike[i_osc] = 0;
					}

					/* Add spiking information */
					addSpikingInfo(i_osc, global_t, head_spikingArrays, tail_spikingArrays, N_i + N_e);

					/* If there is a neuron fire, add outputs. */
					addSpkArrivalArrays(N_LE, v_perts, 0.0, i_osc, global_t, head_inputArrays, tail_inputArrays, post_syn_Mat, syn_p, network_p);

					wb_V_pk = V[0];
					int j;
#if (INTN_TYPE == WANG_BUZSAKI)
					for (j = 0; j < wb_N_states; j++)
#endif
#if (INTN_TYPE == BORGER_WALKER)
					for (j = 0; j < bw_N_states; j++)
#endif
					{
						wb_Icells_pk[j] = wb_Icells[j];
					}
				}

				state_wb_Icells[i_osc] = tmp_state;
			}
			else
			{
				int tmp_state = state_mp_Ecells[i_osc - N_i];

				if (isHHSpk(&tmp_state, pre_m[i_osc - N_i], nw_cal_x_inf(V[i_osc], -37, 5), m_mp_HH_th, m_mp_go_down_HH_th) == 1)
				{
					/* Record the first spiking time */
					if (is_first_spike[i_osc] == 1)
					{
						first_spike_time[i_osc] = global_t;
						is_first_spike[i_osc] = 0;
					}

					/* Add spiking information */
					addSpikingInfo(i_osc, global_t, head_spikingArrays, tail_spikingArrays, N_i + N_e);

					/* If there is a neuron fire, add outputs. */
					addSpkArrivalArrays(N_LE, v_perts, 0.0, i_osc, global_t, head_inputArrays, tail_inputArrays, post_syn_Mat, syn_p, network_p);

					mp_V_pk = V[N_i];
					int j;
					for (j = 0; j < nw_N_states; j++)
					{
						mp_Ecells_pk[j] = mp_Ecells[j];
					}
				}

				state_mp_Ecells[i_osc - N_i] = tmp_state;
			}

			tmp_factor_refrac_ptr = tmp_factor_refrac_ptr->next;
			tmp_S = tmp_S->next;
			tmp_X = tmp_X->next;
		}

		/* This function has to do after checking spiking because we can use 0-ms delay. */
		inc_X_and_update_S_X_pert_after_spike_arrival(V, head_S, syn_p, ext_input_props, is_connected_with_E, is_connected_with_I, head_factor_refrac_ptr, N_LE, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, global_t, N_i, N_e, head_X, head_inputArrays, tail_inputArrays, &RmIe_arrays[0]);

		/* Gram-Schmidt Reorthogonalization. Generate orthogonal vectors that are not necessary to be orthonormal */
		double norm[N_LE];
		GSR(N_LE, N_i + N_e, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, is_connected_with_E, is_connected_with_I, norm); // No change direction implemented yet!!!!

		/* Compute LEs */
		N_of_Samples_of_LEs = cal_LEs(N_of_Samples_of_LEs, N_LE, N_i + N_e, dt, global_t, t_Begin_Rec_LEs, pert_size, LEs, avg_LEs, norm, avg_avg_LEs, &N_of_Samples_of_avg_avg_LEs, t_Begin_cal_avg_avg_LEs);

		/* Normalize the perturbation vectors */
		if (is_rescaling_pert == 1)
		{
			double rescaling_factor[N_LE];

			norm_perts(rescaling_factor, N_LE, N_i + N_e, pert_size, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, is_connected_with_E, is_connected_with_I, norm);

			/* Rescaling the delay perturbations */
			int i;
			for (i = 0; i < N_i + N_e; i++)
			{
				input_struct *tmp_input	= head_inputArrays[i];
				while (tmp_input != NULL)
				{
					if (tmp_input->from < N_i + N_e)
					{
						v_pert_elem_struct *tmp_v_pert = tmp_input->v_pert;
						int j;
						for (j = 0; j < N_LE; j++)
						{
							tmp_v_pert->dv = tmp_v_pert->dv*rescaling_factor[j];
							tmp_v_pert = tmp_v_pert->next;
						}
					}
					tmp_input = tmp_input->next;
				}
			}
		}

		/* Apply the actual perturbation on Voltages */
		if (is_turn_on_pert == 1)
		{
			if ((global_t > t_to_apply_the_actual_perturbation) && (is_Already_Pert == 0))
			{
				is_Already_Pert = 1;
			}
		}

		/* Add the Poisson external inputs */
		addSpkArrivalArrays_Ext_inputs(global_t, head_inputArrays, tail_inputArrays, network_p, ext_input_props, sim_p);

		/* Record LEs */
		if (is_rec_LEs == 1)
		{
			cnt_LEs = rec_LEs(global_t, t_Begin_Rec_LEs, N_LE, cnt_LEs, LEs, pFileSpkProfile2);
		}

		/* Record trajectories */
		if (is_rec_traj == 1)
		{
			cnt_traj = rec_Traj(N_i + N_e, V,
					S_max_connected_I, g_syn_E_max_CI, g_syn_I_max_CI, g_syn_E_ext_max_CI, g_syn_I_ext_max_CI, V_syn_E_max_CI, V_syn_I_max_CI, V_syn_E_ext_max_CI, V_syn_I_ext_max_CI,
					S_min_connected_I, g_syn_E_min_CI, g_syn_I_min_CI, g_syn_E_ext_min_CI, g_syn_I_ext_min_CI, V_syn_E_min_CI, V_syn_I_min_CI, V_syn_E_ext_min_CI, V_syn_I_ext_min_CI,
					S_max_connected_E, g_syn_E_max_CE, g_syn_I_max_CE, g_syn_E_ext_max_CE, g_syn_I_ext_max_CE, V_syn_E_max_CE, V_syn_I_max_CE, V_syn_E_ext_max_CE, V_syn_I_ext_max_CE,
					S_min_connected_E, g_syn_E_min_CE, g_syn_I_min_CE, g_syn_E_ext_min_CE, g_syn_I_ext_min_CE, V_syn_E_min_CE, V_syn_I_min_CE, V_syn_E_ext_min_CE, V_syn_I_ext_min_CE,
					pFileSpkProfile1, max_connected_I, min_connected_I, max_connected_E, min_connected_E, cnt_traj);
		}

		/* Shall we exit the execution? */
		if (is_exit_execution == 1)
		{
			t_i = Nt + 1;
		}
	}

	printf("[%s] Finish loop\n", getCurrentTime());
	/********************************************* After execution processes *********************************************/
	if (is_rec_LEs_End == 1)
	{
		for (i = 0; i < N_LE; i++)
		{
			cnt_LEs_End = recResults(LEs[i], pFileSpkProfile5, cnt_LEs_End, 0);
			cnt_LEs_End = recResults(avg_LEs[i], pFileSpkProfile5, cnt_LEs_End, 0);
			cnt_LEs_End = recResults(avg_avg_LEs[i], pFileSpkProfile5, cnt_LEs_End, 0);
		}

		fclose(pFileSpkProfile5);
	}

	if (is_rescaling_pert == 1)
	{
		if (is_show_end_LEs == 1)
		{
			for (i = 0; i < N_LE; i++)
			{
				printf("LE %d = (%1.16lf, %1.16lf, %1.16lf)\n", i, LEs[i], avg_LEs[i], avg_avg_LEs[i]);
			}
		}

		/* Print out avg_LEs and associated directions */
		if (is_show_end_pert == 1)
		{
			printf("[%s] >> LEs and associated directions \n", getCurrentTime());

			for (i = 0; i < N_LE; i++)
			{
				printf("LE %d = (%1.16lf, %1.16lf, %1.16lf)\n", i, LEs[i], avg_LEs[i], avg_avg_LEs[i]);

				int j;
				for (j = 0; j < N_i + N_e; j++)
				{
					printf("v_%d, %1.16lf", j, getAEle(v_perts, i, j, N_LE, N_i + N_e));

					if (is_connected_with_E[j] == 1)
					{
						printf(" s_E_%d, %1.16lf", j, getAEle(s_E_perts, i, j, N_LE, N_i + N_e));
						printf(" x_E_%d, %1.16lf", j, getAEle(x_E_perts, i, j, N_LE, N_i + N_e));
					}

					if (is_connected_with_I[j] == 1)
					{
						printf(" s_I_%d, %1.16lf", j, getAEle(s_I_perts, i, j, N_LE, N_i + N_e));
						printf(" x_I_%d, %1.16lf", j, getAEle(x_I_perts, i, j, N_LE, N_i + N_e));
					}

					printf("\n");
				}
			}

			printf("\n");
		}
	}

	if ((is_rescaling_pert == 0) && (is_apply_the_actual_perturbation_on_pert_vectors == 0))
	{
		/* Copy ending perturbations */
		for (i = 0; i < N_i + N_e; i++)
		{
			int j;
			for (j = 0; j < N_LE; j++)
			{
				setAEle(v_perts_end, getAEle(v_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
				setAEle(s_E_perts_end, getAEle(s_E_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
				setAEle(x_E_perts_end, getAEle(x_E_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
				setAEle(s_I_perts_end, getAEle(s_I_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
				setAEle(x_I_perts_end, getAEle(x_I_perts, j, i, N_LE, N_i + N_e), j, i, N_LE, N_i + N_e);
			}
		}

		double norm_pert_beg[N_LE], norm_pert_end[N_LE];
		int i_LE;
		for (i_LE = 0; i_LE < N_LE; i_LE++)
		{
			norm_GSR(i_LE, N_LE, N_i + N_e, v_perts_end, s_E_perts_end, x_E_perts_end, s_I_perts_end, x_I_perts_end, is_connected_with_E, is_connected_with_I, norm_pert_end);
			norm_GSR(i_LE, N_LE, N_i + N_e, v_perts_beg, s_E_perts_beg, x_E_perts_beg, s_I_perts_beg, x_I_perts_beg, is_connected_with_E, is_connected_with_I, norm_pert_beg);

			if (is_show_end_LEs == 1)
			{
				printf("No rescaling LE %d: %1.16lf\n", i_LE, 1/(Nt*dt)*M_LOG(norm_pert_end[i_LE]/norm_pert_beg[i_LE]));
			}
		}
	}

	// Print first spike time
//	printf("------------------------------------------------------------------------------------------------\n");
//	for (i = 0; i < N_i + N_e; i++)
//	{
//		if (i < N_i)
//		{
//			printf("%1.16lf\n", first_spike_time[i]);
//		}
//	}


	if (is_rec_spkp == 1)
	{
		printf("[%s] >> Record spiking\n", getCurrentTime());

		/* Record simulation properties */
		int		cnt=0;

		cnt = rec_Params(wb_V_pk, wb_Icells_pk, mp_V_pk, mp_Ecells_pk, N_e, N_i, mp_Ecells, wb_Icells, V, pFileSpkProfile, RmIe_E, RmIe_I, syn_p, E_neuron_p, I_neuron_p, network_p, ext_input_props, sim_p, cnt);

		/* Record the results */
		for (i = 0; i < N_i + N_e; i++)
		{
			spiking_struct	*tmp_spiking = head_spikingArrays[i];
			while (tmp_spiking != NULL)
			{
				cnt = recResults(tmp_spiking->time, pFileSpkProfile, cnt, 0);		// Rec. simulation duration.

				tmp_spiking = tmp_spiking->next;
			}
			cnt = recResults(-1.0, pFileSpkProfile, cnt, 0);					// Rec. separator.
		}

		/* Close the files */
		fclose(pFileSpkProfile);
	}

	if (is_rec_traj == 1)
	{
		fclose(pFileSpkProfile1);
	}

	if (is_rec_LEs == 1)
	{
		fclose(pFileSpkProfile2);
	}

	if (is_rec_LEs_End_Pert == 1)
	{
		printf("[%s] Record perturbations\n", getCurrentTime());

		for (i = 0; i < N_LE; i++)
		{
			int j;
			for (j = 0; j < N_i + N_e; j++)
			{
				cnt_LEs_End_Pert = recResults(getAEle(v_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile3, cnt_LEs_End_Pert, 0);

				if (is_connected_with_E[j] == 1)
				{
					cnt_LEs_End_Pert = recResults(getAEle(s_E_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile3, cnt_LEs_End_Pert, 0);
					cnt_LEs_End_Pert = recResults(getAEle(x_E_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile3, cnt_LEs_End_Pert, 0);
				}

				if (is_connected_with_I[j] == 1)
				{
					cnt_LEs_End_Pert = recResults(getAEle(s_I_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile3, cnt_LEs_End_Pert, 0);
					cnt_LEs_End_Pert = recResults(getAEle(x_I_perts, i, j, N_LE, N_i + N_e), pFileSpkProfile3, cnt_LEs_End_Pert, 0);
				}
			}
		}

		fclose(pFileSpkProfile3);
	}

	if (is_rec_LEs_Pert == 1)
	{
		fclose(pFileSpkProfile4);
	}

	/* Clear memory */
	printf("[%s] >> Free memory\n", getCurrentTime());

	free(syn_p);
	free(E_neuron_p);
	free(I_neuron_p);
	free(network_p);
	free(ext_input_props);
	free(sim_p);

	gsl_rng_free(r_white_noise_E);
	gsl_rng_free(r_white_noise_I);

	free_struct_arrays((void *)head_inputArrays		,(void *)tail_inputArrays	,N_i+N_e,0);
	free_struct_arrays((void *)head_spikingArrays	,(void *)tail_spikingArrays	,N_i+N_e,1);
	free_struct_arrays((void *)post_syn_Mat			,(void *)NULL				,N_i+N_e,2);
	free_struct_arrays((void *)pre_syn_Mat			,(void *)NULL				,N_i+N_e,2);

	tmp_factor_refrac_ptr = head_factor_refrac_ptr;
	tmp_S = head_S;
	tmp_X = head_X;
	for (i = 0; i < N_i + N_e; i++)
	{
		factor_refrac_st_struct	*tmp = tmp_factor_refrac_ptr;
		tmp_factor_refrac_ptr = tmp_factor_refrac_ptr->next;

		S_elem_struct *tmp_tmp = tmp_S;
		tmp_S = tmp_S->next;

		X_elem_struct *tmp_tmp_tmp = tmp_X;
		tmp_X = tmp_X->next;

		free(tmp);
		free(tmp_tmp);
		free(tmp_tmp_tmp);
	}

	free(ode_solver_params);

	printf("[%s] End\n", getCurrentTime());fflush(stdout);

	if (isTellDone == 1)
	{
		/* Tell the script that I'm done */
		FILE *pFile = fopen(DONE, "wb");
		if (pFile == NULL)
		{
			printf("[main] cannot open [%s]\n", DONE);
		}
		else
		{
			fclose(pFile);
		}
	}

	return EXIT_SUCCESS;
}

#if (IS_DUMP_NETWORK_EXITS == 1)
static void dump_network(FILE *pFile, int N_e, int N_i, double V[], post_syn_osc_struct *gj_Mat[], post_syn_osc_struct *post_syn_Mat[])
{
	int N = N_e + N_i;
	int pre, post;
	int cnt = 0;

	cnt = recResults(N_e, pFile, cnt, 1);
	cnt = recResults(N_i, pFile, cnt, 1);

	cnt = recResults(1000000.0, pFile, cnt, 0);

	/* Dump the initial voltages */
	for (pre = 0; pre < N; pre++)
	{
		cnt = recResults(V[pre], pFile, cnt, 1);
	}

	cnt = recResults(1000000.0, pFile, cnt, 0);

	/* Dump the connections of the chemical synapses */
	for (pre = 0; pre < N; pre++)
	{
		post_syn_osc_struct *head_post_syn_osc = post_syn_Mat[pre];
		while (head_post_syn_osc != NULL)
		{
			post = head_post_syn_osc->id;

			cnt = recResults(post, pFile, cnt, 1);

			head_post_syn_osc = head_post_syn_osc->next;
		}

		cnt = recResults(1000000.0, pFile, cnt, 0);
	}

	/* Dump the connections of the electrical synapses */
	for (pre = 0; pre < N; pre++)
	{
		post_syn_osc_struct *head_gj_osc = gj_Mat[pre];
		while (head_gj_osc != NULL)
		{
			post = head_gj_osc->id;

			cnt = recResults(post, pFile, cnt, 1);

			head_gj_osc = head_gj_osc->next;
		}

		cnt = recResults(1000000.0, pFile, cnt, 0);
	}
}
#endif

static void init_global_variables()
{
	global_t = 0.0;
	t_before_call_gsl = 0.0;
	osc_E_ext_id = 0;
	osc_I_ext_id = 0;
	n_periodic_AMPA_external_inputs_to_E = 0;
	n_periodic_AMPA_external_inputs_to_I = 0;
	n_periodic_GABA_external_inputs_to_E = 0;
	n_periodic_GABA_external_inputs_to_I = 0;
	is_all_delays_the_same = 0;
}

static double update_V_and_PertVectors(ode_solver_params_struct *ode_solver_params, double epsabs, double epsrel, int N_i, int N_osc, int N_LE, double dt,
		double V[], double wb_Icells[], double mp_Ecells[],
		double v_perts[])
{
	int N_e = N_osc - N_i;
	int i_LE, i_osc, total_N_states;
	double prev_V[N_osc];

#if (INTN_TYPE == WANG_BUZSAKI)
	int intn_N_states = wb_N_states;
#endif
#if (INTN_TYPE == BORGER_WALKER)
	int intn_N_states = bw_N_states;
#endif

	total_N_states = 2*N_osc + intn_N_states*N_i + nw_N_states*N_e;

	double perts_gsl[total_N_states];
	double prev_mp_Ecells[nw_N_states*N_e];

#if (INTN_TYPE == WANG_BUZSAKI)
	double prev_wb_Icells[wb_N_states*N_i];
#endif
#if (INTN_TYPE == BORGER_WALKER)
	double prev_wb_Icells[bw_N_states*N_i];
#endif

	/* Store previous voltages */
	for (i_osc = 0; i_osc < N_osc; i_osc++)
	{
		prev_V[i_osc] = V[i_osc];
	}

	/* Store additional interneuron states excluding V */
	for (i_osc = 0; i_osc < N_i; i_osc++)
	{
		int j;
#if (INTN_TYPE == WANG_BUZSAKI)
		for (j = 0; j < wb_N_states; j++)
		{
			prev_wb_Icells[wb_N_states*i_osc + j] = wb_Icells[wb_N_states*i_osc + j];
		}
#endif
#if (INTN_TYPE == BORGER_WALKER)
		for (j = 0; j < bw_N_states; j++)
		{
			prev_wb_Icells[bw_N_states*i_osc + j] = wb_Icells[bw_N_states*i_osc + j];
		}
#endif
	}

	/* Store additional E-cells states excluding V */
	for (i_osc = 0; i_osc < N_e; i_osc++)
	{
		int j;
		for (j = 0; j < nw_N_states; j++)
		{
			prev_mp_Ecells[nw_N_states*i_osc + j] = mp_Ecells[nw_N_states*i_osc + j];
		}
	}

	/* Evolution of the perturbed trajectory */
	for (i_LE = 0; i_LE < N_LE; i_LE++)
	{
		double tmp_t = t_before_call_gsl;

		/* Initialize perturbed voltages */
		for (i_osc = 0; i_osc < N_osc; i_osc++)
		{
			perts_gsl[2*i_osc] 		= prev_V[i_osc];
			perts_gsl[2*i_osc + 1] 	= getAEle(v_perts, i_LE, i_osc, N_LE, N_osc);
		}

		/* Initialize additional states for interneurons */
		for (i_osc = 0; i_osc < N_i; i_osc++)
		{
			int j;
#if (INTN_TYPE == WANG_BUZSAKI)
			for (j = 0; j < wb_N_states; j++)
			{
				perts_gsl[2*N_osc + wb_N_states*i_osc + j] = prev_wb_Icells[wb_N_states*i_osc + j];
			}
#endif
#if (INTN_TYPE == BORGER_WALKER)
			for (j = 0; j < bw_N_states; j++)
			{
				perts_gsl[2*N_osc + bw_N_states*i_osc + j] = prev_wb_Icells[bw_N_states*i_osc + j];
			}
#endif
		}

		/* Initialize additional states for E-cells */
		for (i_osc = 0; i_osc < N_e; i_osc++)
		{
			int j;
			for (j = 0; j < nw_N_states; j++)
			{
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*i_osc + j] = prev_mp_Ecells[nw_N_states*i_osc + j];
			}
		}

		ode_solver_params->LEs_i 		= i_LE;

		if (is_use_GSL_2_solve_ODE == 1)
		{
			/* Initialize GSL ode solver */
			gsl_odeiv2_system 	tmp_sys_gsl	= {RHS_v_v_perts, NULL, total_N_states, ode_solver_params};
			gsl_odeiv2_driver 	*tmp_d_gsl 	= gsl_odeiv2_driver_alloc_y_new(&tmp_sys_gsl, GSL_TYPE, dt, epsabs, epsrel);

			/* Evolve one step further */
			int s_gsl = gsl_odeiv2_driver_apply(tmp_d_gsl, &tmp_t, tmp_t + dt, perts_gsl); gsl_odeiv2_driver_free(tmp_d_gsl);
			if (s_gsl != GSL_SUCCESS)
			{
				printf ("error: driver returned %d and exit(0)\n", s_gsl);
				exit(0);
			}
		}
		else
		{
			my_ode_solver(tmp_t, dt, perts_gsl, N_i, N_osc, ode_solver_params);
			tmp_t = tmp_t + dt;
		}

		/* Evolve perturbed S */
		for (i_osc = 0; i_osc < N_osc; i_osc++)
		{
			int N_i = ode_solver_params->network_props->N_i;
			syn_props_struct 		*syn_props 			= ode_solver_params->syn_props;
			ext_input_props_struct 	*ext_input_props 	= ode_solver_params->ext_input_props;
			int 					*is_connected_with_E= ode_solver_params->is_connected_with_E;
			int 					*is_connected_with_I= ode_solver_params->is_connected_with_I;
			double 					*s_E_perts			= ode_solver_params->s_E_perts_ptr;
			double 					*x_E_perts			= ode_solver_params->x_E_perts_ptr;
			double 					*s_I_perts			= ode_solver_params->s_I_perts_ptr;
			double 					*x_I_perts			= ode_solver_params->x_I_perts_ptr;

			double dummy;
			double tau_r_E, tau_r_I;
			double tau_d_E, tau_d_I;
			get_synapse_kinetics(i_osc, N_i,
					&dummy, &dummy, &dummy, &dummy,
					&dummy, &dummy, &dummy, &dummy,
					&dummy, &dummy, &dummy, &dummy,
					&tau_r_E, &tau_r_I, &dummy, &dummy,
					&tau_d_E, &tau_d_I, &dummy, &dummy,
					syn_props, ext_input_props);

			if (getIntPtrEle(is_connected_with_E, 0, i_osc, 1, N_osc) == 1)
			{
				double ds_E_0 = getDoublePtrEle(s_E_perts, i_LE, i_osc, N_LE, N_osc);
				double dx_E_0 = getDoublePtrEle(x_E_perts, i_LE, i_osc, N_LE, N_osc);
				double tmp = cal_s_syn(dt, 0, tau_r_E, 	tau_d_E, ds_E_0, dx_E_0);

				setDoublePtrEle(s_E_perts, tmp, i_LE, i_osc, N_LE, N_osc);
			}

			if (getIntPtrEle(is_connected_with_I, 0, i_osc, 1, N_osc) == 1)
			{
				double ds_I_0 = getDoublePtrEle(s_I_perts, i_LE, i_osc, N_LE, N_osc);
				double dx_I_0 = getDoublePtrEle(x_I_perts, i_LE, i_osc, N_LE, N_osc);
				double tmp = cal_s_syn(dt, 0, tau_r_I, 	tau_d_I, ds_I_0, dx_I_0);

				setDoublePtrEle(s_I_perts, tmp, i_LE, i_osc, N_LE, N_osc);
			}
		}

		/* Evolve perturbed X */
		for (i_osc = 0; i_osc < N_osc; i_osc++)
		{
			int N_i = ode_solver_params->network_props->N_i;
			syn_props_struct 		*syn_props 			= ode_solver_params->syn_props;
			ext_input_props_struct 	*ext_input_props 	= ode_solver_params->ext_input_props;
			int 					*is_connected_with_E= ode_solver_params->is_connected_with_E;
			int 					*is_connected_with_I= ode_solver_params->is_connected_with_I;
			double 					*x_E_perts			= ode_solver_params->x_E_perts_ptr;
			double 					*x_I_perts			= ode_solver_params->x_I_perts_ptr;

			double tau_r_E, tau_r_I;
			double dummy;
			get_synapse_kinetics(i_osc, N_i,
					&dummy, &dummy, &dummy, &dummy,
					&dummy, &dummy, &dummy, &dummy,
					&dummy, &dummy, &dummy, &dummy,
					&tau_r_E, &tau_r_I, &dummy, &dummy,
					&dummy, &dummy, &dummy, &dummy,
					syn_props, ext_input_props);

			if (getIntPtrEle(is_connected_with_E, 0, i_osc, 1, N_osc) == 1)
			{
				setDoublePtrEle(x_E_perts, M_EXP(-dt/tau_r_E)*getDoublePtrEle(x_E_perts, i_LE, i_osc, N_LE, N_osc), i_LE, i_osc, N_LE, N_osc);
			}

			if (getIntPtrEle(is_connected_with_I, 0, i_osc, 1, N_osc) == 1)
			{
				setDoublePtrEle(x_I_perts, M_EXP(-dt/tau_r_I)*getDoublePtrEle(x_I_perts, i_LE, i_osc, N_LE, N_osc), i_LE, i_osc, N_LE, N_osc);
			}
		}

		/* Evolve perturbed voltages */
		for (i_osc = 0; i_osc < N_osc; i_osc++)
		{
			setAEle(v_perts, perts_gsl[2*i_osc + 1], i_LE, i_osc, N_LE, N_osc);
		}
	}

	/* Update the voltages */
	for (i_osc = 0; i_osc < N_osc; i_osc++)
	{
		V[i_osc] = perts_gsl[2*i_osc];
	}

	/* Update additional states for interneurons */
	for (i_osc = 0; i_osc < N_i; i_osc++)
	{
		int j;
#if (INTN_TYPE == WANG_BUZSAKI)
		for (j = 0; j < wb_N_states; j++)
		{
			wb_Icells[wb_N_states*i_osc + j] = perts_gsl[2*N_osc + wb_N_states*i_osc + j];
		}
#endif
#if (INTN_TYPE == BORGER_WALKER)
		for (j = 0; j < bw_N_states; j++)
		{
			wb_Icells[bw_N_states*i_osc + j] = perts_gsl[2*N_osc + bw_N_states*i_osc + j];
		}
#endif
	}

	/* Update additional states for E-cells */
	for (i_osc = 0; i_osc < N_e; i_osc++)
	{
		int j;
		for (j = 0; j < nw_N_states; j++)
		{
			mp_Ecells[nw_N_states*i_osc + j] = perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*i_osc + j];
		}
	}

	return t_before_call_gsl + dt;
}


static int display(int show_Every_i, int show_Every, double cur_t, double end_t, int isShowTime)
{
	if (show_Every_i > show_Every)
	{
		if (isShowTime == 1)
		{
			printf("[%s] %1.0lf/%1.0lf\n", getCurrentTime(), cur_t, end_t);fflush(stdout);
		}

		show_Every_i=0;
	}
	else
	{
		show_Every_i++;
	}

	return show_Every_i;
}

static int isLinearlyDependent(double x1[], double x2[], int N, double tol_eq)
{
	double x_dot = 0.0;
	int i;

	for (i = 0; i < N; i++)
	{
		x_dot += x1[i]*x2[i];
	}

	double norm_x1 = cal_norm(x1, N);
	double norm_x2 = cal_norm(x2, N);

	if (cmp(fabs(x_dot/(norm_x1*norm_x2)), 1.0, tol_eq) == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

static double cal_norm(double x[], int N)
{
	double tmp = 0.0;
	int i;

	for (i = 0; i < N; i++)
	{
		tmp += x[i]*x[i];
	}

	return sqrt(tmp);
}


static void norm_perts(double rescaling_factor[], int N_LE, int N_osc, double pert_size, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int is_connected_with_E[], int is_connected_with_I[], double norm[])
{
	int i;
	for (i = 0; i < N_LE; i++)
	{
		int j;
		double factor;
		if (beg_rescaling_pert_after <= global_t)
		{
			factor = pert_size/norm[i];
		}
		else
		{
			factor = 1.0;
		}

		rescaling_factor[i] = factor;

		/* Resize the perturbation */
		for (j = 0; j < N_osc; j++)
		{
			setAEle(v_perts, factor*getAEle(v_perts, i, j, N_LE, N_osc), i, j, N_LE, N_osc);

			if (is_connected_with_E[j] == 1)
			{
				setAEle(s_E_perts, factor*getAEle(s_E_perts, i, j, N_LE, N_osc), i, j, N_LE, N_osc);
				setAEle(x_E_perts, factor*getAEle(x_E_perts, i, j, N_LE, N_osc), i, j, N_LE, N_osc);
			}

			if (is_connected_with_I[j] == 1)
			{
				setAEle(s_I_perts, factor*getAEle(s_I_perts, i, j, N_LE, N_osc), i, j, N_LE, N_osc);
				setAEle(x_I_perts, factor*getAEle(x_I_perts, i, j, N_LE, N_osc), i, j, N_LE, N_osc);
			}
		}
	}
}

static int rec_LEs(double t, double t_Begin_Rec_LEs, int N_LE, int cnt_LEs, double LEs[], FILE *pFileSpkProfile2)
{
	if (t_Begin_Rec_LEs <= t)
	{
		int i;
		for (i = 0; i < N_LE; i++)
		{
			cnt_LEs = recResults(LEs[i], pFileSpkProfile2, cnt_LEs, 0);
		}
		return cnt_LEs;
	}
	else
	{
		return cnt_LEs;
	}
}

static int cal_LEs(int N_of_Samples_of_LEs, int N_LE, int N_osc, double dt, double t, double t_Begin_Rec_LEs, double pert_size, double LEs[], double avg_LEs[], double norm[], double avg_avg_LEs[], int *N_of_Samples_of_avg_avg_LEs, double t_Begin_cal_avg_avg_LEs)
{
	int ret_N_of_Samples_of_LEs = N_of_Samples_of_LEs;

	if (t_Begin_Rec_LEs <= t)
	{
		int i;
		for (i = 0; i < N_LE; i++)
		{
			if (is_rescaling_pert == 1)
			{
				double local_LE = (1/dt)*(M_LOG(norm[i]) - M_LOG(pert_size));

				/* Update avg_LEs */
				avg_LEs[i] = (((double)N_of_Samples_of_LEs)*avg_LEs[i] + local_LE)/((double)(N_of_Samples_of_LEs + 1));

				/* It's more useful to keep the actual value of LE */
				LEs[i] = local_LE;

			}
			else
			{
				/* Check if norm[i] is too small or too big? */
				if ((fabs(norm[i]) <= TOO_SMALL_NUMBER) && (TOO_BIG_NUMBER <= fabs(norm[i])))
				{
					printf("[%s]: Pert. vector (%1.16lf) of LE %d goes too big or too small\n", getCurrentTime(), norm[i], i);
					is_exit_execution = 1;
				}
				else
				{
					double local_LE = (1/global_t)*(M_LOG(norm[i]) - M_LOG(pert_size));

					/* Update avg_LEs */
					avg_LEs[i] = local_LE;

					/* It's more useful to keep the actual value of LE */
					LEs[i] = local_LE;
				}
			}
		}

		if (is_exit_execution == 0)
		{
			ret_N_of_Samples_of_LEs += 1;
		}


		if (t_Begin_cal_avg_avg_LEs <= t)
		{
			if (is_exit_execution == 0)
			{
				int tmp_N_of_Samples_of_avg_avg_LEs = (*N_of_Samples_of_avg_avg_LEs);

				for (i = 0; i < N_LE; i++)
				{
					if (is_rescaling_pert == 1)
					{
						/* Update avg_avg_LEs */
						avg_avg_LEs[i] = (((double)tmp_N_of_Samples_of_avg_avg_LEs)*avg_avg_LEs[i] + avg_LEs[i])/((double)(tmp_N_of_Samples_of_avg_avg_LEs + 1));
					}
					else
					{
						avg_avg_LEs[i] = avg_LEs[i];
					}
				}

				*N_of_Samples_of_avg_avg_LEs = *N_of_Samples_of_avg_avg_LEs + 1;
			}
		}
	}

	return ret_N_of_Samples_of_LEs;
}

static void GSR(int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int is_connected_with_E[], int is_connected_with_I[], double norm[])
{
	int i;

	/* Cal. the first norm */
	norm_GSR(0, N_LE, N_osc, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, is_connected_with_E, is_connected_with_I, norm);

	/* GSR */
	for (i = 1; i < N_LE; i++)
	{
		/* Check for implementation later !!!!!!!!!!!!!!!!!!!!!!!!!! */
		if (is_change_pert_directions == 1)
		{
//			double dot[i], sum_RHS[N_osc];
//
//			dot_GSR(i, N_LE, N_osc, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, is_connected_with_E, is_connected_with_I, dot);
//
//			sum_GSR(i, N_LE, N_osc, v_perts, dot, norm, sum_RHS);
//
//			subtract_GSR(i, N_LE, N_osc, v_perts, sum_RHS);
		}

		norm_GSR(i, N_LE, N_osc, v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts, is_connected_with_E, is_connected_with_I, norm);
	}
}

static void norm_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int is_connected_with_E[], int is_connected_with_I[], double norm[])
{
	int i;

	norm[LE_i] = 0.0;
	for (i = 0; i < N_osc; i++)
	{
		norm[LE_i] += getAEle(v_perts, LE_i, i, N_LE, N_osc)*getAEle(v_perts, LE_i, i, N_LE, N_osc);

		if (is_connected_with_E[i] == 1)
		{
			norm[LE_i] += getAEle(s_E_perts, LE_i, i, N_LE, N_osc)*getAEle(s_E_perts, LE_i, i, N_LE, N_osc);
			norm[LE_i] += getAEle(x_E_perts, LE_i, i, N_LE, N_osc)*getAEle(x_E_perts, LE_i, i, N_LE, N_osc);
		}

		if (is_connected_with_I[i] == 1)
		{
			norm[LE_i] += getAEle(s_I_perts, LE_i, i, N_LE, N_osc)*getAEle(s_I_perts, LE_i, i, N_LE, N_osc);
			norm[LE_i] += getAEle(x_I_perts, LE_i, i, N_LE, N_osc)*getAEle(x_I_perts, LE_i, i, N_LE, N_osc);
		}
	}


	norm[LE_i] = sqrt(norm[LE_i]);
}

//static void dot_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], double is_connected_with_E[], double is_connected_with_I[], double dot[])
//{
//	int i;
//	for (i = 0; i <= (LE_i - 1); i++)
//	{
//		dot[i] = 0.0;
//
//		int j;
//		for (j = 0; j < N_osc; j++)
//		{
//			dot[i] += getAEle(v_perts, LE_i, j, N_LE, N_osc)*getAEle(v_perts, i, j, N_LE, N_osc);
//
//			if (is_connected_with_E[j] == 1)
//			{
//				dot[i] += getAEle(s_E_perts, LE_i, j, N_LE, N_osc)*getAEle(s_E_perts, i, j, N_LE, N_osc);
//				dot[i] += getAEle(x_E_perts, LE_i, j, N_LE, N_osc)*getAEle(x_E_perts, i, j, N_LE, N_osc);
//			}
//
//			if (is_connected_with_I[j] == 1)
//			{
//				dot[i] += getAEle(s_I_perts, LE_i, j, N_LE, N_osc)*getAEle(s_I_perts, i, j, N_LE, N_osc);
//				dot[i] += getAEle(x_I_perts, LE_i, j, N_LE, N_osc)*getAEle(x_I_perts, i, j, N_LE, N_osc);
//			}
//		}
//	}
//}

//static void sum_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], double is_connected_with_E[], double is_connected_with_I[], double dot[], double norm[], double sum_RHS[])
//{
//	int j;
//	for (j = 0; j < N_osc; j++)
//	{
//		sum_RHS[j] = 0.0;
//
//		int i;
//		for (i = 0; i <= (LE_i - 1); i++)
//		{
//			sum_RHS[j] += dot[i]*getAEle(v_perts, i, j, N_LE, N_osc)/(norm[i]*norm[i]);
//
//			if (is_connected_with_E[j] == 1)
//			{
//				sum_RHS[j] += getAEle(s_E_perts, LE_i, j, N_LE, N_osc)*getAEle(s_E_perts, i, j, N_LE, N_osc);
//				sum_RHS[j] += getAEle(x_E_perts, LE_i, j, N_LE, N_osc)*getAEle(x_E_perts, i, j, N_LE, N_osc);
//			}
//
//			if (is_connected_with_I[j] == 1)
//			{
//				sum_RHS[j] += getAEle(s_I_perts, LE_i, j, N_LE, N_osc)*getAEle(s_I_perts, i, j, N_LE, N_osc);
//				sum_RHS[j] += getAEle(x_I_perts, LE_i, j, N_LE, N_osc)*getAEle(x_I_perts, i, j, N_LE, N_osc);
//			}
//		}
//	}
//}

//static void subtract_GSR(int LE_i, int N_LE, int N_osc, double v_perts[], double sum_RHS[])
//{
//	int i;
//	for (i = 0; i < N_osc; i++)
//	{
//		setAEle(v_perts, getAEle(v_perts, LE_i, i, N_LE, N_osc) - sum_RHS[i], LE_i, i, N_LE, N_osc);
//	}
//}

static void init_perts(
		/* Input */
		int N_LE,
		int N_osc,
		double pert_size,
		int is_connected_with_E[], int is_connected_with_I[],
		/* Output */
		double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], int *N_states_of_pert_vector)
{
	if (is_use_rnd_init_perts_seed == 1)
	{
		srand48(rnd_init_perts_seed);
	}

	int i, j, isGoodCand;

	/* Determine the number of states */
	int tmp_N_states = N_osc;
	for (i = 0; i < N_osc; i++)
	{
		tmp_N_states += 2*is_connected_with_E[i] + 2*is_connected_with_I[i];
	}

	(*N_states_of_pert_vector) = tmp_N_states;

	double cand[tmp_N_states];
	double storages[N_LE * tmp_N_states];

	/* Determine linearly independent perturbations */
	for (i = 0; i < N_LE; i++)
	{
		isGoodCand = 0;
		while (isGoodCand == 0)
		{
			isGoodCand = 1;

			/* Initailize random */
			get_cand_pert(pert_size, tmp_N_states, cand);

			/* Check if the current perturbations are linearly dependent with the previous perturbations */
			for (j = 0; j <= i - 1; j++)
			{
				int k;
				double prev_perts[tmp_N_states];
				for (k = 0; k < tmp_N_states; k++)
				{
					prev_perts[k] = getAEle(storages, j, k, N_LE, tmp_N_states);
				}

				if (isLinearlyDependent(prev_perts, cand, tmp_N_states, 1e-14) == 1)
				{
					isGoodCand = 0;
					break;
				}
			}

			if (isGoodCand == 1)
			{
				/* Fill in the good candidate */
				for (j = 0; j < tmp_N_states; j++)
				{
					setAEle(storages, cand[j], i, j, N_LE, tmp_N_states);
				}
			}
		}
	}

	/* Distribute the storage of perturbations */
	for (i = 0; i < N_LE; i++)
	{
		int i_storages = 0;
		for (j = 0; j < N_osc; j++)
		{
			setAEle(v_perts, getAEle(storages, i, i_storages, N_LE, tmp_N_states), i, j, N_LE, N_osc);
			i_storages++;

			if (is_connected_with_E[j] == 1)
			{
				setAEle(s_E_perts, getAEle(storages, i, i_storages, N_LE, tmp_N_states), i, j, N_LE, N_osc);
				i_storages++;

				setAEle(x_E_perts, getAEle(storages, i, i_storages, N_LE, tmp_N_states), i, j, N_LE, N_osc);
				i_storages++;
			}

			if (is_connected_with_I[j] == 1)
			{
				setAEle(s_I_perts, getAEle(storages, i, i_storages, N_LE, tmp_N_states), i, j, N_LE, N_osc);
				i_storages++;

				setAEle(x_I_perts, getAEle(storages, i, i_storages, N_LE, tmp_N_states), i, j, N_LE, N_osc);
				i_storages++;
			}
		}
	}
}

static void get_cand_pert(double pert_size, int N_osc, double cand[])
{
	int j;

	for (j = 0; j < N_osc; j++)
	{
		cand[j] =  drand48();
	}

	double norm_cand = cal_norm(cand, N_osc);

	for (j = 0; j < N_osc; j++)
	{
		cand[j] =  cand[j]/norm_cand*pert_size;
	}
}

static int rec_Traj(int N_osc, double V[],
		S_elem_struct *S_max_connected_I, double g_syn_E_max_CI, double g_syn_I_max_CI, double g_syn_E_ext_max_CI, double g_syn_I_ext_max_CI, double V_syn_E_max_CI, double V_syn_I_max_CI, double V_syn_E_ext_max_CI, double V_syn_I_ext_max_CI,
		S_elem_struct *S_min_connected_I, double g_syn_E_min_CI, double g_syn_I_min_CI, double g_syn_E_ext_min_CI, double g_syn_I_ext_min_CI, double V_syn_E_min_CI, double V_syn_I_min_CI, double V_syn_E_ext_min_CI, double V_syn_I_ext_min_CI,
		S_elem_struct *S_max_connected_E, double g_syn_E_max_CE, double g_syn_I_max_CE, double g_syn_E_ext_max_CE, double g_syn_I_ext_max_CE, double V_syn_E_max_CE, double V_syn_I_max_CE, double V_syn_E_ext_max_CE, double V_syn_I_ext_max_CE,
		S_elem_struct *S_min_connected_E, double g_syn_E_min_CE, double g_syn_I_min_CE, double g_syn_E_ext_min_CE, double g_syn_I_ext_min_CE, double V_syn_E_min_CE, double V_syn_I_min_CE, double V_syn_E_ext_min_CE, double V_syn_I_ext_min_CE,
		FILE *pFileSpkProfile1, int max_connected_I, int min_connected_I, int max_connected_E, int min_connected_E, int cnt_traj)
{
	double V_max_CI, 			V_min_CI, 			V_max_CE, 			V_min_CE;
	double I_syn_E_max_CI,		I_syn_E_min_CI, 	I_syn_E_max_CE, 	I_syn_E_min_CE;
	double I_syn_I_max_CI,		I_syn_I_min_CI, 	I_syn_I_max_CE, 	I_syn_I_min_CE;
	double I_syn_E_ext_max_CI,	I_syn_E_ext_min_CI, I_syn_E_ext_max_CE, I_syn_E_ext_min_CE;
	double I_syn_I_ext_max_CI,	I_syn_I_ext_min_CI, I_syn_I_ext_max_CE, I_syn_I_ext_min_CE;
	double I_syn_max_CI, 		I_syn_min_CI, 		I_syn_max_CE, 		I_syn_min_CE;

	if (max_connected_I != NO_I_NEURONS)
	{
		V_max_CI			= V[max_connected_I];
		V_min_CI			= V[min_connected_I];

//		I_syn_E_max_CI		= cal_I_syn(g_syn_E_max_CI, S_max_connected_I->s_E, V_max_CI, V_syn_E_max_CI);
//		I_syn_E_min_CI		= cal_I_syn(g_syn_E_min_CI, S_min_connected_I->s_E, V_min_CI, V_syn_E_min_CI);
//		I_syn_I_max_CI		= cal_I_syn(g_syn_I_max_CI, S_max_connected_I->s_I, V_max_CI, V_syn_I_max_CI);
//		I_syn_I_min_CI		= cal_I_syn(g_syn_I_min_CI, S_min_connected_I->s_I, V_min_CI, V_syn_I_min_CI);
//		I_syn_E_ext_max_CI	= cal_I_syn(g_syn_E_ext_max_CI, S_max_connected_I->s_E_ext, V_max_CI, V_syn_E_ext_max_CI);
//		I_syn_E_ext_min_CI	= cal_I_syn(g_syn_E_ext_min_CI, S_min_connected_I->s_E_ext, V_min_CI, V_syn_E_ext_min_CI);
//		I_syn_I_ext_max_CI	= cal_I_syn(g_syn_I_ext_max_CI, S_max_connected_I->s_I_ext, V_max_CI, V_syn_I_ext_max_CI);
//		I_syn_I_ext_min_CI	= cal_I_syn(g_syn_I_ext_min_CI, S_min_connected_I->s_I_ext, V_min_CI, V_syn_I_ext_min_CI);
//		I_syn_max_CI		= I_syn_E_max_CI + I_syn_I_max_CI + I_syn_E_ext_max_CI + I_syn_I_ext_max_CI;
//		I_syn_min_CI		= I_syn_E_min_CI + I_syn_I_min_CI + I_syn_E_ext_min_CI + I_syn_I_ext_min_CI;

		I_syn_E_max_CI		= g_syn_E_max_CI*S_max_connected_I->s_E;
		I_syn_E_min_CI		= g_syn_E_min_CI*S_min_connected_I->s_E;
		I_syn_I_max_CI		= g_syn_I_max_CI*S_max_connected_I->s_I;
		I_syn_I_min_CI		= g_syn_I_min_CI*S_min_connected_I->s_I;
		I_syn_E_ext_max_CI	= g_syn_E_ext_max_CI*S_max_connected_I->s_E_ext;
		I_syn_E_ext_min_CI	= g_syn_E_ext_min_CI*S_min_connected_I->s_E_ext;
		I_syn_I_ext_max_CI	= g_syn_I_ext_max_CI*S_max_connected_I->s_I_ext;
		I_syn_I_ext_min_CI	= g_syn_I_ext_min_CI*S_min_connected_I->s_I_ext;
		I_syn_max_CI		= I_syn_E_max_CI + I_syn_I_max_CI + I_syn_E_ext_max_CI + I_syn_I_ext_max_CI;
		I_syn_min_CI		= I_syn_E_min_CI + I_syn_I_min_CI + I_syn_E_ext_min_CI + I_syn_I_ext_min_CI;
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		V_max_CE			= V[max_connected_E];
		V_min_CE			= V[min_connected_E];

//		I_syn_E_max_CE		= cal_I_syn(g_syn_E_max_CE, S_max_connected_E->s_E, V_max_CE, V_syn_E_max_CE);
//		I_syn_E_min_CE		= cal_I_syn(g_syn_E_min_CE, S_min_connected_E->s_E, V_min_CE, V_syn_E_min_CE);
//		I_syn_I_max_CE		= cal_I_syn(g_syn_I_max_CE, S_max_connected_E->s_I, V_max_CE, V_syn_I_max_CE);
//		I_syn_I_min_CE		= cal_I_syn(g_syn_I_min_CE, S_min_connected_E->s_I, V_min_CE, V_syn_I_min_CE);
//		I_syn_E_ext_max_CE	= cal_I_syn(g_syn_E_ext_max_CE, S_max_connected_E->s_E_ext, V_max_CE, V_syn_E_ext_max_CE);
//		I_syn_E_ext_min_CE	= cal_I_syn(g_syn_E_ext_min_CE, S_min_connected_E->s_E_ext, V_min_CE, V_syn_E_ext_min_CE);
//		I_syn_I_ext_max_CE	= cal_I_syn(g_syn_I_ext_max_CE, S_max_connected_E->s_I_ext, V_max_CE, V_syn_I_ext_max_CE);
//		I_syn_I_ext_min_CE	= cal_I_syn(g_syn_I_ext_min_CE, S_min_connected_E->s_I_ext, V_min_CE, V_syn_I_ext_min_CE);
//		I_syn_max_CE		= I_syn_E_max_CE + I_syn_I_max_CE + I_syn_E_ext_max_CE + I_syn_I_ext_max_CE;
//		I_syn_min_CE		= I_syn_E_min_CE + I_syn_I_min_CE + I_syn_E_ext_min_CE + I_syn_I_ext_min_CE;

		I_syn_E_max_CE		= g_syn_E_max_CE*S_max_connected_E->s_E;
		I_syn_E_min_CE		= g_syn_E_min_CE*S_min_connected_E->s_E;
		I_syn_I_max_CE		= g_syn_I_max_CE*S_max_connected_E->s_I;
		I_syn_I_min_CE		= g_syn_I_min_CE*S_min_connected_E->s_I;
		I_syn_E_ext_max_CE	= g_syn_E_ext_max_CE*S_max_connected_E->s_E_ext;
		I_syn_E_ext_min_CE	= g_syn_E_ext_min_CE*S_min_connected_E->s_E_ext;
		I_syn_I_ext_max_CE	= g_syn_I_ext_max_CE*S_max_connected_E->s_I_ext;
		I_syn_I_ext_min_CE	= g_syn_I_ext_min_CE*S_min_connected_E->s_I_ext;
		I_syn_max_CE		= I_syn_E_max_CE + I_syn_I_max_CE + I_syn_E_ext_max_CE + I_syn_I_ext_max_CE;
		I_syn_min_CE		= I_syn_E_min_CE + I_syn_I_min_CE + I_syn_E_ext_min_CE + I_syn_I_ext_min_CE;
	}

	// Voltages
	if (max_connected_I != NO_I_NEURONS)
	{
		cnt_traj=recResults(V_max_CI, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(V_min_CI, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		cnt_traj=recResults(V_max_CE, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(V_min_CE, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	// All synapses including the external ones
	if (max_connected_I != NO_I_NEURONS)
	{

		cnt_traj=recResults(I_syn_max_CI, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_min_CI, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		cnt_traj=recResults(I_syn_max_CE, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_min_CE, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	// Inputs from the I-cells
	if (max_connected_I != NO_I_NEURONS)
	{
		cnt_traj=recResults(I_syn_I_max_CI, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_I_min_CI, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		cnt_traj=recResults(I_syn_I_max_CE, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_I_min_CE, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	// Inputs from the E-cells
	if (max_connected_I != NO_I_NEURONS)
	{
		cnt_traj=recResults(I_syn_E_max_CI, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_E_min_CI, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		cnt_traj=recResults(I_syn_E_max_CE, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_E_min_CE, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	// External E inputs
	if (max_connected_I != NO_I_NEURONS)
	{
		cnt_traj=recResults(I_syn_E_ext_max_CI, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_E_ext_min_CI, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		cnt_traj=recResults(I_syn_E_ext_max_CE, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_E_ext_min_CE, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	// External I inputs
	if (max_connected_I != NO_I_NEURONS)
	{
		cnt_traj=recResults(I_syn_I_ext_max_CI, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_I_ext_min_CI, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	if (max_connected_E != NO_E_NEURONS)
	{
		cnt_traj=recResults(I_syn_I_ext_max_CE, pFileSpkProfile1, cnt_traj, 0);
		cnt_traj=recResults(I_syn_I_ext_min_CE, pFileSpkProfile1, cnt_traj, 0);
	}
	else
	{
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
		cnt_traj=recResults(0.0,pFileSpkProfile1,cnt_traj,0);
	}

	return cnt_traj;
}

static void max_min_connected_Osc(
		S_elem_struct *head_S,
		S_elem_struct **S_max_connected_I, double *g_syn_E_max_CI, double *g_syn_I_max_CI, double *g_syn_E_ext_max_CI, double *g_syn_I_ext_max_CI, double *V_syn_E_max_CI, double *V_syn_I_max_CI, double *V_syn_E_ext_max_CI, double *V_syn_I_ext_max_CI,
		S_elem_struct **S_min_connected_I, double *g_syn_E_min_CI, double *g_syn_I_min_CI, double *g_syn_E_ext_min_CI, double *g_syn_I_ext_min_CI, double *V_syn_E_min_CI, double *V_syn_I_min_CI, double *V_syn_E_ext_min_CI, double *V_syn_I_ext_min_CI,
		S_elem_struct **S_max_connected_E, double *g_syn_E_max_CE, double *g_syn_I_max_CE, double *g_syn_E_ext_max_CE, double *g_syn_I_ext_max_CE, double *V_syn_E_max_CE, double *V_syn_I_max_CE, double *V_syn_E_ext_max_CE, double *V_syn_I_ext_max_CE,
		S_elem_struct **S_min_connected_E, double *g_syn_E_min_CE, double *g_syn_I_min_CE, double *g_syn_E_ext_min_CE, double *g_syn_I_ext_min_CE, double *V_syn_E_min_CE, double *V_syn_I_min_CE, double *V_syn_E_ext_min_CE, double *V_syn_I_ext_min_CE,
		int *max_connected_I, int *min_connected_I, int *max_connected_E, int *min_connected_E,
		network_props_struct *network_p, post_syn_osc_struct *post_syn_Mat[],
		syn_props_struct *syn_props, ext_input_props_struct *ext_input_props)
{
	int i;
	int	N_e = network_p->N_e;
	int N_i = network_p->N_i;
	S_elem_struct *tmp_S;

	int N_pre[N_i + N_e];
	for (i = 0; i < N_i + N_e; i++)
	{
		N_pre[i] = 0;

		post_syn_osc_struct	*head_post_syn_osc = post_syn_Mat[i];
		while (head_post_syn_osc != NULL)
		{
			N_pre[i]++;
			head_post_syn_osc = head_post_syn_osc->next;
		}
	}

	(*max_connected_I) = 0;					// Osc. I ID that has the most number of presynaptic neurons.
	(*min_connected_I) = N_i - 1;			// Osc. I ID that has the least number of presynaptic neurons.
	(*max_connected_E) = N_i;				// Osc. E ID that has the most number of presynaptic neurons.
	(*min_connected_E) = N_i + N_e - 1;		// Osc. E ID that has the least number of presynaptic neurons.

	tmp_S = head_S;
	for (i = 0; i < N_i + N_e; i++)
	{
		if (i == 0)
		{
			(*S_max_connected_I) = tmp_S;
		}

		if (i == (N_i - 1))
		{
			(*S_min_connected_I) = tmp_S;
		}

		if (i == N_i)
		{
			(*S_max_connected_E) = tmp_S;
		}

		if (i == (N_i + N_e - 1))
		{
			(*S_min_connected_E) = tmp_S;
		}

		tmp_S = tmp_S->next;
	}

	tmp_S = head_S;
	for (i = 0; i < N_i; i++)
	{
		if (N_pre[(*max_connected_I)] < N_pre[i])
		{
			(*max_connected_I) = i;
			(*S_max_connected_I) = tmp_S;
		}

		if (N_pre[i] < N_pre[(*min_connected_I)])
		{
			(*min_connected_I) = i;
			(*S_min_connected_I) = tmp_S;
		}

		tmp_S = tmp_S->next;
	}

	if (N_i == 0)
	{
		(*max_connected_I) 	= NO_I_NEURONS;
		(*min_connected_I) 	= NO_I_NEURONS;
		(*S_max_connected_I)= NULL;
		(*S_min_connected_I)= NULL;
	}

	for (i = N_i; i < N_i + N_e; i++)
	{
		if (N_pre[(*max_connected_E)] < N_pre[i])
		{
			(*max_connected_E) = i;
			(*S_max_connected_E) = tmp_S;
		}

		if (N_pre[i] < N_pre[(*min_connected_E)])
		{
			(*min_connected_E) = i;
			(*S_min_connected_E) = tmp_S;
		}

		tmp_S = tmp_S->next;
	}

	if (N_e == 0)
	{
		(*max_connected_E) 	= NO_E_NEURONS;
		(*min_connected_E) 	= NO_E_NEURONS;
		(*S_max_connected_E)= NULL;
		(*S_min_connected_E)= NULL;
	}

	/* Get the kinetics of synpase */
	double dummy;
	if (N_i != 0)
	{
		get_synapse_kinetics(*max_connected_I, N_i,
				g_syn_E_max_CI, g_syn_I_max_CI, g_syn_E_ext_max_CI, g_syn_I_ext_max_CI,
				V_syn_E_max_CI, V_syn_I_max_CI, V_syn_E_ext_max_CI, V_syn_I_ext_max_CI,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				syn_props, ext_input_props);

		get_synapse_kinetics(*min_connected_I, N_i,
				g_syn_E_min_CI, g_syn_I_min_CI, g_syn_E_ext_min_CI, g_syn_I_ext_min_CI,
				V_syn_E_min_CI, V_syn_I_min_CI, V_syn_E_ext_min_CI, V_syn_I_ext_min_CI,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				syn_props, ext_input_props);
	}

	if (N_e != 0)
	{
		get_synapse_kinetics(*max_connected_E, N_i,
				g_syn_E_max_CE, g_syn_I_max_CE, g_syn_E_ext_max_CE, g_syn_I_ext_max_CE,
				V_syn_E_max_CE, V_syn_I_max_CE, V_syn_E_ext_max_CE, V_syn_I_ext_max_CE,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				syn_props, ext_input_props);

		get_synapse_kinetics(*min_connected_E, N_i,
				g_syn_E_min_CE, g_syn_I_min_CE, g_syn_E_ext_min_CE, g_syn_I_ext_min_CE,
				V_syn_E_min_CE, V_syn_I_min_CE, V_syn_E_ext_min_CE, V_syn_I_ext_min_CE,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				syn_props, ext_input_props);
	}
}

static int rec_Params(double wb_V_pk, double wb_Icells_pk[], double mp_V_pk, double mp_Ecells_pk[], int N_e, int N_i, double mp_Ecells[], double wb_Icells[], double v[], FILE *pFileSpkProfile,double RmIe_E, double RmIe_I, syn_props_struct *syn_p,E_neuron_props_struct *E_neuron_p,I_neuron_props_struct *I_neuron_p,network_props_struct *network_p,ext_input_props_struct *ext_input_props,sim_props_struct *sim_p,int cnt)
{

	cnt=recResults(syn_p->g_syn_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->g_syn_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->g_syn_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->g_syn_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(syn_p->v_syn_AMPA,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->v_syn_GABA,pFileSpkProfile,cnt,1);

	cnt=recResults(syn_p->tau_m_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_m_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(syn_p->tau_l_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_l_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_l_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_l_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(syn_p->tau_r_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_r_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_r_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_r_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(syn_p->tau_d_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_d_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_d_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(syn_p->tau_d_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(E_neuron_p->V_rest,pFileSpkProfile,cnt,1);
	cnt=recResults(E_neuron_p->V_th,pFileSpkProfile,cnt,1);
	cnt=recResults(E_neuron_p->V_reset,pFileSpkProfile,cnt,1);
	cnt=recResults(E_neuron_p->V_refrac,pFileSpkProfile,cnt,1);
	cnt=recResults(E_neuron_p->tau_m,pFileSpkProfile,cnt,1);
	cnt=recResults(E_neuron_p->C_m,pFileSpkProfile,cnt,1);

	cnt=recResults(I_neuron_p->V_rest,pFileSpkProfile,cnt,1);
	cnt=recResults(I_neuron_p->V_th,pFileSpkProfile,cnt,1);
	cnt=recResults(I_neuron_p->V_reset,pFileSpkProfile,cnt,1);
	cnt=recResults(I_neuron_p->V_refrac,pFileSpkProfile,cnt,1);
	cnt=recResults(I_neuron_p->tau_m,pFileSpkProfile,cnt,1);
	cnt=recResults(I_neuron_p->C_m,pFileSpkProfile,cnt,1);

	cnt=recResults((double)network_p->N_i,pFileSpkProfile,cnt,1);
	cnt=recResults((double)network_p->N_e,pFileSpkProfile,cnt,1);
	cnt=recResults((double)network_p->N_ext,pFileSpkProfile,cnt,1);

	cnt=recResults(network_p->p_EE,pFileSpkProfile,cnt,1);
	cnt=recResults(network_p->p_EI,pFileSpkProfile,cnt,1);
	cnt=recResults(network_p->p_IE,pFileSpkProfile,cnt,1);
	cnt=recResults(network_p->p_II,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->lambda_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->lambda_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->lambda_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->lambda_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->v_syn_AMPA,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->v_syn_GABA,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->g_ext_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->g_ext_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->g_ext_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->g_ext_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->tau_m_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_m_on_E,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->tau_l_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_l_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_l_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_l_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->tau_r_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_r_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_r_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_r_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(ext_input_props->tau_d_AMPA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_d_AMPA_on_I,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_d_GABA_on_E,pFileSpkProfile,cnt,1);
	cnt=recResults(ext_input_props->tau_d_GABA_on_I,pFileSpkProfile,cnt,1);

	cnt=recResults(RmIe_E,pFileSpkProfile,cnt,1);
	cnt=recResults(RmIe_I,pFileSpkProfile,cnt,1);

	cnt=recResults(sim_p->dt,pFileSpkProfile,cnt,1);
	cnt=recResults((double)sim_p->Nt,pFileSpkProfile,cnt,1);
	cnt=recResults((double)sim_p->seed,pFileSpkProfile,cnt,1);
	cnt=recResults((double)sim_p->isUsingSeed,pFileSpkProfile,cnt,1);

//	// Keep last state
//	if (0 < N_i)
//	{
//		cnt=recResults(v[0],pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < wb_N_states; j++)
//		{
//			cnt=recResults(wb_Icells[j],pFileSpkProfile,cnt,1);
//		}
//	}
//	else
//	{
//		cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < wb_N_states; j++)
//		{
//			cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//		}
//	}
//
//	if (0 < N_e)
//	{
//		cnt=recResults(v[N_i],pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < nw_N_states; j++)
//		{
//			cnt=recResults(mp_Ecells[j],pFileSpkProfile,cnt,1);
//		}
//	}
//	else
//	{
//		cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < nw_N_states; j++)
//		{
//			cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//		}
//	}
//

//	// Keep peak
//	if (0 < N_i)
//	{
//		cnt=recResults(wb_V_pk,pFileSpkProfile,cnt,1);
//
//#if (INTN_TYPE == WANG_BUZSAKI)
//		int j;
//		for (j = 0; j < wb_N_states; j++)
//		{
//			cnt=recResults(wb_Icells_pk[j],pFileSpkProfile,cnt,1);
//		}
//#endif
//#if (INTN_TYPE == BORGER_WALKER)
//		int j;
//		for (j = 0; j < bw_N_states; j++)
//		{
//			cnt=recResults(wb_Icells_pk[j],pFileSpkProfile,cnt,1);
//		}
//#endif
//	}


//	// Keep peak
//	if (0 < N_i)
//	{
//		cnt=recResults(wb_V_pk,pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < wb_N_states; j++)
//		{
//			cnt=recResults(wb_Icells_pk[j],pFileSpkProfile,cnt,1);
//		}
//	}
//	else
//	{
//		cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < wb_N_states; j++)
//		{
//			cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//		}
//	}
//
//	if (0 < N_e)
//	{
//		cnt=recResults(mp_V_pk,pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < nw_N_states; j++)
//		{
//			cnt=recResults(mp_Ecells_pk[j],pFileSpkProfile,cnt,1);
//		}
//	}
//	else
//	{
//		cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//
//		int j;
//		for (j = 0; j < nw_N_states; j++)
//		{
//			cnt=recResults(0.0,pFileSpkProfile,cnt,1);
//		}
//	}


	return cnt;
}

static void determineConnections(int is_connected_with_E[], int is_connected_with_I[],
		post_syn_osc_struct *post_syn_Mat[], post_syn_osc_struct *pre_syn_Mat[],
		double p_EE,double p_EI,double p_IE,double p_II,
		double p_EE_gap_junction, double p_EI_gap_junction, double p_IE_gap_junction, double p_II_gap_junction,
		int N_i,int N_e,
		int isGJExactNumberOfConns)
{
	if (is_use_rnd_determineConnections == 1)
	{
		srand48(rnd_determineConnections);
	}

	int n_ee = 0, n_ei = 0, n_ie = 0, n_ii = 0;

	/* Chemical synapse, post_syn_Mat contains lists of postsynaptic neurons */
	int post, pre;
	for (pre = 0; pre < N_i + N_e; pre++)				// Pre
	{
		post_syn_osc_struct *head_post_syn_osc = post_syn_Mat[pre];

		for (post = 0; post < N_i + N_e; post++)		// Post
		{
			if (pre < N_i) // From I-cells
			{
				if (!((allow_II_connections == IS_SELF_II_NOTALLOW) && (pre == post)))
				{
					double rand_number = drand48();

					int isConnected=0;
					if (post < N_i)
					{
						if (rand_number < p_II)
							isConnected=1;
					}
					else
					{
						if (rand_number < p_IE)
							isConnected=1;
					}

					if (isConnected==1)
					{
						if (post < N_i)
						{
							n_ii += 1;
						}
						else
						{
							n_ie += 1;
						}

						is_connected_with_I[post] = 1;

						post_syn_osc_struct	*tmp	=(post_syn_osc_struct *)malloc(sizeof(post_syn_osc_struct));
						tmp->pre					=NULL;
						tmp->id						=post;
						tmp->next					=NULL;

						if (head_post_syn_osc == NULL)
						{
							post_syn_Mat[pre]	=tmp;
							head_post_syn_osc	=tmp;
						}
						else
						{
							head_post_syn_osc->next			=tmp;
							tmp->pre			=head_post_syn_osc;
							head_post_syn_osc				=tmp;
						}
					}
				}
			}
			else // From E-cells
			{
				if (!((allow_EE_connections == IS_SELF_EE_NOTALLOW) && (pre == post)))
				{
					double rand_number = drand48();

					int isConnected=0;
					if (post < N_i)
					{
						if (rand_number < p_EI)
							isConnected=1;
					}
					else
					{
						if (rand_number < p_EE)
							isConnected=1;
					}

					if (isConnected==1)
					{
						if (post < N_i)
						{
							n_ei += 1;
						}
						else
						{
							n_ee += 1;
						}

						is_connected_with_E[post] = 1;

						post_syn_osc_struct	*tmp	=(post_syn_osc_struct *)malloc(sizeof(post_syn_osc_struct));
						tmp->pre					=NULL;
						tmp->id						=post;
						tmp->next					=NULL;

						if (head_post_syn_osc == NULL)
						{
							post_syn_Mat[pre]	=tmp;
							head_post_syn_osc	=tmp;
						}
						else
						{
							head_post_syn_osc->next			=tmp;
							tmp->pre			=head_post_syn_osc;
							head_post_syn_osc				=tmp;
						}
					}
				}
			}
		}
	}

	/* Gap junction, pre_syn_Mat contains lists of presynaptic neurons.  Gap junction is symmetric then the connections are bi-directional. When p_II_gap_junction == 1.0, a neuron connects to N_i - 1 neurons. */
	int N_GJ_II_NeighborNeurons_realization[N_i];
	int i;

	for (i = 0; i < N_i; i++)
	{
		N_GJ_II_NeighborNeurons_realization[i] = 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (isGJExactNumberOfConns == 1)
	{
		printf("[%s] Exact GJ connections\n", getCurrentTime());

		double N_GJ_II_NeighborNeurons = ceil(p_II_gap_junction*(N_i - 1));
		double Quota_2_NeighborNeurons[N_i];

		for (i = 0; i < N_i; i++)
		{
			Quota_2_NeighborNeurons[i] = N_GJ_II_NeighborNeurons;
		}


		for (post = 0; post < N_i; post++)
		{
			// Construct an available list
			post_syn_osc_struct *curr_available_neurons = NULL, *head_available_neurons = NULL;
			int j, N_available_NeighborNeurons = 0;
			for (j = 0; j < N_i; j++)
			{
				if ((post != j) && (cmp(0.0, Quota_2_NeighborNeurons[j], 1e-6) < 0))
				{
					N_available_NeighborNeurons += 1;

					post_syn_osc_struct	*tmp	=(post_syn_osc_struct *)malloc(sizeof(post_syn_osc_struct));
					tmp->pre					=NULL;
					tmp->id						=j;
					tmp->next					=NULL;

					if (curr_available_neurons == NULL)
					{
						head_available_neurons = tmp;
						curr_available_neurons = tmp;
					}
					else
					{
						curr_available_neurons->next 	=tmp;
						tmp->pre						=curr_available_neurons;
						curr_available_neurons			=tmp;
					}
				}
			}

			if (0 < N_available_NeighborNeurons)
			{
				// Random priority of the neighbor neurons
				double 	rand_double[N_available_NeighborNeurons];
				int		rand_int[N_available_NeighborNeurons];

				int ii;
				for (ii = 0; ii < N_available_NeighborNeurons; ii++)
				{
					rand_double[ii] = drand48();
					rand_int[ii] = ii;
				}

				// Sort the neighbor neurons from low to high
				int jj;
				for (ii = 0; ii < N_available_NeighborNeurons; ii++)
				{
					for(jj = ii + 1; jj < N_available_NeighborNeurons; jj++)
					{
						if(rand_double[ii] > rand_double[jj])
						{
							double temp = rand_double[ii];
							rand_double[ii] = rand_double[jj];
							rand_double[jj] = temp;

							int tmp_int = rand_int[ii];
							rand_int[ii] = rand_int[jj];
							rand_int[jj] = tmp_int;
						}
					}
				}

				// Pick up the last N_available_NeighborNeurons neurons
				if (cmp(Quota_2_NeighborNeurons[post], (double)N_available_NeighborNeurons, 1e-6) < 0)
				{
					int iii = N_available_NeighborNeurons - 1;
					while (cmp(0.0, Quota_2_NeighborNeurons[post], 1e-3) < 0)
					{
						// Select a desired neuron
						post_syn_osc_struct *tmp_available_neurons = head_available_neurons;
						int iiii;
						for (iiii = 0; iiii < rand_int[iii]; iiii++)
						{
							tmp_available_neurons = tmp_available_neurons->next;
						}

						connect_Neurons_via_GJ(pre_syn_Mat, post, tmp_available_neurons->id);
						Quota_2_NeighborNeurons[post] -= 1.0;
						Quota_2_NeighborNeurons[tmp_available_neurons->id] -= 1.0;

						N_GJ_II_NeighborNeurons_realization[post] += 1;
						N_GJ_II_NeighborNeurons_realization[tmp_available_neurons->id] += 1;

						iii -= 1;
					}
				}
				else
				{
					post_syn_osc_struct *tmp_available_neurons = head_available_neurons;

					while (tmp_available_neurons != NULL)
					{
						connect_Neurons_via_GJ(pre_syn_Mat, post, tmp_available_neurons->id);
						Quota_2_NeighborNeurons[post] -= 1.0;
						Quota_2_NeighborNeurons[tmp_available_neurons->id] -= 1.0;

						N_GJ_II_NeighborNeurons_realization[post] += 1;
						N_GJ_II_NeighborNeurons_realization[tmp_available_neurons->id] += 1;

						tmp_available_neurons = tmp_available_neurons->next;
					}

				}

				// Free memory
				while (head_available_neurons != NULL)
				{
					post_syn_osc_struct	*tmp = head_available_neurons->next;

					free(head_available_neurons);

					head_available_neurons = tmp;
				}
			}
		}
	}
	else
	{
		printf("[%s] Approximate GJ connections\n", getCurrentTime());

		for (post = 0; post < N_i; post++)
		{
			for (pre = (post + 1); pre < N_i; pre++)
			{
				if ((drand48() <= p_II_gap_junction) && (connect_Neurons_via_GJ(pre_syn_Mat, post, pre) == 1))
				{
					N_GJ_II_NeighborNeurons_realization[post] += 1;
					N_GJ_II_NeighborNeurons_realization[pre] += 1;
				}
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double total_NGJ_II = 0;
	for (i = 0; i < N_i; i++)
	{
		total_NGJ_II += (double)N_GJ_II_NeighborNeurons_realization[i];
	}
	double avg_number_of_branch_for_each_neuron = total_NGJ_II/((double)N_i);
	total_NGJ_II = total_NGJ_II/2.0;

//	double total_number_of_nodes = total_NGJ_II*2.0;
//	double pred_avg_number_of_branch_for_each_neuron = total_number_of_nodes/((double)N_i);

//	for (i = 0; i < N_i; i++)
//	{
//		post_syn_osc_struct *tmp_available_neurons = pre_syn_Mat[i];
//
//		while (tmp_available_neurons != NULL)
//		{
//			printf("A(%d,%d) = A(%d,%d) + 1;\n", i + 1, tmp_available_neurons->id + 1, i + 1, tmp_available_neurons->id + 1);
//			tmp_available_neurons = tmp_available_neurons->next;
//		}
//	}

	printf("Chemical synapses	:EE=%d (p=%1.16lf), EI=%d (p=%1.16lf), IE=%d (p=%1.16lf), II=%d (p=%1.16lf)\n",
			n_ee, (((double)n_ee)/((double)N_e))/((double)N_e),
			n_ei, (((double)n_ei)/((double)N_e))/((double)N_i),
			n_ie, (((double)n_ie)/((double)N_i))/((double)N_e),
			n_ii, (((double)n_ii)/((double)N_i))/((double)N_i));

	printf("Gap junctions		:II=%1.1lf (pred=%1.1lf), #PreCells=%1.1lf (exp=%1.1lf)\n",
			total_NGJ_II, p_II_gap_junction*(((double)N_i)*(((double)N_i) - 1.0))/2.0,
			avg_number_of_branch_for_each_neuron, p_II_gap_junction*((double)N_i - 1.0));
//	printf("%1.16lf, %1.16lf\n", avg_number_of_branch_for_each_neuron, pred_avg_number_of_branch_for_each_neuron/(((double)N_i) - 1.0));
}

static int is_GJ_connected(int post, int pre, post_syn_osc_struct *pre_syn_Mat[])
{
	post_syn_osc_struct *head_pre_syn_osc = pre_syn_Mat[post];

	while (head_pre_syn_osc != NULL)
	{
		if (head_pre_syn_osc->id == pre)
		{
			return 1;
		}

		head_pre_syn_osc = head_pre_syn_osc->next;
	}

	return 0;
}

static int connect_Neurons_via_GJ(post_syn_osc_struct *pre_syn_Mat[], int post, int pre)
{
	if (is_GJ_connected(post, pre, pre_syn_Mat) == 0)
	{
		post_syn_osc_struct *head_pre_syn_osc = NULL, *tmp = NULL;

		// Connect post to pre neurons
		head_pre_syn_osc = pre_syn_Mat[post];

		tmp	=(post_syn_osc_struct *)malloc(sizeof(post_syn_osc_struct));
		tmp->pre	=NULL;
		tmp->id		=pre;
		tmp->next	=NULL;

		if (head_pre_syn_osc == NULL)
		{
			pre_syn_Mat[post] = tmp;
		}
		else
		{
			while (head_pre_syn_osc->next != NULL)
			{
				head_pre_syn_osc = head_pre_syn_osc->next;
			}

			head_pre_syn_osc->next	=tmp;
			tmp->pre				=head_pre_syn_osc;
		}

		// Connect pre to post neurons
		head_pre_syn_osc = pre_syn_Mat[pre];

		tmp	=(post_syn_osc_struct *)malloc(sizeof(post_syn_osc_struct));
		tmp->pre	=NULL;
		tmp->id		=post;
		tmp->next	=NULL;

		if (head_pre_syn_osc == NULL)
		{
			pre_syn_Mat[pre] = tmp;
		}
		else
		{
			while (head_pre_syn_osc->next != NULL)
			{
				head_pre_syn_osc = head_pre_syn_osc->next;
			}

			head_pre_syn_osc->next			=tmp;
			tmp->pre						=head_pre_syn_osc;
		}

		return 1;
	}
	else
	{
		return 0;
	}
}

int	recResults(double data,FILE *pFile,int cnt,int isAppendSeparator)
{
	double	datum[1];

	datum[0] = data;

	if (fseek(pFile, sizeof(double)*cnt, SEEK_SET) != 0)
	{
		printf("[recResults] cannot seek and exit(EXIT_FAILURE)\n");
		exit(EXIT_FAILURE);
	}

	if (fwrite(datum, sizeof(double), 1, pFile) != 1)
	{
		printf("[recResults] cannot write and exit(EXIT_FAILURE)\n");
		exit(EXIT_FAILURE);
	}

	cnt++;

	if (isAppendSeparator==1)
	{
		datum[0] = -1.0;	// Rec. separation.

		if (fseek(pFile, sizeof(double)*cnt, SEEK_SET) != 0)
		{
			printf("[recResults] cannot seek and exit(EXIT_FAILURE)\n");
			exit(EXIT_FAILURE);
		}

		if (fwrite(datum, sizeof(double), 1, pFile) != 1)
		{
			printf("[recResults] cannot write and exit(EXIT_FAILURE)\n");
			exit(EXIT_FAILURE);
		}

		cnt++;
	}

	return cnt;
}

void free_struct_arrays(void *head_Arrays[],void *tail_Arrays[],int N,int type)
{
	int i;
	for (i=0;i<N;i++)
	{
		void *tmp=head_Arrays[i];

		while (tmp!=NULL)
		{
			void *tmp1=tmp;

			switch(type)
			{
			case 0:
				tmp=(void *)(((input_struct *)tmp)->next);
				break;
			case 1 :
				tmp=(void *)(((spiking_struct *)tmp)->next);
				break;
			case 2 :
				tmp=(void *)(((post_syn_osc_struct *)tmp)->next);
				break;
			}

			free(tmp1);
		}

		head_Arrays[i]=NULL;
		if (tail_Arrays!=NULL)
		{
			tail_Arrays[i]=NULL;
		}
	}
}

static void addSpikingInfo(int i, double t, spiking_struct *head_spikingArrays[], spiking_struct *tail_spikingArrays[], int N)
{
	/* Get memory of spiking variable. */
	spiking_struct *tmp	=(spiking_struct *)malloc(sizeof(spiking_struct));
	tmp->pre			=NULL;
	tmp->time			=t;
	tmp->next			=NULL;

	if (tail_spikingArrays[i] == NULL)
	{
		head_spikingArrays[i] = tmp;
		tail_spikingArrays[i] = tmp;
	}
	else
	{
		tmp->pre					=tail_spikingArrays[i];
		tail_spikingArrays[i]->next	=tmp;

		tail_spikingArrays[i]		=tmp;
	}
}

static void addSpkArrivalArrays(int N_LE, double v_perts[], double dvdt, int pre, double t, input_struct *head_inputArrays[], input_struct *tail_inputArrays[], post_syn_osc_struct *post_syn_Mat[], syn_props_struct *syn_p, network_props_struct *network_p)
{
	int N_i = network_p->N_i;
	int post;

	post_syn_osc_struct *head_post_syn_osc = post_syn_Mat[pre];
	while (head_post_syn_osc != NULL)
	{
		post = head_post_syn_osc->id;

		int is_dist_latency;
		double tau_l, sigma_tau_l;
		if (pre < N_i)
		{
			if (post < N_i)
			{
				tau_l = syn_p->tau_l_GABA_on_I;
				sigma_tau_l = sigma_II_latency;
				is_dist_latency = is_dist_II_latency;
			}
			else
			{
				tau_l = syn_p->tau_l_GABA_on_E;
				sigma_tau_l = sigma_IE_latency;
				is_dist_latency = is_dist_IE_latency;
			}
		}
		else
		{
			if (post < N_i)
			{
				tau_l = syn_p->tau_l_AMPA_on_I;
				sigma_tau_l = sigma_EI_latency;
				is_dist_latency = is_dist_EI_latency;
			}
			else
			{
				tau_l = syn_p->tau_l_AMPA_on_E;
				sigma_tau_l = sigma_EE_latency;
				is_dist_latency = is_dist_EE_latency;
			}
		}

		if (is_dist_latency == 1)
		{
			tau_l = get_effective_latency(tau_l, sigma_tau_l);
		}

		/* Store perturbations of voltage at the spiking time */
		int j, N_osc = network_p->N_i + network_p->N_e;

		v_pert_elem_struct *head_v_pert_elem = (v_pert_elem_struct *)malloc(sizeof(v_pert_elem_struct));
		head_v_pert_elem->prev	=NULL;
		head_v_pert_elem->dv	=getAEle(v_perts, 0, pre, N_LE, N_osc);
		head_v_pert_elem->next	=NULL;

		v_pert_elem_struct *tmp_head_v_pert_elem = head_v_pert_elem;
		for (j = 1; j < N_LE; j++)
		{
			v_pert_elem_struct *tmp = (v_pert_elem_struct *)malloc(sizeof(v_pert_elem_struct));
			tmp->prev	=NULL;
			tmp->dv		=getAEle(v_perts, j, pre, N_LE, N_osc);
			tmp->next	=NULL;

			tmp_head_v_pert_elem->next = tmp;
			tmp->prev = tmp_head_v_pert_elem;

			tmp_head_v_pert_elem = tmp;
		}

		add_Input(head_v_pert_elem, dvdt, t, 1.0, head_inputArrays, tail_inputArrays, tau_l, pre, post, is_all_delays_the_same);

		head_post_syn_osc = head_post_syn_osc->next;
	}
}

static void addSpkArrivalArrays_Ext_inputs(double t, input_struct *head_inputArrays[], input_struct *tail_inputArrays[], network_props_struct *network_p, ext_input_props_struct *ext_input_props, sim_props_struct *sim_p)
{
	double dt			= sim_p->dt;
	int N_i				= network_p->N_i;
	int N_e				= network_p->N_e;
	int N_ext			= network_p->N_ext;

	int post, pre;
	for (post = 0; post < N_i + N_e; post++)
	{
		double N_input_AMPA = 0.0	, N_input_GABA = 0.0;
		double tau_l_AMPA = 0.0		, tau_l_GABA = 0.0;
		double lambda_AMPA = 0.0	, lambda_GABA = 0.0;

		int	   n_periodic_AMPA_external_inputs = 0, n_periodic_GABA_external_inputs = 0;

		/* Pick parameters */
		if (post < N_i)
		{
			// Excitation
			tau_l_AMPA	= ext_input_props->tau_l_AMPA_on_I;
			lambda_AMPA	= ext_input_props->lambda_AMPA_on_I;
			n_periodic_AMPA_external_inputs = n_periodic_AMPA_external_inputs_to_I;

			// Inhibition
			tau_l_GABA 	= ext_input_props->tau_l_GABA_on_I;
			lambda_GABA	= ext_input_props->lambda_GABA_on_I;
			n_periodic_GABA_external_inputs = n_periodic_GABA_external_inputs_to_I;
		}
		else
		{
			// Excitation
			tau_l_AMPA	= ext_input_props->tau_l_AMPA_on_E;
			lambda_AMPA	= ext_input_props->lambda_AMPA_on_E;
			n_periodic_AMPA_external_inputs = n_periodic_AMPA_external_inputs_to_E;

			// Inhibition
			tau_l_GABA 	= ext_input_props->tau_l_GABA_on_E;
			lambda_GABA	= ext_input_props->lambda_GABA_on_E;
			n_periodic_GABA_external_inputs = n_periodic_GABA_external_inputs_to_E;
		}

		/* Choose how to generate the external inputs */
		switch(type_of_external_inputs)
		{
		    case POISSON_EXTERNAL_INPUTS:
		    	for (pre = 0; pre < N_ext; pre++)
		    	{
		    		// Excitation
		    		double rand_nb = drand48();
		    		if (rand_nb < lambda_AMPA*dt)
		    		{
		    			N_input_AMPA = 1.0;
		    		}

		    		// Inhibition
		    		rand_nb = drand48();
		    		if (rand_nb < lambda_GABA*dt)
		    		{
		    			N_input_GABA = 1.0;
		    		}
		    	}

		        break;
		    case PERIODIC_EXTERNAL_INPUTS:
		    	// Excitation
				if (time_shift + n_periodic_AMPA_external_inputs/lambda_AMPA <= global_t)
				{
					N_input_AMPA = N_ext;
				}

				// Inhibition
				if (time_shift + n_periodic_GABA_external_inputs/lambda_GABA <= global_t)
				{
					N_input_GABA = N_ext;
				}

		        break;
		    case PERIODIC_POISSON_EXTERNAL_INPUTS:
		    	// Excitation
				lambda_AMPA += periodic_poisson_A*sin(2*PI*periodic_poisson_f*(global_t + periodic_poisson_time_shift)/1000);
				for (pre = 0; pre < N_ext; pre++)
				{
					double rand_nb = drand48();
					if (rand_nb < lambda_AMPA*dt)
					{
						N_input_AMPA = 1.0;
					}
				}

				// Inhibition
				lambda_GABA += periodic_poisson_A*sin(2*PI*periodic_poisson_f*(global_t + periodic_poisson_time_shift)/1000);
				for (pre = 0; pre < N_ext; pre++)
				{
					double rand_nb = drand48();
					if (rand_nb < lambda_GABA*dt)
					{
						N_input_GABA = 1.0;
					}
				}

		        break;
		    default:
		        printf("addSpkArrivalArrays_Ext_inputs: meet non-option and exit(0)\n");
		        exit(0);
		        break;
		}

		/* Add the generated external inputs to the queues */
		// Excitation
		if (cmp(0.0, N_input_AMPA, 1e-14) < 0)
		{
			add_Input(NULL, 0, t, N_input_AMPA, head_inputArrays, tail_inputArrays, tau_l_AMPA, osc_E_ext_id, post, is_all_delays_the_same);
		}

		// Inhibition
		if (cmp(0.0, N_input_GABA, 1e-14) < 0)
		{
			add_Input(NULL, 0, t, N_input_GABA, head_inputArrays, tail_inputArrays, tau_l_GABA, osc_I_ext_id, post, is_all_delays_the_same);
		}
	}

	if (type_of_external_inputs == PERIODIC_EXTERNAL_INPUTS)
	{
		// Excitation
		if (time_shift + n_periodic_AMPA_external_inputs_to_E/ext_input_props->lambda_AMPA_on_E <= global_t)
		{
			n_periodic_AMPA_external_inputs_to_E += 1;
		}

		if (time_shift + n_periodic_AMPA_external_inputs_to_I/ext_input_props->lambda_AMPA_on_I <= global_t)
		{
			n_periodic_AMPA_external_inputs_to_I += 1;
		}

		// Inhibition
		if (time_shift + n_periodic_GABA_external_inputs_to_E/ext_input_props->lambda_GABA_on_E <= global_t)
		{
			n_periodic_GABA_external_inputs_to_E += 1;
		}

		if (time_shift + n_periodic_GABA_external_inputs_to_I/ext_input_props->lambda_GABA_on_I <= global_t)
		{
			n_periodic_GABA_external_inputs_to_I += 1;
		}
	}
}

static void init_struct_arrays(void *A[], int N)
{
	int i;
	for (i=0; i<N; i++)
	{
		A[i] = NULL;
	}
}

static void update_factor_refrac(factor_refrac_st_struct *factor_refrac_ptr, int is_just_after_refrac[], double factor_counter[], sim_props_struct* sim_p, network_props_struct *network_p)
{
	double 	dt	=sim_p->dt;
	int 	N	=network_p->N_i + network_p->N_e;

	int i;
	for (i = 0; i < N; i++)
	{
		/* Keep the state before */
		int f_ref_bef = factor_refrac_ptr->is_not_in_refrac;

		/* Update the refractory period variables */
		factor_counter[i] = factor_counter[i] - dt;

		if (factor_counter[i] <= 0)
		{
			/* Not during refractory period */
			factor_counter[i] = 0.0;

			factor_refrac_ptr->is_not_in_refrac = 1;
		}
		else
		{
			/* During refractory period */
			factor_refrac_ptr->is_not_in_refrac = 0;
		}

		/* Check if it just after passes the refractory period */
		if ((f_ref_bef == 0) && (factor_refrac_ptr->is_not_in_refrac == 1))
		{
			is_just_after_refrac[i] = 1;
		}
		else
		{
			is_just_after_refrac[i] = 0;
		}

		factor_refrac_ptr = factor_refrac_ptr->next;
	}
}

static void init_arrays(double A[],double x,int N)
{
	int i;
	for (i=0;i<N;i++)
	{
		A[i]=x;
	}
}

static void init_arrays_int(int A[], int x, int N)
{
	int i;
	for (i = 0; i < N; i++)
	{
		A[i] = x;
	}
}

static double get_RmIe(int i, int N_i, double RmIe_E, double RmIe_I)
{
	if (i < N_i)
	{
		return RmIe_I;
	}
	else
	{
		return RmIe_E;
	}
}

static double get_Rm(int i, int N_i, double Rm_E, double Rm_I)
{
	if (i < N_i)
	{
		return Rm_I;
	}
	else
	{
		return Rm_E;
	}
}


static int chk_is_all_delays_the_same(double *tau_l_AMPA_on_E, double *tau_l_AMPA_on_I, double *tau_l_GABA_on_E, double *tau_l_GABA_on_I, double *tau_l_ext_AMPA_on_E, double *tau_l_ext_AMPA_on_I, double *tau_l_ext_GABA_on_E, double *tau_l_ext_GABA_on_I)
{
	double tol_eq = 1e-10;

	if ((cmp(*tau_l_AMPA_on_E, *tau_l_AMPA_on_I, tol_eq) == 0)
		&& (cmp(*tau_l_AMPA_on_E, *tau_l_GABA_on_E, tol_eq) == 0)
		&& (cmp(*tau_l_AMPA_on_E, *tau_l_GABA_on_I, tol_eq) == 0)
		&& (cmp(*tau_l_AMPA_on_E, *tau_l_ext_AMPA_on_E, tol_eq) == 0)
		&& (cmp(*tau_l_AMPA_on_E, *tau_l_ext_AMPA_on_I, tol_eq) == 0)
		&& (cmp(*tau_l_AMPA_on_E, *tau_l_ext_GABA_on_E, tol_eq) == 0)
		&& (cmp(*tau_l_AMPA_on_E, *tau_l_ext_GABA_on_I, tol_eq) == 0)
		&& (is_dist_EE_latency == 0)
		&& (is_dist_EI_latency == 0)
		&& (is_dist_IE_latency == 0)
		&& (is_dist_II_latency == 0))
	{
		*tau_l_AMPA_on_I = *tau_l_AMPA_on_E;
		*tau_l_GABA_on_E = *tau_l_AMPA_on_E;
		*tau_l_GABA_on_I = *tau_l_AMPA_on_E;
		*tau_l_ext_AMPA_on_E = *tau_l_AMPA_on_E;
		*tau_l_ext_AMPA_on_I = *tau_l_AMPA_on_E;
		*tau_l_ext_GABA_on_E = *tau_l_AMPA_on_E;
		*tau_l_ext_GABA_on_I = *tau_l_AMPA_on_E;

		printf("[%s] all delays are the same\n", getCurrentTime());fflush(stdout);
		return 1;
	}
	else
	{
		printf("[%s] all delays are NOT the same or delays are random\n", getCurrentTime());fflush(stdout);
		return 0;
	}
}

static double get_R_m(int id_osc, int N_i)
{
	if (id_osc < N_i)
	{
		return R_m_i;
	}
	else
	{
		return R_m_e;
	}
}

static double get_tau_m(int id_osc, int N_i)
{
	if (id_osc < N_i)
	{
		return tau_m_i;
	}
	else
	{
		return tau_m_e;
	}
}

/**************** GSL related functions ****************/
static int RHS_v_v_perts(double t, const double V_dV[], double f_V_dV[], void *params)
{
	int i;
	int N_i = ((ode_solver_params_struct *)params)->network_props->N_i;
	int N_e = ((ode_solver_params_struct *)params)->network_props->N_e;
	int	is_initV_mode = ((ode_solver_params_struct *)params)->is_initV_mode;

	factor_refrac_st_struct *factor_refrac_ptr 	= ((ode_solver_params_struct *)params)->factor_refrac_ptr;
	S_elem_struct 			*head_S 			= ((ode_solver_params_struct *)params)->head_S;
	X_elem_struct 			*head_X 			= ((ode_solver_params_struct *)params)->head_X;
	syn_props_struct 		*syn_props 			= ((ode_solver_params_struct *)params)->syn_props;
	ext_input_props_struct 	*ext_input_props 	= ((ode_solver_params_struct *)params)->ext_input_props;
	post_syn_osc_struct		**pre_syn_Mat_ptr	= ((ode_solver_params_struct *)params)->pre_syn_Mat_ptr;

	double 					*RmIe_arrays		= ((ode_solver_params_struct *)params)->RmIe_arrays_ptr;


	for (i = 0; i < N_i + N_e; i++)
	{
		double R_m 			= get_R_m		(i, N_i);
		double RmIe 		= getDoublePtrEle(RmIe_arrays, 0, i, 1, N_i + N_e);

		double tau_m_neuron = get_tau_m		(i, N_i);

		post_syn_osc_struct	*pre_syn_list = getPre_syn_MatPtrEle(pre_syn_Mat_ptr, 0, i, 1, N_i + N_e);

		double g_syn_E, g_syn_I, g_syn_E_ext, g_syn_I_ext;
		double V_syn_E, V_syn_I, V_syn_E_ext, V_syn_I_ext;
		double tau_m_E, tau_m_I, tau_m_E_ext, tau_m_I_ext;
		double tau_r_E, tau_r_I, tau_r_E_ext, tau_r_I_ext;
		double tau_d_E, tau_d_I, tau_d_E_ext, tau_d_I_ext;
		get_synapse_kinetics(i, N_i,
				&g_syn_E, &g_syn_I, &g_syn_E_ext, &g_syn_I_ext,
				&V_syn_E, &V_syn_I, &V_syn_E_ext, &V_syn_I_ext,
				&tau_m_E, &tau_m_I, &tau_m_E_ext, &tau_m_I_ext,
				&tau_r_E, &tau_r_I, &tau_r_E_ext, &tau_r_I_ext,
				&tau_d_E, &tau_d_I, &tau_d_E_ext, &tau_d_I_ext,
				syn_props, ext_input_props);

		double s_syn_E 		= cal_s_syn(t, t_before_call_gsl, tau_r_E, 		tau_d_E, 		head_S->s_E, 		head_X->x_E);
		double s_syn_I 		= cal_s_syn(t, t_before_call_gsl, tau_r_I, 		tau_d_I, 		head_S->s_I, 		head_X->x_I);
		double s_syn_E_ext 	= cal_s_syn(t, t_before_call_gsl, tau_r_E_ext,	tau_d_E_ext,	head_S->s_E_ext,	head_X->x_E_ext);
		double s_syn_I_ext 	= cal_s_syn(t, t_before_call_gsl, tau_r_I_ext,	tau_d_I_ext,	head_S->s_I_ext,	head_X->x_I_ext);

		/********************** dVdt **********************/
		double V = V_dV[2*i];
		double I_syn_E, I_syn_I;
		double I_syn_E_ext, I_syn_I_ext;

		if (is_Cond_synapse == 1)
		{
			I_syn_E 		= R_m*g_syn_E		*s_syn_E	*(V - V_syn_E);
			I_syn_I 		= R_m*g_syn_I		*s_syn_I	*(V - V_syn_I);
			I_syn_E_ext 	= R_m*g_syn_E_ext	*s_syn_E_ext*(V - V_syn_E_ext);
			I_syn_I_ext 	= R_m*g_syn_I_ext	*s_syn_I_ext*(V - V_syn_I_ext);
		}
		else
		{
			I_syn_E 		= R_m*g_syn_E		*s_syn_E	*(V_i_th - V_syn_E);
			I_syn_I 		= R_m*g_syn_I		*s_syn_I	*(V_i_th - V_syn_I);
			I_syn_E_ext 	= R_m*g_syn_E_ext	*s_syn_E_ext*(V_i_th - V_syn_E_ext);
			I_syn_I_ext 	= R_m*g_syn_I_ext	*s_syn_I_ext*(V_i_th - V_syn_I_ext);
		}

		double RmI_period_ext_drive;
		if (is_Turn_On_Cont_Sin_Ext_I == 1)
		{
			double phi_i;
			if (is_All_phi_the_same == 1)
			{
				phi_i = 0.0;
			}
			else
			{
				phi_i = 2*PI*i/(N_i + N_e);
			}

			RmI_period_ext_drive = tau_m_neuron*A_cos_ext_input*cos(2*PI*f_cos_ext_input*t + phi_i);
		}
		else
		{
			RmI_period_ext_drive = 0.0;
		}

		if (i < N_i)
		{
			// Interneurons
			switch(INTN_TYPE)
			{
			case WANG_BUZSAKI:
			{
				double wb_h = V_dV[2*(N_i + N_e) + wb_N_states*i + 0];
				double wb_n = V_dV[2*(N_i + N_e) + wb_N_states*i + 1];

#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
				double I_Na = wb_get_I_Na(w_g_Na, V, wb_h);
				double I_K 	= wb_get_I_K(w_g_K, V, wb_n);
				double I_L 	= wb_get_I_L(w_g_L, V);
#else
				double I_Na = wb_get_I_Na(wb_g_Na, V, wb_h);
				double I_K 	= wb_get_I_K(wb_g_K, V, wb_n);
				double I_L 	= wb_get_I_L(wb_g_L, V);
#endif

				double I_gap;
				if (is_initV_mode == 1)
				{
					I_gap = 0.0;
				}
				else
				{
					I_gap = wb_get_I_gap(V, V_dV, pre_syn_list);
				}

#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
				f_V_dV[2*i] = (-I_gap -I_Na - I_K - I_L - w_g_L*I_syn_E - w_g_L*I_syn_I - w_g_L*I_syn_E_ext - w_g_L*I_syn_I_ext + w_g_L*RmIe + w_g_L*RmI_period_ext_drive)/w_C_m;
				f_V_dV[2*(N_i + N_e) + wb_N_states*i + 0] = wb_get_f_h(V, wb_h);
				f_V_dV[2*(N_i + N_e) + wb_N_states*i + 1] = wb_get_f_n(V, wb_n);
#else
				f_V_dV[2*i] = (-I_gap -I_Na - I_K - I_L - wb_g_L*I_syn_E - wb_g_L*I_syn_I - wb_g_L*I_syn_E_ext - wb_g_L*I_syn_I_ext + wb_g_L*RmIe + wb_g_L*RmI_period_ext_drive)/wb_C_m;
				f_V_dV[2*(N_i + N_e) + wb_N_states*i + 0] = wb_get_f_h(V, wb_h);
				f_V_dV[2*(N_i + N_e) + wb_N_states*i + 1] = wb_get_f_n(V, wb_n);
#endif

				/********************** dV_pertsdt **********************/
				f_V_dV[2*i + 1] = 0.0;

				break;
			}
			case BORGER_WALKER:
			{
				double bw_h = V_dV[2*(N_i + N_e) + bw_N_states*i + 0];
				double bw_n = V_dV[2*(N_i + N_e) + bw_N_states*i + 1];

				double I_Na = bw_get_I_Na(bw_g_Na, V, bw_h);
				double I_K 	= bw_get_I_K(bw_g_K, V, bw_n);
				double I_L 	= bw_get_I_L(bw_g_L, V);

				double I_gap;
				if (is_initV_mode == 1)
				{
					I_gap = 0.0;
				}
				else
				{
					I_gap = wb_get_I_gap(V, V_dV, pre_syn_list);
				}

				f_V_dV[2*i] = (-I_gap -I_Na - I_K - I_L - bw_g_L*I_syn_E - bw_g_L*I_syn_I - bw_g_L*I_syn_E_ext - bw_g_L*I_syn_I_ext + bw_g_L*RmIe + bw_g_L*RmI_period_ext_drive)/bw_C_m;
				f_V_dV[2*(N_i + N_e) + bw_N_states*i + 0] = bw_f_h(V, bw_h);
				f_V_dV[2*(N_i + N_e) + bw_N_states*i + 1] = bw_f_n(V, bw_n);

				/********************** dV_pertsdt **********************/
				f_V_dV[2*i + 1] = 0.0;

				break;
			}
			}
		}
		else
		{
			// Pyramidal cells
			double h_NaT = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 0];
			double m_CaT = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 1];
			double h_CaT = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 2];
			double m_CaH = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 3];
			double h_CaH = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 4];
			double m_KDR = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 5];
			double h_KDR = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 6];
			double m_KM  = V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 7];

			double m_NaT = nw_cal_x_inf(V, -37, 5);
			double m_NaP = nw_cal_x_inf(V, -47, 3);

			double tau_h_NaT = 0.2 + 0.007*exp(exp(-(V - 40.6)/51.4));

			double I_NaT = nw_g_NaT*(m_NaT*m_NaT*m_NaT)*h_NaT*(V - nw_E_Na);
			double I_NaP = nw_g_NaP*m_NaP*(V - nw_E_Na);
			double I_CaT = nw_g_CaT*(m_CaT*m_CaT)*h_CaT*(V - nw_E_Ca);
			double I_CaH = nw_g_CaH*(m_CaH*m_CaH)*h_CaH*(V - nw_E_Ca);
			double I_KDR = nw_g_KDR*m_KDR*h_KDR*(V - nw_E_K);
			double I_KM  = nw_g_KM*m_KM*(V - nw_E_K);
			double I_L   = nw_g_L*(V - nw_E_L);

			f_V_dV[2*i] = (-I_NaT -I_NaP -I_CaT -I_CaH -I_KDR -I_KM -I_L - nw_g_L*I_syn_E - nw_g_L*I_syn_I - nw_g_L*I_syn_E_ext - nw_g_L*I_syn_I_ext + nw_g_L*RmIe + nw_g_L*RmI_period_ext_drive)/nw_C_m;

			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 0]	= nw_rhs_gating_variable(V, -75.0, -7.0, h_NaT, tau_h_NaT);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 1]	= nw_rhs_gating_variable(V, -54.0, 5.0, m_CaT, 2.0);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 2] = nw_rhs_gating_variable(V, -65.0, -8.5, h_CaT, 32.0);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 3] = nw_rhs_gating_variable(V, -15.0, 5.0, m_CaH, 0.08);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 4] = nw_rhs_gating_variable(V, -60.0, -7.0, h_CaH, 300.0);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 5] = nw_rhs_gating_variable(V, -5.8, 11.4, m_KDR, 1.0);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 6] = nw_rhs_gating_variable(V, -68.0, -9.7, h_KDR, 1400.0);
			f_V_dV[2*(N_i + N_e) + wb_N_states*N_i + nw_N_states*(i - N_i) + 7]	= nw_rhs_gating_variable(V, -30.0, 10.0, m_KM, 75.0);

			/********************** dV_pertsdt **********************/
			f_V_dV[2*i + 1] = 0.0;
		}

		factor_refrac_ptr = factor_refrac_ptr->next;
		head_S = head_S->next;
		head_X = head_X->next;
	}

	return GSL_SUCCESS;
}

static void get_synapse_kinetics(
		int osc_id, int N_i,
		double *g_syn_E, double *g_syn_I, double *g_syn_E_ext, double *g_syn_I_ext,
		double *V_syn_E, double *V_syn_I, double *V_syn_E_ext, double *V_syn_I_ext,
		double *tau_m_E, double *tau_m_I, double *tau_m_E_ext, double *tau_m_I_ext,
		double *tau_r_E, double *tau_r_I, double *tau_r_E_ext, double *tau_r_I_ext,
		double *tau_d_E, double *tau_d_I, double *tau_d_E_ext, double *tau_d_I_ext,
		syn_props_struct *syn_props, ext_input_props_struct *ext_input_props)
{
	if (osc_id < N_i)
	{
		*g_syn_E		= syn_props->g_syn_AMPA_on_I;
		*g_syn_I		= syn_props->g_syn_GABA_on_I;
		*g_syn_E_ext	= ext_input_props->g_ext_AMPA_on_I;
		*g_syn_I_ext	= ext_input_props->g_ext_GABA_on_I;

		*V_syn_E		= syn_props->v_syn_AMPA;
		*V_syn_I		= syn_props->v_syn_GABA;
		*V_syn_E_ext	= ext_input_props->v_syn_AMPA;
		*V_syn_I_ext	= ext_input_props->v_syn_GABA;

		*tau_m_E		= syn_props->tau_m_on_I;
		*tau_m_I		= syn_props->tau_m_on_I;
		*tau_m_E_ext	= ext_input_props->tau_m_on_I;
		*tau_m_I_ext	= ext_input_props->tau_m_on_I;

		*tau_r_E		= syn_props->tau_r_AMPA_on_I;
		*tau_r_I		= syn_props->tau_r_GABA_on_I;
		*tau_r_E_ext	= ext_input_props->tau_r_AMPA_on_I;
		*tau_r_I_ext	= ext_input_props->tau_r_GABA_on_I;

		*tau_d_E		= syn_props->tau_d_AMPA_on_I;
		*tau_d_I		= syn_props->tau_d_GABA_on_I;
		*tau_d_E_ext	= ext_input_props->tau_d_AMPA_on_I;
		*tau_d_I_ext	= ext_input_props->tau_d_GABA_on_I;
	}
	else
	{
		*g_syn_E		= syn_props->g_syn_AMPA_on_E;
		*g_syn_I		= syn_props->g_syn_GABA_on_E;
		*g_syn_E_ext	= ext_input_props->g_ext_AMPA_on_E;
		*g_syn_I_ext	= ext_input_props->g_ext_GABA_on_E;

		*V_syn_E		= syn_props->v_syn_AMPA;
		*V_syn_I		= syn_props->v_syn_GABA;
		*V_syn_E_ext	= ext_input_props->v_syn_AMPA;
		*V_syn_I_ext	= ext_input_props->v_syn_GABA;

		*tau_m_E		= syn_props->tau_m_on_E;
		*tau_m_I		= syn_props->tau_m_on_E;
		*tau_m_E_ext	= ext_input_props->tau_m_on_E;
		*tau_m_I_ext	= ext_input_props->tau_m_on_E;

		*tau_r_E		= syn_props->tau_r_AMPA_on_E;
		*tau_r_I		= syn_props->tau_r_GABA_on_E;
		*tau_r_E_ext	= ext_input_props->tau_r_AMPA_on_E;
		*tau_r_I_ext	= ext_input_props->tau_r_GABA_on_E;

		*tau_d_E		= syn_props->tau_d_AMPA_on_E;
		*tau_d_I		= syn_props->tau_d_GABA_on_E;
		*tau_d_E_ext	= ext_input_props->tau_d_AMPA_on_E;
		*tau_d_I_ext	= ext_input_props->tau_d_GABA_on_E;
	}
}


static double cal_s_syn(double t, double t0, double tau_r, double tau_d, double s0, double x0)
{
	return M_EXP(-(t - t0)/tau_d)*s0 + (tau_r*x0)/(tau_r - tau_d)*(M_EXP(-(t - t0)/tau_r) - M_EXP(-(t - t0)/tau_d));
}

static void update_S_X(double dt, int N_i, int N_osc, S_elem_struct *head_S, X_elem_struct *head_X, syn_props_struct *syn_props, ext_input_props_struct *ext_input_props)
{
	int i;
	for (i = 0; i < N_osc; i++)
	{
		double dummy;
		double tau_r_E, tau_r_I, tau_r_E_ext, tau_r_I_ext;
		double tau_d_E, tau_d_I, tau_d_E_ext, tau_d_I_ext;
		get_synapse_kinetics(i, N_i,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				&dummy, &dummy, &dummy, &dummy,
				&tau_r_E, &tau_r_I, &tau_r_E_ext, &tau_r_I_ext,
				&tau_d_E, &tau_d_I, &tau_d_E_ext, &tau_d_I_ext,
				syn_props, ext_input_props);

		head_S->s_E		= cal_s_syn(dt, 0.0, tau_r_E, 		tau_d_E, 	head_S->s_E, 	head_X->x_E);
		head_S->s_I		= cal_s_syn(dt, 0.0, tau_r_I, 		tau_d_I, 	head_S->s_I, 	head_X->x_I);
		head_S->s_E_ext	= cal_s_syn(dt, 0.0, tau_r_E_ext,	tau_d_E_ext,head_S->s_E_ext,head_X->x_E_ext);
		head_S->s_I_ext	= cal_s_syn(dt, 0.0, tau_r_I_ext,	tau_d_I_ext,head_S->s_I_ext,head_X->x_I_ext);

		head_X->x_E 	= M_EXP(-dt/tau_r_E)		*head_X->x_E;
		head_X->x_I 	= M_EXP(-dt/tau_r_I)		*head_X->x_I;
		head_X->x_E_ext = M_EXP(-dt/tau_r_E_ext)	*head_X->x_E_ext;
		head_X->x_I_ext = M_EXP(-dt/tau_r_I_ext)	*head_X->x_I_ext;

		head_S = head_S->next;
		head_X = head_X->next;
	}
}

static void init_SX_link_list(
		int N_osc,
		factor_refrac_st_struct **head_factor_refrac_ptr, factor_refrac_st_struct **tmp_factor_refrac_ptr,
		S_elem_struct **head_S, S_elem_struct **tmp_S,
		X_elem_struct **head_X, X_elem_struct **tmp_X)
{
	int i;
	for (i = 0; i < N_osc; i++)
	{
		/* Initialize factor refractory link list */
		factor_refrac_st_struct *tmp = (factor_refrac_st_struct *)malloc(sizeof(factor_refrac_st_struct));
		tmp->is_not_in_refrac = 1;
		tmp->next = NULL;

		if (*head_factor_refrac_ptr == NULL)
		{
			(*head_factor_refrac_ptr) 	= tmp;
			(*tmp_factor_refrac_ptr) 	= tmp;
		}
		else
		{
			(*tmp_factor_refrac_ptr)->next = tmp;
			(*tmp_factor_refrac_ptr) = tmp;
		}

		/* Initialize S link list */
		S_elem_struct *tmp_tmp = (S_elem_struct *)malloc(sizeof(S_elem_struct));
		tmp_tmp->s_E = 0.0;
		tmp_tmp->s_I = 0.0;
		tmp_tmp->s_E_ext = 0.0;
		tmp_tmp->s_I_ext = 0.0;
		tmp_tmp->next = NULL;

		if (*head_S == NULL)
		{
			(*head_S) = tmp_tmp;
			(*tmp_S) = tmp_tmp;
		}
		else
		{
			(*tmp_S)->next = tmp_tmp;
			(*tmp_S) = tmp_tmp;
		}

		/* Initialize X link list */
		X_elem_struct *tmp_tmp_tmp = (X_elem_struct *)malloc(sizeof(X_elem_struct));
		tmp_tmp_tmp->x_E = 0.0;
		tmp_tmp_tmp->x_I = 0.0;
		tmp_tmp_tmp->x_E_ext = 0.0;
		tmp_tmp_tmp->x_I_ext = 0.0;
		tmp_tmp_tmp->next = NULL;

		if (*head_X == NULL)
		{
			(*head_X) = tmp_tmp_tmp;
			(*tmp_X) = tmp_tmp_tmp;
		}
		else
		{
			(*tmp_X)->next = tmp_tmp_tmp;
			(*tmp_X) = tmp_tmp_tmp;
		}
	}
}

static void print_meta_files(char *self_FN)
{
	struct stat fileStat;

	if(stat(self_FN, &fileStat) < 0)
	{
		printf("[print_meta_files] Didn't success to open [%s]", self_FN);
		return;
	}

	printf("%s\n", self_FN);
	printf("Last status change:       %s", ctime(&fileStat.st_ctime));
	printf("Last file access:         %s", ctime(&fileStat.st_atime));
	printf("Last file modification:   %s", ctime(&fileStat.st_mtime));
}


static void openFiles(
		char *abs_file_name_Str, 	int is_rec_spkp, 		FILE **pFileSpkProfile,
		char *abs_file_name_Str1, 	int is_rec_traj, 		FILE **pFileSpkProfile1,
		char *abs_file_name_Str2, 	int is_rec_LEs, 		FILE **pFileSpkProfile2,
		char *abs_file_name_Str3, 	int is_rec_LEs_End_Pert,FILE **pFileSpkProfile3,
		char *abs_file_name_Str4, 	int is_rec_LEs_Pert, 	FILE **pFileSpkProfile4,
		char *abs_file_name_Str5, 	int is_rec_LEs_End, 	FILE **pFileSpkProfile5)
{
	/* Open the files */
	if (is_rec_LEs_End == 1)
	{
		*pFileSpkProfile5		=fopen(abs_file_name_Str5, 	"wb");
		if (*pFileSpkProfile5 == NULL)
		{
			printf("Cannot open [%s] and exit(EXIT_FAILURE)\n", abs_file_name_Str5);
			exit(EXIT_FAILURE);
		}
	}

	if (is_rec_spkp == 1)
	{
		*pFileSpkProfile		=fopen(abs_file_name_Str, 	"wb");
		if (*pFileSpkProfile == NULL)
		{
			printf("Cannot open [%s] and exit(EXIT_FAILURE)\n", abs_file_name_Str);
			exit(EXIT_FAILURE);
		}
	}

	if (is_rec_traj == 1)
	{
		*pFileSpkProfile1 	=fopen(abs_file_name_Str1, 	"wb");
		if (*pFileSpkProfile1 == NULL)
		{
			printf("Cannot open [%s] and exit(EXIT_FAILURE)\n", abs_file_name_Str1);
			exit(EXIT_FAILURE);
		}
	}

	if (is_rec_LEs == 1)
	{
		*pFileSpkProfile2	=fopen(abs_file_name_Str2, 	"wb");
		if (*pFileSpkProfile2 == NULL)
		{
			printf("Cannot open [%s] and exit(EXIT_FAILURE)\n", abs_file_name_Str2);
			exit(EXIT_FAILURE);
		}
	}

	if (is_rec_LEs_End_Pert == 1)
	{
		*pFileSpkProfile3	=fopen(abs_file_name_Str3, 	"wb");
		if (*pFileSpkProfile3 == NULL)
		{
			printf("Cannot open [%s] and exit(EXIT_FAILURE)\n", abs_file_name_Str3);
			exit(EXIT_FAILURE);
		}
	}

	if (is_rec_LEs_Pert == 1)
	{
		*pFileSpkProfile4	=fopen(abs_file_name_Str4, 	"wb");
		if (*pFileSpkProfile4 == NULL)
		{
			printf("Cannot open [%s] and exit(EXIT_FAILURE)\n", abs_file_name_Str4);
			exit(EXIT_FAILURE);
		}
	}
}

static double cal_S_kernel(double t, double tau_r, double tau_d)
{
	return (tau_r/(tau_r - tau_d))*(M_EXP(-t/tau_r) - M_EXP(-t/tau_d));
}

static void inc_X_and_update_S_X_pert_after_spike_arrival(double V[], S_elem_struct *head_S, syn_props_struct *syn_props, ext_input_props_struct
		*ext_input_props, int is_connected_with_E[], int is_connected_with_I[], factor_refrac_st_struct *head_factor_refrac_ptr, int N_LE,
		double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[], double t, int N_i, int N_e,
		X_elem_struct *head_X, input_struct *head_inputArrays[], input_struct *tail_inputArrays[], double *RmIe_arrays)
{
	int i_Osc; // i_Osc is running from 0 to represent a post-synaptic neuron.

	/* Check each post-synaptic neuron for incoming inputs. */
	for (i_Osc = 0; i_Osc < N_i + N_e; i_Osc++)
	{
		int 	N_inputs_E = 0, N_inputs_I = 0, N_inputs_EI = 0;
		double 	A_E = 0.0, A_I = 0.0;
		int 	is_not_in_refrac = head_factor_refrac_ptr->is_not_in_refrac;

		/* Store sorted inputs from E and I neurons */
		input_struct 	*ref_inputs[N_LE], *head_inputArrays_from_EI[N_LE], *tail_inputArrays_from_EI[N_LE];
		int 			N_inputs_E_bef_or_eq[N_LE], N_inputs_E_aft[N_LE];
		int 			N_inputs_I_bef_or_eq[N_LE], N_inputs_I_aft[N_LE];
		int 			N_inputs_EI_bef_or_eq[N_LE], N_inputs_EI_aft[N_LE];

		/* Initialize the storage of the sorted inputs */
		int j;
		for (j = 0; j < N_LE; j++)
		{
			head_inputArrays_from_EI[j] = NULL;
			tail_inputArrays_from_EI[j] = NULL;
			ref_inputs[j] = NULL;

			N_inputs_EI_bef_or_eq[j] = 0;	N_inputs_EI_aft[j] = 0;
			N_inputs_E_bef_or_eq[j] = 0;	N_inputs_E_aft[j] = 0;
			N_inputs_I_bef_or_eq[j] = 0;	N_inputs_I_aft[j] = 0;
		}

		/* Add a reference time */
		for (j = 0; j < N_LE; j++)
		{
			add_Input_v2(NULL, NULL, NULL, NULL, NULL, 0.0, t, 0.0, head_inputArrays_from_EI, tail_inputArrays_from_EI, 0.0, -1, j);
			ref_inputs[j] = head_inputArrays_from_EI[j];
		}

		/* Check for each osc. if there are inputs */
		input_struct *tmp	= head_inputArrays[i_Osc];
		while ((tmp != NULL) && (cmp(tmp->arr_time, t, 1e-10) <= 0))
		{
			/**/
//			if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (tmp->from < (N_i + N_e)))
			if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (tmp->from < (N_i + N_e)) && (osc_obs == i_Osc))
			{
				printf("osc. %d get inputs from osc. %d at %1.16lf\n", i_Osc, tmp->from, global_t);
			}
			/**/

			/* Get the synaptic constants */
			double tau_r, tau_m = 1.0, numberOfSpk;

			numberOfSpk	= tmp->N_input;

			if (tmp->from < N_i)
			{	// From I
				N_inputs_I++;
				N_inputs_EI++;

				if (i_Osc < N_i)
				{
					tau_r = syn_props->tau_r_GABA_on_I;// On I, i.e. I->I
					tau_m = 1.0;

					switch(syn_type_I2I)
					{
					case PULSE_SYNAPSE: // Dirac delta pulse synapse
					{
						double g_syn_E, g_syn_I, g_syn_E_ext, g_syn_I_ext;
						double V_syn_E, V_syn_I, V_syn_E_ext, V_syn_I_ext;
						double tau_m_E, tau_m_I, tau_m_E_ext, tau_m_I_ext;
						double tau_r_E, tau_r_I, tau_r_E_ext, tau_r_I_ext;
						double tau_d_E, tau_d_I, tau_d_E_ext, tau_d_I_ext;
						get_synapse_kinetics(i_Osc, N_i, &g_syn_E, &g_syn_I, &g_syn_E_ext, &g_syn_I_ext, &V_syn_E, &V_syn_I, &V_syn_E_ext, &V_syn_I_ext, &tau_m_E, &tau_m_I, &tau_m_E_ext, &tau_m_I_ext, &tau_r_E, &tau_r_I, &tau_r_E_ext, &tau_r_I_ext, &tau_d_E, &tau_d_I, &tau_d_E_ext, &tau_d_I_ext, syn_props, ext_input_props);

						double R_m = get_R_m(i_Osc, N_i);
						double g_syn = g_syn_I;
						double E_syn = V_syn_I;
						double g_L = 0.0;
						double C_m = 0.0;

						switch(INTN_TYPE)
						{
						case WANG_BUZSAKI:
						{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
							g_L = w_g_L;
							C_m = w_C_m;
#else
							g_L = wb_g_L;
							C_m = wb_C_m;
#endif
							break;
						}
						case BORGER_WALKER:
						{
							g_L = bw_g_L;
							C_m = bw_C_m;
							break;
						}
						}

						double pulse_strength = -g_L*R_m*g_syn*(V[i_Osc] - E_syn)/C_m;
						V[i_Osc] = V[i_Osc] + pulse_strength;

						break;
					}
					default: // Bi-exp synapse
					{
						A_I =  numberOfSpk*tau_m/tau_r;
						head_X->x_I = head_X->x_I + A_I;

						break;
					}
					}
				}
				else
				{
					tau_r = syn_props->tau_r_GABA_on_E;// On E, i.e. I->E
					tau_m = 1.0;

					switch(syn_type_I2E)
					{
					case PULSE_SYNAPSE: // Dirac delta pulse synapse
					{
						double g_syn_E, g_syn_I, g_syn_E_ext, g_syn_I_ext;
						double V_syn_E, V_syn_I, V_syn_E_ext, V_syn_I_ext;
						double tau_m_E, tau_m_I, tau_m_E_ext, tau_m_I_ext;
						double tau_r_E, tau_r_I, tau_r_E_ext, tau_r_I_ext;
						double tau_d_E, tau_d_I, tau_d_E_ext, tau_d_I_ext;
						get_synapse_kinetics(i_Osc, N_i, &g_syn_E, &g_syn_I, &g_syn_E_ext, &g_syn_I_ext, &V_syn_E, &V_syn_I, &V_syn_E_ext, &V_syn_I_ext, &tau_m_E, &tau_m_I, &tau_m_E_ext, &tau_m_I_ext, &tau_r_E, &tau_r_I, &tau_r_E_ext, &tau_r_I_ext, &tau_d_E, &tau_d_I, &tau_d_E_ext, &tau_d_I_ext, syn_props, ext_input_props);

						double R_m = get_R_m(i_Osc, N_i);
						double g_syn = g_syn_I;
						double E_syn = V_syn_I;
						double g_L = 0.0;
						double C_m = 0.0;

						g_L = nw_g_L;
						C_m = nw_C_m;

						double pulse_strength = -g_L*R_m*g_syn*(V[i_Osc] - E_syn)/C_m;

						V[i_Osc] = V[i_Osc] + pulse_strength;

						break;
					}
					default: // Bi-exp synapse
					{
						A_I =  numberOfSpk*tau_m/tau_r;
						head_X->x_I = head_X->x_I + A_I;

						break;
					}
					}
				}

				/* Store sorted inputs I neurons */
				int k;
				v_pert_elem_struct *head_v_pert_elem = tmp->v_pert;
				for (k = 0; k < N_LE; k++)
				{
					v_pert_elem_struct *tmp_v_p_e_s = (v_pert_elem_struct *)malloc(sizeof(v_pert_elem_struct));
					tmp_v_p_e_s->prev = NULL;
					tmp_v_p_e_s->dv	= head_v_pert_elem->dv;
					tmp_v_p_e_s->is_from_E = 0;
					tmp_v_p_e_s->next = NULL;

					double delta_k = -tmp_v_p_e_s->dv/tmp->dvdt;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (tmp->from < (N_i + N_e)) && (i_Osc == osc_obs) && (k == 0))
					{
						printf("osc. %d get inputs from osc. %d, %1.16lf, %1.16lf, %1.16lf, %1.16lf, %1.16lf, %1.16lf\n", i_Osc, tmp->from, t + delta_k, t, delta_k, tmp_v_p_e_s->dv, tmp->dvdt, A_I);
					}
					/**/

					add_Input_v2(N_inputs_I_bef_or_eq, N_inputs_I_aft, N_inputs_EI_bef_or_eq, N_inputs_EI_aft, tmp_v_p_e_s, tmp->dvdt, t, numberOfSpk, head_inputArrays_from_EI, tail_inputArrays_from_EI, delta_k, tmp->from, k);

					head_v_pert_elem = head_v_pert_elem->next;
				}
			}
			else if (tmp->from < N_i + N_e)
			{	// From E
				N_inputs_E++;
				N_inputs_EI++;

				if (i_Osc < N_i)
				{
					tau_r = syn_props->tau_r_AMPA_on_I;// On I, i.e. E->I
					tau_m = 1.0;

					switch(syn_type_E2I)
					{
					case PULSE_SYNAPSE: // Dirac delta pulse synapse
					{
						double g_syn_E, g_syn_I, g_syn_E_ext, g_syn_I_ext;
						double V_syn_E, V_syn_I, V_syn_E_ext, V_syn_I_ext;
						double tau_m_E, tau_m_I, tau_m_E_ext, tau_m_I_ext;
						double tau_r_E, tau_r_I, tau_r_E_ext, tau_r_I_ext;
						double tau_d_E, tau_d_I, tau_d_E_ext, tau_d_I_ext;
						get_synapse_kinetics(i_Osc, N_i, &g_syn_E, &g_syn_I, &g_syn_E_ext, &g_syn_I_ext, &V_syn_E, &V_syn_I, &V_syn_E_ext, &V_syn_I_ext, &tau_m_E, &tau_m_I, &tau_m_E_ext, &tau_m_I_ext, &tau_r_E, &tau_r_I, &tau_r_E_ext, &tau_r_I_ext, &tau_d_E, &tau_d_I, &tau_d_E_ext, &tau_d_I_ext, syn_props, ext_input_props);

						double R_m = get_R_m(i_Osc, N_i);
						double g_syn = g_syn_E;
						double E_syn = V_syn_E;
						double g_L = 0.0;
						double C_m = 0.0;

						switch(INTN_TYPE)
						{
						case WANG_BUZSAKI:
						{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
							g_L = w_g_L;
							C_m = w_C_m;
#else
							g_L = wb_g_L;
							C_m = wb_C_m;
#endif
							break;
						}
						case BORGER_WALKER:
						{
							g_L = bw_g_L;
							C_m = bw_C_m;
							break;
						}
						}

						double pulse_strength = -g_L*R_m*g_syn*(V[i_Osc] - E_syn)/C_m;

						if (is_E2I_suprathreshold == 1)
						{
							is_I_spike = 1;
							V[i_Osc] = E2I_suprathreshold_new_V;
						}else
						{
							V[i_Osc] = V[i_Osc] + pulse_strength;
						}

						break;
					}
					default: // Bi-exp synapse
					{
						A_E = numberOfSpk*tau_m/tau_r;
						head_X->x_E = head_X->x_E + A_E;

						break;
					}
					}
				}
				else
				{
					tau_r = syn_props->tau_r_AMPA_on_E;// On E, i.e. E->E
					tau_m = 1.0;

					switch(syn_type_E2E)
					{
					case PULSE_SYNAPSE: // Dirac delta pulse synapse
					{
						double g_syn_E, g_syn_I, g_syn_E_ext, g_syn_I_ext;
						double V_syn_E, V_syn_I, V_syn_E_ext, V_syn_I_ext;
						double tau_m_E, tau_m_I, tau_m_E_ext, tau_m_I_ext;
						double tau_r_E, tau_r_I, tau_r_E_ext, tau_r_I_ext;
						double tau_d_E, tau_d_I, tau_d_E_ext, tau_d_I_ext;
						get_synapse_kinetics(i_Osc, N_i, &g_syn_E, &g_syn_I, &g_syn_E_ext, &g_syn_I_ext, &V_syn_E, &V_syn_I, &V_syn_E_ext, &V_syn_I_ext, &tau_m_E, &tau_m_I, &tau_m_E_ext, &tau_m_I_ext, &tau_r_E, &tau_r_I, &tau_r_E_ext, &tau_r_I_ext, &tau_d_E, &tau_d_I, &tau_d_E_ext, &tau_d_I_ext, syn_props, ext_input_props);

						double R_m = get_R_m(i_Osc, N_i);
						double g_syn = g_syn_E;
						double E_syn = V_syn_E;
						double g_L = 0.0;
						double C_m = 0.0;

						g_L = nw_g_L;
						C_m = nw_C_m;

						double pulse_strength = -g_L*R_m*g_syn*(V[i_Osc] - E_syn)/C_m;

						V[i_Osc] = V[i_Osc] + pulse_strength;

						break;
					}
					default: // Bi-exp synapse
					{
						A_E = numberOfSpk*tau_m/tau_r;
						head_X->x_E = head_X->x_E + A_E;

						break;
					}
					}
				}

				/* Store sorted inputs E neurons */
				int k;
				v_pert_elem_struct *head_v_pert_elem = tmp->v_pert;
				for (k = 0; k < N_LE; k++)
				{
					v_pert_elem_struct *tmp_v_p_e_s = (v_pert_elem_struct *)malloc(sizeof(v_pert_elem_struct));
					tmp_v_p_e_s->prev = NULL;
					tmp_v_p_e_s->dv	= head_v_pert_elem->dv;
					tmp_v_p_e_s->is_from_E = 1;
					tmp_v_p_e_s->next = NULL;

					double delta_k = -tmp_v_p_e_s->dv/tmp->dvdt;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (tmp->from < (N_i + N_e)) && (i_Osc == osc_obs) && (k == 0))
					{
						printf("osc. %d get inputs from osc. %d, %1.16lf, %1.16lf, %1.16lf, %1.16lf, %1.16lf, %1.16lf\n", i_Osc, tmp->from, t + delta_k, t, delta_k, tmp_v_p_e_s->dv, tmp->dvdt, A_E);
					}
					/**/

					add_Input_v2(N_inputs_E_bef_or_eq, N_inputs_E_aft, N_inputs_EI_bef_or_eq, N_inputs_EI_aft, tmp_v_p_e_s, tmp->dvdt, t, numberOfSpk, head_inputArrays_from_EI, tail_inputArrays_from_EI, delta_k, tmp->from, k);

					head_v_pert_elem = head_v_pert_elem->next;
				}
			}
			else if (tmp->from == osc_E_ext_id)
			{	// From the E external input
				if (i_Osc < N_i)
				{
					tau_r = ext_input_props->tau_r_AMPA_on_I;// On I
					tau_m = 1.0;
				}
				else
				{
					tau_r = ext_input_props->tau_r_AMPA_on_E;// On E
					tau_m = 1.0;
				}

				head_X->x_E_ext = head_X->x_E_ext + (numberOfSpk*tau_m/tau_r);
			}
			else
			{	// From the I external input
				if (i_Osc < N_i) {
					tau_r = ext_input_props->tau_r_GABA_on_I;// On I
					tau_m = 1.0;
				}else{
					tau_r = ext_input_props->tau_r_GABA_on_E;// On E
					tau_m = 1.0;
				}

				head_X->x_I_ext = head_X->x_I_ext + (numberOfSpk*tau_m/tau_r);
			}

			/* Delete the used arrival of input event */
			input_struct *tmp1 = tmp;
			tmp = tmp->next;

			input_struct *previous = tmp1->pre;
			if (previous != NULL)
			{										// tmp1 is not the head of the list.
				previous->next = tmp;
				if (tmp != NULL)
				{									// tmp1 is not the tail of the list.
					tmp->pre = previous;
				}
				else
				{
					tail_inputArrays[i_Osc]=previous;	// tmp1 is the tail of the list.
				}
			}
			else
			{										// tmp1 is the head of the list.
				if (tmp != NULL)
				{									// The list contains more than tmp1.
					tmp->pre = NULL;
					head_inputArrays[i_Osc] = tmp;
				}
				else
				{
					head_inputArrays[i_Osc] = NULL;		// The list is empty.
					tail_inputArrays[i_Osc] = NULL;
				}
			}

			/* Free memory */
			v_pert_elem_struct *head_v_pert_elem = tmp1->v_pert;
			while (head_v_pert_elem != NULL)
			{
				v_pert_elem_struct *tmp_head = head_v_pert_elem;
				head_v_pert_elem = head_v_pert_elem->next;

				free(tmp_head);
			}

			free(tmp1);
		}

		/* Delete a reference time */
		for (j = 0; j < N_LE; j++)
		{
			input_struct *pre_t_ref = ref_inputs[j]->pre;
			input_struct *post_t_ref = ref_inputs[j]->next;

			if (pre_t_ref != NULL)
			{
				pre_t_ref->next = post_t_ref;
			}

			if (post_t_ref != NULL)
			{
				post_t_ref->pre = pre_t_ref;
			}

			free(ref_inputs[j]);
			ref_inputs[j] = NULL;

			if ((pre_t_ref == NULL) && (post_t_ref == NULL))
			{// Only t reference time
				head_inputArrays_from_EI[j] = NULL;
				tail_inputArrays_from_EI[j] = NULL;
			}
			else if (pre_t_ref == NULL)
			{// t reference time is the head
				head_inputArrays_from_EI[j] = post_t_ref;
			}
			else if (post_t_ref == NULL)
			{// t reference time is the tail
				tail_inputArrays_from_EI[j] = pre_t_ref;
			}
		}

		/* Update the perturbation vectors V, S, and X */
		if (0 < N_inputs_EI)
		{
			double tau_m_I = 0.0, tau_r_I = 0.0, tau_d_I = 0.0;
			double tau_m_E = 0.0, tau_r_E = 0.0, tau_d_E = 0.0;

			if (i_Osc < N_i)
			{
				tau_r_I = syn_props->tau_r_GABA_on_I;// On I
				tau_d_I = syn_props->tau_d_GABA_on_I;
				tau_m_I = syn_props->tau_m_on_I;
			}
			else
			{
				tau_r_I = syn_props->tau_r_GABA_on_E;// On E
				tau_d_I = syn_props->tau_d_GABA_on_E;
				tau_m_I = syn_props->tau_m_on_E;
			}

			if (i_Osc < N_i)
			{
				tau_r_E = syn_props->tau_r_AMPA_on_I;// On I
				tau_d_E = syn_props->tau_d_AMPA_on_I;
				tau_m_E = syn_props->tau_m_on_I;
			}
			else
			{
				tau_r_E = syn_props->tau_r_AMPA_on_E;// On E
				tau_d_E = syn_props->tau_d_AMPA_on_E;
				tau_m_E = syn_props->tau_m_on_E;
			}

			double Cv = 0.0, Ce = 0.0, Ci = 0.0;
			cal_Cv_Ce_Ci(i_Osc, N_i, V[i_Osc], head_S->s_E, head_S->s_I, head_S->s_E_ext, head_S->s_I_ext, syn_props, ext_input_props, &Cv, &Ce, &Ci);

			update_V_perts(
					is_connected_with_E[i_Osc], is_connected_with_I[i_Osc],
					Cv, Ce, Ci,
					is_not_in_refrac,
					tau_m_E, tau_r_E, tau_d_E, tau_m_I, tau_r_I, tau_d_I,
					t, N_LE, N_i + N_e, i_Osc,
					N_inputs_EI, N_inputs_E, N_inputs_I,
					A_E, A_I,
					v_perts, s_E_perts, x_E_perts, s_I_perts, x_I_perts,
					head_inputArrays_from_EI,
					N_inputs_EI_bef_or_eq, N_inputs_EI_aft,
					N_inputs_E_bef_or_eq, N_inputs_E_aft,
					N_inputs_I_bef_or_eq, N_inputs_I_aft);
		}

		/* Free the storage of the sorted inputs */
		for (j = 0; j < N_LE; j++)
		{
			input_struct *tmp = head_inputArrays_from_EI[j];
			while (tmp != NULL)
			{
				v_pert_elem_struct *tmp_v_p = tmp->v_pert;
				while (tmp_v_p != NULL)
				{
					v_pert_elem_struct *tmp_v_p1 = tmp_v_p;

					tmp_v_p = tmp_v_p->next;
					free(tmp_v_p1);
				}

				input_struct *tmp1 = tmp;
				tmp = tmp->next;
				free(tmp1);
			}

			head_inputArrays_from_EI[j] = NULL;
			tail_inputArrays_from_EI[j] = NULL;
		}

		head_S = head_S->next;
		head_factor_refrac_ptr = head_factor_refrac_ptr->next;
		head_X = head_X->next;
	}
}

static void update_V_perts(
		int is_connected_with_E, int is_connected_with_I,
		double Cv, double Ce, double Ci,
		int is_not_in_refrac,
		double tau_m_E, double tau_r_E, double tau_d_E, double tau_m_I, double tau_r_I, double tau_d_I,
		double t, int N_LE, int N_osc, int i_Osc,
		int N_inputs_EI, int N_inputs_E, int N_inputs_I,
		double A_E, double A_I,
		double v_perts[], double s_E_perts[], double x_E_perts[], double s_I_perts[], double x_I_perts[],
		input_struct *head_inputArrays_from_EI[],
		int N_inputs_EI_bef_or_eq[], int N_inputs_EI_aft[],
		int N_inputs_E_bef_or_eq[], int N_inputs_E_aft[],
		int N_inputs_I_bef_or_eq[], int N_inputs_I_aft[])
{
	int i_LEs;
	for (i_LEs = 0; i_LEs < N_LE; i_LEs++)
	{
		double dV 	= getAEle(v_perts, i_LEs, i_Osc, N_LE, N_osc);
		double dS_E = getAEle(s_E_perts, i_LEs, i_Osc, N_LE, N_osc);
		double dX_E = getAEle(x_E_perts, i_LEs, i_Osc, N_LE, N_osc);
		double dS_I = getAEle(s_I_perts, i_LEs, i_Osc, N_LE, N_osc);
		double dX_I = getAEle(x_I_perts, i_LEs, i_Osc, N_LE, N_osc);

		/**/
		if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
		{
			printf("1:%1.16lf,%1.16lf\n", dX_E, dX_I);
			printf("1 bef:%d,%d,%d\n", N_inputs_EI_bef_or_eq[i_LEs], N_inputs_E_bef_or_eq[i_LEs], N_inputs_I_bef_or_eq[i_LEs]);
			printf("1 aft:%d,%d,%d\n", N_inputs_EI_aft[i_LEs], N_inputs_E_aft[i_LEs], N_inputs_I_aft[i_LEs]);
		}
		/**/

		double t_last_EI_upd = -1.0, t_last_E_upd = -1.0, t_last_I_upd = -1.0;

		input_struct *tmp_inputArrays = head_inputArrays_from_EI[i_LEs];

		/* Update events before and equal t */
		int i_EI_bef;
		int is_the_last_input_E;
		for (i_EI_bef = 0; i_EI_bef < N_inputs_EI_bef_or_eq[i_LEs]; i_EI_bef++)
		{
			double curr_time = tmp_inputArrays->arr_time;

			/**/
			if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
			{
				printf("begin bef:%1.16lf,%d\n", curr_time, tmp_inputArrays->v_pert->is_from_E);
			}
			/**/

			if (t_last_EI_upd < 0)
			{
				t_last_EI_upd = curr_time;
			}
			else
			{
				if (is_not_in_refrac == 1)
				{
					double dt = curr_time - t_last_EI_upd;

					if (is_the_last_input_E == 1)
					{
						dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E + A_E, dS_I, dX_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
					}
					else
					{
						dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E, dS_I, dX_I + A_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
					}
				}

				t_last_EI_upd = curr_time;
			}

			if (tmp_inputArrays->v_pert->is_from_E == 1)
			{
				if (t_last_E_upd < 0)
				{
					t_last_E_upd = curr_time;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("11:%1.16lf\n", curr_time);
					}
					/**/
				}
				else
				{
					double dt = curr_time - t_last_E_upd;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("1 bef:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					dS_E = update_dS(dS_E, dX_E, dt, tau_r_E, tau_d_E, A_E);
					dX_E = update_dX(dX_E, dt, tau_r_E, A_E);

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("1 aft:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					t_last_E_upd = curr_time;
				}

				if (0 < t_last_I_upd)
				{
					double dt = curr_time - t_last_I_upd;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("1-bef:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					dS_I = update_dS(dS_I, dX_I, dt, tau_r_I, tau_d_I, 0.0);
					dX_I = update_dX(dX_I, dt, tau_r_I, 0.0);

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("1-aft:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					t_last_I_upd = curr_time;
				}

				is_the_last_input_E = 1;
			}
			else
			{
				if (t_last_I_upd < 0)
				{
					t_last_I_upd = curr_time;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("11:%1.16lf\n", curr_time);
					}
					/**/
				}
				else
				{
					double dt = curr_time - t_last_I_upd;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("11-bef:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					dS_I = update_dS(dS_I, dX_I, dt, tau_r_I, tau_d_I, A_I);
					dX_I = update_dX(dX_I, dt, tau_r_I, A_I);

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("11-aft:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					t_last_I_upd = curr_time;
				}

				if (0 < t_last_E_upd)
				{
					double dt = curr_time - t_last_E_upd;

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("11 bef:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					dS_E = update_dS(dS_E, dX_E, dt, tau_r_E, tau_d_E, 0.0);
					dX_E = update_dX(dX_E, dt, tau_r_E, 0.0);

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("11 aft:%1.16lf,%1.16lf\n", dX_E, dX_I);
					}
					/**/

					t_last_E_upd = curr_time;
				}

				is_the_last_input_E = 0;
			}

			tmp_inputArrays = tmp_inputArrays->next;
		}

		/**/
		if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
		{
			printf("2:%1.16lf,%1.16lf\n", dX_E, dX_I);
		}
		/**/

		if (1 <= N_inputs_EI_bef_or_eq[i_LEs])
		{
			if (is_not_in_refrac == 1)
			{
				double dt = t - t_last_EI_upd;

				if (is_the_last_input_E == 1)
				{
					dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E + A_E, dS_I, dX_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
				}
				else
				{
					dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E, dS_I, dX_I + A_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
				}
			}

			t_last_EI_upd = t;
		}

		if (1 <= N_inputs_E_bef_or_eq[i_LEs])
		{
			double dt = t - t_last_E_upd;

			dS_E = update_dS(dS_E, dX_E, dt, tau_r_E, tau_d_E, A_E);
			dX_E = update_dX(dX_E, dt, tau_r_E, A_E);

			/**/
			if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
			{
				printf("21:%1.16lf,%1.16lf,%1.16lf\n", dt, t, t_last_E_upd);
			}
			/**/

			t_last_E_upd = t;
		}

		if (1 <= N_inputs_I_bef_or_eq[i_LEs])
		{
			double dt = t - t_last_I_upd;

			dS_I = update_dS(dS_I, dX_I, dt, tau_r_I, tau_d_I, A_I);
			dX_I = update_dX(dX_I, dt, tau_r_I, A_I);

			/**/
			if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
			{
				printf("21 Inh:%1.16lf,%1.16lf,%1.16lf\n", dt, t, t_last_E_upd);
			}
			/**/

			t_last_I_upd = t;
		}

		/**/
		if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
		{
			printf("3:%1.16lf,%1.16lf\n", dX_E, dX_I);
		}
		/**/

		/* Update events at t */
		if (0 < N_inputs_EI)
		{
			t_last_EI_upd = t;
		}

		if (0 < N_inputs_E)
		{
			dX_E = dX_E - (N_inputs_E - 1)*A_E;
			t_last_E_upd = t;
		}
		if (0 < N_inputs_I)
		{
			dX_I = dX_I - (N_inputs_I - 1)*A_I;
			t_last_I_upd = t;
		}

		/**/
		if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
		{
			printf("4:%1.16lf,%1.16lf\n", dX_E, dX_I);
		}
		/**/

		/* Update events after t */
		int i_EI_aft, is_E_first_input = 1, is_I_first_input = 1, is_EI_first_input = 1;
		int is_E_last_input = 0, is_I_last_input = 0;
		int cnt_E = 0, cnt_I = 0;
		for (i_EI_aft = 0; i_EI_aft < N_inputs_EI_aft[i_LEs]; i_EI_aft++)
		{
			double curr_time = tmp_inputArrays->arr_time;

			/**/
			if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
			{
				printf("begin aft:%1.16lf,%d\n", curr_time, tmp_inputArrays->v_pert->is_from_E);
			}
			/**/

			if (is_not_in_refrac == 1)
			{
				if (is_EI_first_input == 1)
				{
					double dt = curr_time - t_last_EI_upd;
					double tmp_A_E = 0.0, tmp_A_I = 0.0;
					if (0 < N_inputs_E)
					{
						tmp_A_E = A_E;
					}
					if (0 < N_inputs_I)
					{
						tmp_A_I = A_I;
					}

					dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E - tmp_A_E, dS_I, dX_I - tmp_A_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
					is_EI_first_input = 0;
				}
				else
				{
					double dt = curr_time - t_last_EI_upd;

					if (is_the_last_input_E == 1)
					{
						dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E + A_E, dS_I, dX_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
					}
					else
					{
						dV = update_dV(dt, Cv, Ce, Ci, dV, dS_E, dX_E, dS_I, dX_I + A_I, tau_r_E, tau_d_E, tau_r_I, tau_d_I, is_connected_with_E, is_connected_with_I);
					}
				}
			}

			t_last_EI_upd = curr_time;

			if (tmp_inputArrays->v_pert->is_from_E == 1)
			{
				cnt_E++;
				if (cnt_E == N_inputs_E_aft[i_LEs])
				{
					is_E_last_input = 1;
				}

				if (is_E_first_input == 1)
				{
					double dt = curr_time - t_last_E_upd;

					dS_E = update_dS(dS_E, dX_E, dt, tau_r_E, tau_d_E, -A_E);
					dX_E = update_dX(dX_E, dt, tau_r_E, -A_E);

					/**/
					if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
					{
						printf("41:%1.16lf\n", dt);
					}
					/**/

					t_last_E_upd = curr_time;
					is_E_first_input = 0;

				}
				else
				{
					double dt = curr_time - t_last_E_upd;

					dS_E = update_dS(dS_E, dX_E, dt, tau_r_E, tau_d_E, A_E);
					dX_E = update_dX(dX_E, dt, tau_r_E, A_E);

					t_last_E_upd = curr_time;
				}

				if ((0 == is_I_first_input) && (is_I_last_input == 0))
				{
					double dt = curr_time - t_last_I_upd;

					dS_I = update_dS(dS_I, dX_I, dt, tau_r_I, tau_d_I, 0.0);
					dX_I = update_dX(dX_I, dt, tau_r_I, 0.0);

					t_last_I_upd = curr_time;
				}

				is_the_last_input_E = 1;
			}
			else
			{
				cnt_I++;
				if (cnt_I == N_inputs_I_aft[i_LEs])
				{
					is_I_last_input = 1;
				}

				if (is_I_first_input == 1)
				{
					double dt = curr_time - t_last_I_upd;

					dS_I = update_dS(dS_I, dX_I, dt, tau_r_I, tau_d_I, -A_I);
					dX_I = update_dX(dX_I, dt, tau_r_I, -A_I);

					t_last_I_upd = curr_time;
					is_I_first_input = 0;
				}
				else
				{
					double dt = curr_time - t_last_I_upd;

					dS_I = update_dS(dS_I, dX_I, dt, tau_r_I, tau_d_I, A_I);
					dX_I = update_dX(dX_I, dt, tau_r_I, A_I);

					t_last_I_upd = curr_time;
				}

				if ((0 == is_E_first_input) && (is_E_last_input == 0))
				{
					double dt = curr_time - t_last_E_upd;

					dS_E = update_dS(dS_E, dX_E, dt, tau_r_E, tau_d_E, 0.0);
					dX_E = update_dX(dX_E, dt, tau_r_E, 0.0);

					t_last_E_upd = curr_time;
				}

				is_the_last_input_E = 0;
			}

			tmp_inputArrays = tmp_inputArrays->next;
		}

		/**/
		if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
		{
			printf("5:%1.16lf,%1.16lf\n", dX_E, dX_I);
		}
		/**/

		if (0 == N_inputs_E_aft[i_LEs])
		{
			dX_E = dX_E - A_E;
		}
		else
		{
			dX_E = dX_E + A_E;
		}

		if (0 == N_inputs_I_aft[i_LEs])
		{
			dX_I = dX_I - A_I;
		}
		else
		{
			dX_I = dX_I + A_I;
		}

		/**/
		if ((t_obs_beg < global_t) && (global_t < t_obs_end) && (i_Osc == osc_obs) && (i_LEs == 0))
		{
			printf("6:%1.16lf,%1.16lf\n", dX_E, dX_I);
		}
		/**/

		setAEle(v_perts, dV, i_LEs, i_Osc, N_LE, N_osc);
		setAEle(s_E_perts, dS_E, i_LEs, i_Osc, N_LE, N_osc);
		setAEle(x_E_perts, dX_E, i_LEs, i_Osc, N_LE, N_osc);
		setAEle(s_I_perts, dS_I, i_LEs, i_Osc, N_LE, N_osc);
		setAEle(x_I_perts, dX_I, i_LEs, i_Osc, N_LE, N_osc);
	}
}

static double update_dV(double dt, double Cv, double Ce, double Ci, double dV, double dS_E, double dX_E, double dS_I, double dX_I, double tau_r_E, double tau_d_E, double tau_r_I, double tau_d_I, int is_connected_with_E, int is_connected_with_I)
{
	if (is_using_linearized_perturbations == 1)
	{
		return dV;
	}
	else
	{
		double integral_terms = 0.0;
		if (is_connected_with_E == 1)
		{
			double a1 = -Cv - 1/tau_d_E;
			double a2 = -Cv - 1/tau_r_E;
			double a3 = -Cv - 1/tau_d_E;
			double A = tau_r_E/(tau_r_E - tau_d_E);
			double integral_E_term = dS_E*(M_EXP(a1*dt) - 1)/a1 + dX_E*A*(M_EXP(a2*dt) - 1)/a2 - dX_E*A*(M_EXP(a3*dt) - 1)/a3;
			integral_terms += Ce*integral_E_term;;
		}

		if (is_connected_with_I == 1)
		{
			double a1 = -Cv - 1/tau_d_I;
			double a2 = -Cv - 1/tau_r_I;
			double a3 = -Cv - 1/tau_d_I;
			double A = tau_r_I/(tau_r_I - tau_d_I);
			double integral_I_term = dS_I*(M_EXP(a1*dt) - 1)/a1 + dX_I*A*(M_EXP(a2*dt) - 1)/a2 - dX_I*A*(M_EXP(a3*dt) - 1)/a3;
			integral_terms += Ci*integral_I_term;;
		}

		return M_EXP(Cv*dt)*dV + M_EXP(Cv*dt)*integral_terms;
	}
}

static double update_dS(double dS_prev, double dX_prev, double dt, double tau_r, double tau_d, double A)
{
	if (is_using_linearized_perturbations == 1)
	{
		return dS_prev + A*dt/tau_d;
	}
	else
	{
		return dS_prev*M_EXP(-dt/tau_d) + dX_prev*cal_S_kernel(dt, tau_r, tau_d) + A*cal_S_kernel(dt, tau_r, tau_d);
	}
}

static double update_dX(double dX_prev, double dt, double tau_r, double A)
{
	if (is_using_linearized_perturbations == 1)
	{
		return dX_prev + A - A*dt/tau_r;
	}
	else
	{
		return dX_prev*M_EXP(-dt/tau_r) + A*M_EXP(-dt/tau_r);
	}
}

static void cal_Cv_Ce_Ci(int i_Osc, int N_i,
		double V, double s_E, double s_I, double s_E_ext, double s_I_ext,
		syn_props_struct *syn_props, ext_input_props_struct *ext_input_props,
		double *Cv, double *Ce, double *Ci)
{
	double R_m 			= get_R_m		(i_Osc, N_i);
	double tau_m_neuron = get_tau_m		(i_Osc, N_i);

	double g_syn_E, g_syn_I, g_syn_E_ext, g_syn_I_ext;
	double V_syn_E, V_syn_I, V_syn_E_ext, V_syn_I_ext;
	double dummy;
	get_synapse_kinetics(i_Osc, N_i,
			&g_syn_E, &g_syn_I, &g_syn_E_ext, &g_syn_I_ext,
			&V_syn_E, &V_syn_I, &V_syn_E_ext, &V_syn_I_ext,
			&dummy, &dummy, &dummy, &dummy,
			&dummy, &dummy, &dummy, &dummy,
			&dummy, &dummy, &dummy, &dummy,
			syn_props, ext_input_props);

	double dI_syn_E_dV 		= R_m*g_syn_E		*s_E;
	double dI_syn_I_dV 		= R_m*g_syn_I		*s_I;
	double dI_syn_E_ext_dV 	= R_m*g_syn_E_ext	*s_E_ext;
	double dI_syn_I_ext_dV 	= R_m*g_syn_I_ext	*s_I_ext;

	*Cv = (-1 - dI_syn_E_dV - dI_syn_I_dV - dI_syn_E_ext_dV - dI_syn_I_ext_dV)/tau_m_neuron;
	*Ce = (-R_m*g_syn_E*(V - V_syn_E))/tau_m_neuron;
	*Ci = (-R_m*g_syn_I*(V - V_syn_I))/tau_m_neuron;
}

static double get_effective_latency(double mean_tau_l, double sigma_tau_l)
{
	double cand_tau_l = gsl_cdf_gaussian_Pinv(drand48(), sigma_tau_l) + mean_tau_l;

	if (cmp(cand_tau_l, 0.0, 1e-10) <= 0)
	{
		cand_tau_l = mean_tau_l + fabs(mean_tau_l - cand_tau_l);
	}

	return cand_tau_l;
}

static void init_RmIe(int N_i, int N_e, double Rm_E, double Rm_I, double RmIe_E, double RmIe_I, double RmIe_arrays[])
{
	if (is_use_rnd_init_RmIe_seed == 1)
	{
		srand48(rnd_init_RmIe_seed);
	}

	int i;
	for (i = 0; i < (N_i + N_e); i++)
	{
		if (i < N_i)
		{
			if (is_hetero_Ie_I == 0)
			{
				RmIe_arrays[i] = RmIe_I;
			}
			else
			{
				if (cmp(sigma_Ie_I, 0.0, 1e-6) == 0)
				{
					RmIe_arrays[i] = RmIe_I;
				}
				else
				{
					double mean_I = get_RmIe(i, N_i, RmIe_E, RmIe_I)/get_Rm(i, N_i, Rm_E, Rm_I);
					RmIe_arrays[i] = (gsl_cdf_gaussian_Pinv(drand48(), sigma_Ie_I) + mean_I)*get_Rm(i, N_i, Rm_E, Rm_I);

					// norminv(drand48(), 0, sigma_Ie_I) = gsl_cdf_gaussian_Pinv(drand48(), sigma_Ie_I)
				}
			}
		}
		else
		{
			if (is_hetero_Ie_E == 0)
			{
				RmIe_arrays[i] = RmIe_E;
			}
			else
			{
				if (cmp(sigma_Ie_E, 0.0, 1e-6) == 0)
				{
					RmIe_arrays[i] = RmIe_E;
				}
				else
				{
					double mean_I = get_RmIe(i, N_i, RmIe_E, RmIe_I)/get_Rm(i, N_i, Rm_E, Rm_I);
					RmIe_arrays[i] = (gsl_cdf_gaussian_Pinv(drand48(), sigma_Ie_E) + mean_I)*get_Rm(i, N_i, Rm_E, Rm_I);

					// norminv(drand48(), 0, sigma_Ie_E) = gsl_cdf_gaussian_Pinv(drand48(), sigma_Ie_E)
				}
			}
		}
	}
}

static void my_ode_solver(double t, double dt, double perts_gsl[], int N_i, int N_osc, ode_solver_params_struct *ode_solver_params)
{
	if (is_use_explicit_Euler == 1)
	{
		explicitEuler(t, dt, perts_gsl, N_i, N_osc, ode_solver_params);
	}
}


static void explicitEuler(double t, double dt, double perts_gsl[], int N_i, int N_osc, ode_solver_params_struct *ode_solver_params)
{
	int N_e = N_osc - N_i;

#if (INTN_TYPE == WANG_BUZSAKI)
	int intn_N_states = wb_N_states;
#endif
#if (INTN_TYPE == BORGER_WALKER)
	int intn_N_states = bw_N_states;
#endif

	double f_V_dV[2*N_osc + intn_N_states*N_i + nw_N_states*N_e];

	/* Compute RHS */
	RHS_v_v_perts(t, perts_gsl, f_V_dV, ode_solver_params);

	/* Explicit Euler */
	int i_osc;
	for (i_osc = 0; i_osc < N_osc; i_osc++)
	{
		// Compute effect of the White noise
		double sigma, tau_m = get_tau_m(i_osc, N_i), w;
		if (i_osc < N_i)
		{
			sigma = sigma_white_noise_I;
			w = gsl_ran_gaussian(r_white_noise_I, 1.0);
		}
		else
		{
			sigma = sigma_white_noise_E;
			w = gsl_ran_gaussian(r_white_noise_E, 1.0);
		}

		// [sigma] = mV, [tau_m] = mS
		double eff_sigma = sigma/sqrt(tau_m);

		/* std(noise) = eff_sigma*sqrt(dt); */
		double noise = eff_sigma*sqrt(dt)*w;

		// Update V with respect to dV/dt = f(V)/C_m + sigma/sqrt(tau_m)*w(t);
		perts_gsl[2*i_osc] 		= perts_gsl[2*i_osc] + dt*f_V_dV[2*i_osc] + noise;

		// Update perturbation of V
		perts_gsl[2*i_osc + 1] 	= perts_gsl[2*i_osc + 1] + dt*f_V_dV[2*i_osc + 1];
	}

	// Update additional states of interneurons
	for (i_osc = 0; i_osc < N_i; i_osc++)
	{
		int j;
		for (j = 0; j < intn_N_states; j++)
		{
			perts_gsl[2*N_osc + intn_N_states*i_osc + j] = perts_gsl[2*N_osc + intn_N_states*i_osc + j] + dt*f_V_dV[2*N_osc + intn_N_states*i_osc + j];
		}
	}

	// Update additional states of E-cells
	for (i_osc = 0; i_osc < N_e; i_osc++)
	{
		int j;
		for (j = 0; j < nw_N_states; j++)
		{
			perts_gsl[2*N_osc + intn_N_states*N_i + nw_N_states*i_osc + j] = perts_gsl[2*N_osc + intn_N_states*N_i + nw_N_states*i_osc + j] + dt*f_V_dV[2*N_osc + intn_N_states*N_i + nw_N_states*i_osc + j];
		}
	}
}

static double wb_get_f_h(double v, double h)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	double h_inf = wb_get_h_inf(v);
	double tau_h = 0.6/(1 + exp(-0.12*(v + 67)));

	return (h_inf - h)/tau_h;
#else
	double alpha_h = 0.07*(M_EXP(-(v + 58)/20));
	double beta_h = 1/(M_EXP(-0.1*(v + 28)) + 1);

	return wb_Phi*(alpha_h*(1 - h) - beta_h*h);
#endif
}

static double wb_get_f_n(double v, double n)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	double n_inf = wb_get_n_inf(v);
	double tau_n = 0.5 + 2.0/(1 + exp(0.045*(v - 50)));

	return (n_inf - n)/tau_n;
#else
	double alpha_n = -0.01*(v + 34)/(M_EXP(-0.1*(v + 34)) - 1);
	double beta_n = 0.125*M_EXP(-(v + 44)/80.0);

	return wb_Phi*(alpha_n*(1 - n) - beta_n*n);
#endif
}

static double wb_get_I_Na(double g_Na, double v, double h)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	double m_inf = wb_get_m_inf(v);

	return g_Na*m_inf*m_inf*m_inf*h*(v - w_E_Na);
#else
	double m_inf = wb_get_m_inf(v);

	return g_Na*m_inf*m_inf*m_inf*h*(v - wb_E_Na);
#endif
}

static double wb_get_m_inf(double v)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	return 1/(1.0 + exp(-0.08*(v + 26)));
#else
   	double alpha_m = -0.1*(v + 35)/(M_EXP(-0.1*(v + 35)) - 1);
   	double beta_m = 4.0*M_EXP(-(v + 60)/18);

   	return alpha_m/(alpha_m + beta_m);
#endif
}

static double wb_get_I_K(double g_K, double v, double n)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	return g_K*n*n*n*n*(v - w_E_K);
#else
	return g_K*n*n*n*n*(v - wb_E_K);
#endif
}

static double wb_get_I_L(double g_L, double v)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	return g_L*(v - w_E_L);
#else
	return g_L*(v - wb_E_L);
#endif
}

static int isHHSpk(int *state_HH, double pre_x, double x, double x_th, double x_go_down_th)
{
	int prev_state_HH = (*state_HH);

	if (((*state_HH) == CELL_NOGEN_SPK_STATE) && (x < x_th))
	{
		(*state_HH) = CELL_NOGEN_SPK_STATE;
	}
	else if (((*state_HH) == CELL_NOGEN_SPK_STATE) && (x >= x_th) && (pre_x >= x))
	{
		(*state_HH) = CELL_NOGEN_SPK_STATE;
	}
	else if (((*state_HH) == CELL_NOGEN_SPK_STATE) && (x >= x_th) && (pre_x < x))
	{
		(*state_HH) = CELL_GEN_SPK_STATE;
	}
	else if (((*state_HH) == CELL_GEN_SPK_STATE) && (x <= x_go_down_th) && (pre_x > x))
	{
		(*state_HH) = CELL_NOGEN_SPK_STATE;
	}
	else
	{
		(*state_HH) = CELL_GEN_SPK_STATE;
	}

	if (((prev_state_HH == CELL_NOGEN_SPK_STATE) && ((*state_HH) == CELL_GEN_SPK_STATE)) || (is_I_spike == 1))
	{
		is_I_spike = 0;
		return 1;
	}
	else
	{
		return 0;
	}
}

static double cal_obs_times(int i_osc, int N_i, double RmIe_arrays[])
{
	if (i_osc < N_i)
	{// Interneurons
		switch(INTN_TYPE)
		{
		case WANG_BUZSAKI:
		{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
			double Iapp = w_g_L*RmIe_arrays[i_osc];

			double tmp = drand48()*interpolate_White_T0(Iapp);
			return tmp;

			break;
#else
			double Iapp = wb_g_L*RmIe_arrays[i_osc];

			if (Iapp < wb_Iapp_th)
			{
				return -1.0;
			}
			else
			{
				return drand48()*interpolate_WB_T0(Iapp);
			}

			break;
#endif
		}
		case BORGER_WALKER:
		{
			double Iapp = bw_g_L*RmIe_arrays[i_osc];

			if (Iapp < bw_Iapp_th)
			{
				return -1.0;
			}
			else
			{
				return drand48()*interpolate_BW_T0(Iapp);
			}

			break;
		}
		}
	}
	else
	{// Pyramidal cells
		double Iapp = nw_g_L*RmIe_arrays[i_osc];


		if (type_Pyr_cells == PYR_CELL_CA1)
		{// CA1 pyr. cells
			if ((Iapp < nw_CA1_Iapp_th) || (nw_CA1_Iapp_Maxth < Iapp))
			{
				return -1.0;
			}
			else
			{
				return drand48()*interpolate_CA1_T0(Iapp);
			}
		}
		else
		{// CA3 pyr. cells
			if ((Iapp < nw_CA3_Iapp_th) || (nw_CA3_Iapp_Maxth < Iapp))
			{
				return -1.0;
			}
			else
			{
				return drand48()*interpolate_CA3_T0(Iapp);
			}
		}
	}
}

static void init_Cells(double v[], double m_wb_HH_th, double wb_Icells[], int state_wb_Icells[], double m_mp_HH_th, double mp_Ecells[], int state_mp_Ecells[], int N_i, int N_e, ode_solver_params_struct *ode_solver_params, double RmIe_arrays[])
{
	double dt = 0.05;
	double v1= -57.5;
	double v2 = -65.0;

	if (is_use_rnd_init_V == 1)
	{
		srand48(rnd_init_V);
	}

	int N_osc = N_i + N_e;
	int i_osc, total_N_states;

#if (INTN_TYPE == WANG_BUZSAKI)
	int intn_N_states = wb_N_states;
#endif
#if (INTN_TYPE == BORGER_WALKER)
	int intn_N_states = bw_N_states;
#endif

	total_N_states = 2*N_osc + intn_N_states*N_i + nw_N_states*N_e;

	double perts_gsl[total_N_states];
	double obs_times[N_osc];
	double longest_obs_times = -2.0;

	/* Generate random times to observe the states */
	for (i_osc = 0; i_osc < N_osc; i_osc++)
	{
		if (is_rand_V_init_X == 1)
		{
			obs_times[i_osc] = cal_obs_times(i_osc, N_i, RmIe_arrays);
		}
		else
		{
			/* Old implementation changed to new on 17 Sep 2015 */
//			obs_times[i_osc] = -1.0;

			/* New implementation changed to new on 17 Sep 2015 */
			if (i_osc < N_i)
			{// I-cells
				obs_times[i_osc] = OBS_TIMES_ICELLS;
			}
			else
			{// E-cells
				obs_times[i_osc] = OBS_TIMES_ECELLS;
			}
		}

		if (longest_obs_times < obs_times[i_osc])
		{
			longest_obs_times = obs_times[i_osc];
		}
	}

	/* Initialize voltages */
	for (i_osc = 0; i_osc < N_osc; i_osc++)
	{
		if (i_osc < N_i)
		{// I-cells
			switch(INTN_TYPE)
			{
			case WANG_BUZSAKI:
			{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
				double Iapp = w_g_L*RmIe_arrays[i_osc];

//				if (Iapp < w_Iapp_th)
//				{
					perts_gsl[2*i_osc] = w_steady_V;
//				}
//				else
//				{
//					perts_gsl[2*i_osc] = interpolate_White_V(Iapp);
//				}

//				if (i_osc == 0)
//				{
//					perts_gsl[2*i_osc] = v1;
//				}
//				else
//				{
//					perts_gsl[2*i_osc] = v2;
//				}

				break;
#else
				double Iapp = wb_g_L*RmIe_arrays[i_osc];

				if (Iapp < wb_Iapp_th)
				{
					perts_gsl[2*i_osc] = wb_steady_V;
				}
				else
				{
					perts_gsl[2*i_osc] = interpolate_WB_V(Iapp);
				}

				break;
#endif
			}
			case BORGER_WALKER:
			{
				double Iapp = bw_g_L*RmIe_arrays[i_osc];

				if (Iapp < bw_Iapp_th)
				{
					perts_gsl[2*i_osc] = bw_steady_V;
				}
				else
				{
					perts_gsl[2*i_osc] = interpolate_BW_V(Iapp);
				}

				break;
			}
			}
		}
		else
		{// E-cells
			double Iapp = nw_g_L*RmIe_arrays[i_osc];

			if (type_Pyr_cells == PYR_CELL_CA1)
			{// CA1 pyr. cells
				if ((Iapp < nw_CA1_Iapp_th) || (nw_CA1_Iapp_Maxth < Iapp))
				{
					perts_gsl[2*i_osc] = nw_steady_V;
				}
				else
				{
					perts_gsl[2*i_osc] = interpolate_CA1_V(Iapp);
				}
			}
			else
			{// CA3 pyr. cells
				if ((Iapp < nw_CA3_Iapp_th) || (nw_CA3_Iapp_Maxth < Iapp))
				{
					perts_gsl[2*i_osc] = nw_steady_V;
				}
				else
				{
					perts_gsl[2*i_osc] = interpolate_CA3_V(Iapp);
				}
			}
		}

		perts_gsl[2*i_osc + 1] 	= 0.0;
	}

	/* Initialize additional states for interneurons */
	for (i_osc = 0; i_osc < N_i; i_osc++)
	{
		switch(INTN_TYPE)
		{
		case WANG_BUZSAKI:
		{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
			double Iapp = w_g_L*RmIe_arrays[i_osc];

//			if (Iapp < w_Iapp_th)
//			{
				perts_gsl[2*N_osc + wb_N_states*i_osc + 0] = wb_get_h_inf(w_steady_V);
				perts_gsl[2*N_osc + wb_N_states*i_osc + 1] = wb_get_n_inf(w_steady_V);
//			}
//			else
//			{
//				perts_gsl[2*N_osc + wb_N_states*i_osc + 0] = interpolate_White_h(Iapp);
//				perts_gsl[2*N_osc + wb_N_states*i_osc + 1] = interpolate_White_n(Iapp);
//			}

//			if (i_osc == 0)
//			{
//				perts_gsl[2*N_osc + wb_N_states*i_osc + 0] = wb_get_h_inf(v1);
//				perts_gsl[2*N_osc + wb_N_states*i_osc + 1] = wb_get_n_inf(v1);
//			}
//			else
//			{
//				perts_gsl[2*N_osc + wb_N_states*i_osc + 0] = wb_get_h_inf(v2);
//				perts_gsl[2*N_osc + wb_N_states*i_osc + 1] = wb_get_n_inf(v2);
//			}

			break;
#else
			double Iapp = wb_g_L*RmIe_arrays[i_osc];

			if (Iapp < wb_Iapp_th)
			{
				perts_gsl[2*N_osc + wb_N_states*i_osc + 0] = wb_get_h_inf(wb_steady_V);
				perts_gsl[2*N_osc + wb_N_states*i_osc + 1] = wb_get_n_inf(wb_steady_V);
			}
			else
			{
				perts_gsl[2*N_osc + wb_N_states*i_osc + 0] = interpolate_WB_h(Iapp);
				perts_gsl[2*N_osc + wb_N_states*i_osc + 1] = interpolate_WB_n(Iapp);
			}

			break;
#endif
		}
		case BORGER_WALKER:
		{
			double Iapp = bw_g_L*RmIe_arrays[i_osc];

			if (Iapp < bw_Iapp_th)
			{
				perts_gsl[2*N_osc + bw_N_states*i_osc + 0] = bw_h_inf_V(bw_steady_V);
				perts_gsl[2*N_osc + bw_N_states*i_osc + 1] = bw_n_inf_V(bw_steady_V);
			}
			else
			{
				perts_gsl[2*N_osc + bw_N_states*i_osc + 0] = interpolate_BW_h(Iapp);
				perts_gsl[2*N_osc + bw_N_states*i_osc + 1] = interpolate_BW_n(Iapp);
			}

			break;
		}
		}
	}

	/* Initialize additional states for E-cells */
	for (i_osc = N_i; i_osc < N_osc; i_osc++)
	{
		double Iapp = nw_g_L*RmIe_arrays[i_osc];

		if (type_Pyr_cells == PYR_CELL_CA1)
		{// CA1 pyr. cells
			if ((Iapp < nw_CA1_Iapp_th) || (nw_CA1_Iapp_Maxth < Iapp))
			{
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 0] = nw_cal_x_inf(nw_steady_V, -75.0, -7.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 1] = nw_cal_x_inf(nw_steady_V, -54.0, 5.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 2] = nw_cal_x_inf(nw_steady_V, -65.0, -8.5);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 3] = nw_cal_x_inf(nw_steady_V, -15.0, 5.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 4] = nw_cal_x_inf(nw_steady_V, -60.0, -7.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 5] = nw_cal_x_inf(nw_steady_V, -5.8, 11.4);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 6] = nw_cal_x_inf(nw_steady_V, -68.0, -9.7);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 7] = nw_cal_x_inf(nw_steady_V, -30.0, 10.0);
			}
			else
			{
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 0] = interpolate_CA1_h_NaT(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 1] = interpolate_CA1_m_CaT(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 2] = interpolate_CA1_h_CaT(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 3] = interpolate_CA1_m_CaH(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 4] = interpolate_CA1_h_CaH(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 5] = interpolate_CA1_m_KDR(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 6] = interpolate_CA1_h_KDR(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 7] = interpolate_CA1_m_KM(Iapp);
			}
		}
		else
		{// CA3 pyr. cells
			if ((Iapp < nw_CA3_Iapp_th) || (nw_CA3_Iapp_Maxth < Iapp))
			{
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 0] = nw_cal_x_inf(nw_steady_V, -75.0, -7.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 1] = nw_cal_x_inf(nw_steady_V, -54.0, 5.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 2] = nw_cal_x_inf(nw_steady_V, -65.0, -8.5);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 3] = nw_cal_x_inf(nw_steady_V, -15.0, 5.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 4] = nw_cal_x_inf(nw_steady_V, -60.0, -7.0);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 5] = nw_cal_x_inf(nw_steady_V, -5.8, 11.4);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 6] = nw_cal_x_inf(nw_steady_V, -68.0, -9.7);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 7] = nw_cal_x_inf(nw_steady_V, -30.0, 10.0);
			}
			else
			{
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 0] = interpolate_CA3_h_NaT(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 1] = interpolate_CA3_m_CaT(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 2] = interpolate_CA3_h_CaT(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 3] = interpolate_CA3_m_CaH(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 4] = interpolate_CA3_h_CaH(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 5] = interpolate_CA3_m_KDR(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 6] = interpolate_CA3_h_KDR(Iapp);
				perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + 7] = interpolate_CA3_m_KM(Iapp);
			}
		}
	}

	ode_solver_params->LEs_i 		= 1;

	double intV_Tmax = 1.0;
	if (intV_Tmax < longest_obs_times)
	{
		intV_Tmax = longest_obs_times;
	}

	int pre_is_initV_mode = ode_solver_params->is_initV_mode;
	ode_solver_params->is_initV_mode = 1;

	double t = 0.0;
	while (t <= (intV_Tmax + 1.0))
	{
		/* Check if a neuron is ready to be initialized */
		for (i_osc = 0; i_osc < N_osc; i_osc++)
		{
			if (obs_times[i_osc] <= t)
			{
				v[i_osc] = perts_gsl[2*i_osc];

				if (i_osc < N_i)
				{
					switch(INTN_TYPE)
					{
					case WANG_BUZSAKI:
					{
						int j;
						for (j = 0; j < wb_N_states; j++)
						{
							wb_Icells[wb_N_states*i_osc + j] = perts_gsl[2*N_osc + wb_N_states*i_osc + j];
						}

						break;
					}
					case BORGER_WALKER:
					{
						int j;
						for (j = 0; j < bw_N_states; j++)
						{
							wb_Icells[bw_N_states*i_osc + j] = perts_gsl[2*N_osc + bw_N_states*i_osc + j];
						}

						break;
					}
					}

				}
				else
				{
					int j;
					for (j = 0; j < nw_N_states; j++)
					{
						mp_Ecells[nw_N_states*(i_osc - N_i) + j] = perts_gsl[2*N_osc + wb_N_states*N_i + nw_N_states*(i_osc - N_i) + j];
					}
				}

				obs_times[i_osc] = intV_Tmax + 10;	// To not visit the same condition.
			}
		}


		/* Evolution of the perturbed trajectory */
		double tmp_t = t;

		/* Initialize GSL ode solver */
		gsl_odeiv2_system 	tmp_sys_gsl	= {RHS_v_v_perts, NULL, total_N_states, ode_solver_params};
		gsl_odeiv2_driver 	*tmp_d_gsl 	= gsl_odeiv2_driver_alloc_y_new(&tmp_sys_gsl, GSL_TYPE, dt, EPSABS, EPSREL);

		/* Evolve one step further */
		int s_gsl = gsl_odeiv2_driver_apply(tmp_d_gsl, &tmp_t, tmp_t + dt, perts_gsl); gsl_odeiv2_driver_free(tmp_d_gsl);
		if (s_gsl != GSL_SUCCESS)
		{
			printf ("error: driver returned %d and exit(0)\n", s_gsl);
			exit(0);
		}

		t = t + dt;
	}

	/* Do they in the spike-generating state? */
	for (i_osc = 0; i_osc < N_osc; i_osc++)
	{
		if (i_osc < N_i)
		{// Interneurons
			switch(INTN_TYPE)
			{
			case WANG_BUZSAKI:
				if (m_wb_HH_th <= wb_get_m_inf(v[i_osc]))
				{
					state_wb_Icells[i_osc] = CELL_GEN_SPK_STATE;
				}
				else
				{
					state_wb_Icells[i_osc] = CELL_NOGEN_SPK_STATE;
				}

				break;
			case BORGER_WALKER:
				if (m_wb_HH_th <= bw_m_inf_V(v[i_osc]))
				{
					state_wb_Icells[i_osc] = CELL_GEN_SPK_STATE;
				}
				else
				{
					state_wb_Icells[i_osc] = CELL_NOGEN_SPK_STATE;
				}

				break;
			}
		}
		else
		{// Pyr. cells
			if  (m_mp_HH_th <= nw_cal_x_inf(v[i_osc], -37, 5))
			{
				state_mp_Ecells[i_osc - N_i] = CELL_GEN_SPK_STATE;
			}
			else
			{
				state_mp_Ecells[i_osc - N_i] = CELL_NOGEN_SPK_STATE;
			}
		}
	}


	ode_solver_params->is_initV_mode = pre_is_initV_mode;
}

post_syn_osc_struct	*getPre_syn_MatPtrEle(post_syn_osc_struct **pre_syn_Mat_ptr, int x, int y, int x_l, int y_l)
{
	return *(pre_syn_Mat_ptr + y+ x*y_l);
}

static double wb_get_I_gap(double V_post, const double V_dV[], post_syn_osc_struct *pre_syn_list)
{
	double I_gap = 0.0;

	while (pre_syn_list != NULL)
	{
		I_gap += wb_g_gap_junction*(V_post - V_dV[2*(pre_syn_list->id)]);
		pre_syn_list = pre_syn_list->next;
	}

	return I_gap;
}

static double nw_rhs_gating_variable(double V, double V_x, double k_x, double x, double tau_x)
{
    return (nw_cal_x_inf(V, V_x, k_x) - x)/tau_x;
}

static double nw_cal_x_inf(double V, double V_x, double k_x)
{
    return 1/(1 + exp(-(V - V_x)/k_x));
}

static double wb_get_n_inf(double v)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	return 1/(1 + exp(-0.045*(v + 10)));
#else
	double alpha_m = -0.01*(v + 34)/(exp(-0.1*(v + 34)) - 1);
	double beta_m = 0.125*exp(-(v + 44)/80.0);

	return alpha_m/(alpha_m + beta_m);
#endif
}

static double wb_get_h_inf(double v)
{
#ifdef WHITE_INTN_SUB_TYPE_OF_WANG_BUZSAKI
	return 1/(1 + exp(0.13*(v + 38)));
#else
	double alpha_m = 0.07*(exp(-(v + 58)/20));
	double beta_m = 1/(exp(-0.1*(v + 28)) + 1);

	return alpha_m/(alpha_m + beta_m);
#endif
}

static double bw_f_h(double V, double h)
{
	double h_inf = bw_h_inf_V(V);

	double alpha_h = bw_alpha_h_V(V);
	double beta_h = bw_beta_h_V(V);

	double tau_h = 1/(alpha_h + beta_h);

	return (h_inf - h)/tau_h;
}

static double bw_f_n(double V, double n)
{
	double n_inf = bw_n_inf_V(V);

	double alpha_n = bw_alpha_n_V(V);
	double beta_n = bw_beta_n_V(V);

	double tau_n = 1/(alpha_n + beta_n);

	return (n_inf - n)/tau_n;
}

static double bw_get_I_L(double g_L, double v)
{
	return g_L*(v - bw_E_L);
}

static double bw_get_I_K(double g_K, double v, double n)
{
	return g_K*n*n*(v - bw_E_K);
}

static double bw_get_I_Na(double g_Na, double v, double h)
{
	double m_inf = bw_m_inf_V(v);

	return g_Na*m_inf*m_inf*m_inf*h*(v - bw_E_Na);
}


static double bw_m_inf_V(double V)
{
	double alpha_m = bw_alpha_m_V(V);
	double beta_m = bw_beta_m_V(V);

	return alpha_m/(alpha_m + beta_m);
}

static double bw_h_inf_V(double V)
{
	double alpha_h = bw_alpha_h_V(V);
	double beta_h = bw_beta_h_V(V);

	return alpha_h/(alpha_h + beta_h);
}

static double bw_n_inf_V(double V)
{
	double alpha_n = bw_alpha_n_V(V);
	double beta_n = bw_beta_n_V(V);

	return alpha_n/(alpha_n + beta_n);
}

static double bw_alpha_m_V(double V)
{
	return 40.0*(75.5 - V)/(exp((75.5 - V)/13.5) - 1.0);
}

static double bw_beta_m_V(double V)
{
	return 1.2262/exp(V/42.248);
}

static double bw_alpha_h_V(double V)
{
	return 0.0035/exp(V/24.186);
}

static double bw_beta_h_V(double V)
{
	if (cmp(V, -51.25, 1e-3) == 0)
	{
		return 0.017*5.2;
	}else
	{
		return -0.017*(V + 51.25)/(exp(-(V + 51.25)/5.2) - 1.0);
	}
}

static double bw_alpha_n_V(double V)
{
	return (95.0 - V)/(exp((95.0 - V)/11.8) - 1.0);
}

static double bw_beta_n_V(double V)
{
	return 0.025/exp(V/22.222);
}

static void adapt_I0_Ecells(double t, int N_i, int N_e, double RmIe_arrays[], E_neuron_props_struct *E_neuron_p)
{
	if (t > (I0_E_t_begin + I0_E_t_int*(I0_E_running_I - 1)))
	{
		int i;
		for (i = N_i; i < (N_i + N_e); i++)
		{
			RmIe_arrays[i] = E_neuron_p->R_m*(I0_E_begin + (I0_E_at_which_I - 1)*(I0_E_end - I0_E_begin)/(I0_E_N - 1));
		}

		if (I0_E_mode == 1)
		{
			// Increasing current mode
			I0_E_at_which_I = I0_E_at_which_I + 1;
		}
		else
		{
			// Decreasing current mode
			I0_E_at_which_I = I0_E_at_which_I - 1;
		}

		I0_E_running_I = I0_E_running_I + 1;

		if (I0_E_at_which_I == (I0_E_N + 1))
		{
			I0_E_mode = I0_E_mode*(-1);
			I0_E_at_which_I = I0_E_N - 1;
		}

		if (I0_E_at_which_I == 0)
		{
			I0_E_mode = I0_E_mode*(-1);
			I0_E_at_which_I = 2;
		}
	}
}

static void adapt_I0_Icells(double t, int N_i, int N_e, double RmIe_arrays[], I_neuron_props_struct	*I_neuron_p)
{
	if (t > (I0_I_t_begin + I0_I_t_int*(I0_I_running_I - 1)))
	{
		double tmp_I0 = (I0_I_begin + (I0_I_at_which_I - 1)*(I0_I_end - I0_I_begin)/(I0_I_N - 1));
		printf("%d, %1.6lf\n", I0_I_at_which_I, tmp_I0);

		int i;
		for (i = 0; i < N_i; i++)
		{
			RmIe_arrays[i] = I_neuron_p->R_m*tmp_I0;
		}

		if (I0_I_mode == 1)
		{
			// Increasing current mode
			I0_I_at_which_I = I0_I_at_which_I + 1;
		}
		else
		{
			// Decreasing current mode
			I0_I_at_which_I = I0_I_at_which_I - 1;
		}

		I0_I_running_I = I0_I_running_I + 1;

		if (I0_I_at_which_I == (I0_I_N + 1))
		{
			I0_I_mode = I0_I_mode*(-1);
			I0_I_at_which_I = I0_I_N - 1;
		}

		if (I0_I_at_which_I == 0)
		{
			I0_I_mode = I0_I_mode*(-1);
			I0_I_at_which_I = 2;
		}
	}
}
