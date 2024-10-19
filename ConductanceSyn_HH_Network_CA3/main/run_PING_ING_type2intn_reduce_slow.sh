#!/usr/bin/perl
use Net::Domain qw(hostname hostfqdn hostdomain);
use Math::Complex;
my $host = hostname();

$POISSON_EXTERNAL_INPUTS 			=0;
$PERIODIC_EXTERNAL_INPUTS 			=1;
$PERIODIC_POISSON_EXTERNAL_INPUTS 	=2;

$PYR_CELL_CA1 = 0;
$PYR_CELL_CA3 = 1;

$IS_SELF_II_ALLOW 		=1;
$IS_SELF_EE_ALLOW 		=1;
$IS_SELF_II_NOTALLOW 	=0;
$IS_SELF_EE_NOTALLOW 	=0;

$REG_CA1 = 0;
$REG_CA3 = 1;

$nw_A_e = (21590*1e-8); # [cm^2]
$wb_A_i = (18069*1e-8); # [cm^2]

system("echo \"[`date`] Start script\"");

####################### Default parameters #######################
# Store the results #
$FN="Profile";

# Deterministically random, Not used anymore!!! #
$seed		=0;
$isUsingSeed=1;

# Network properties, 0,...,N_i-1=Inh.; N_i,...,N_i+N_e-1=Exc.; N_i+N_e=Ext. #
$N_e	=4000;  # The number of pyramidal neurons
$N_i	=1000;	# The number of interneurons
$N_ext	=1;  	# The number of external oscillators that generate the external inputs

# Typical pyramidal properties #
$V_e_rest    =-70;       # Resting membrane potential [mV]
$V_e_th      =-52;       # The spike threshold [mV]
$V_e_reset   =-59;       # The reset potential [mV]
$V_e_refrac  =2;         # Absolute refractory period [ms]
$tau_m_e     =50;        # !! no more for HH neurons,Membrane time constants [ms]. Need change to an effective time constants.
$C_m_e       =215.9;       # Capacitance [pF]

# Typical interneuron properties #
$V_i_rest    =-70;       # Resting membrane potential [mV]
$V_i_th      =-52;       # The spike threshold [mV]
$V_i_reset   =-59;       # The reset potential [mV]
$V_i_refrac  =1;         # Absolute refractory period [ms]
$tau_m_i     =2;        # !! no more for HH neurons, Membrane time constants [ms]
$C_m_i       =180.69;       # Capacitance [pF]

# Typical synaptic properties #
$g_syn_AMPA_on_E=8.598;		# [nS]
$g_syn_AMPA_on_I=6.150;    	# [nS]
$g_syn_GABA_on_E=24.06;    	# [nS]
$g_syn_GABA_on_I=11.18;    	# [nS]

$v_syn_AMPA=0;           	# [mV]
$v_syn_GABA=-75;         	# [mV], GABA_A

$tau_l_AMPA_on_E=2.50;         # CA3, Debanne, [ms],Synaptic latency
$tau_l_AMPA_on_I=1.30;         # CA1, [ms],Synaptic latency
$tau_l_GABA_on_E=0.95;         # CA1, [ms],Synaptic latency
$tau_l_GABA_on_I=0.60;         # CA1, [ms],Synaptic latency

$tau_r_AMPA_on_E=0.50;       # [ms]
$tau_r_AMPA_on_I=0.45;       # [ms]
$tau_r_GABA_on_E=0.25;       # [ms]
$tau_r_GABA_on_I=0.30;       # [ms]

$tau_d_AMPA_on_E=2.50;         # [ms]
$tau_d_AMPA_on_I=1.00;         # [ms]
$tau_d_GABA_on_E=4.00;         # [ms]
$tau_d_GABA_on_I=2.00;         # [ms]

# Typical external input properties #
$lambda_AMPA_on_E=0;            	# firing rate. Lambda spikes per ms.
$lambda_AMPA_on_I=0;            	# firing rate. Lambda spikes per ms.
$lambda_GABA_on_E=0;            	# firing rate. Lambda spikes per ms.
$lambda_GABA_on_I=0;            	# firing rate. Lambda spikes per ms.

$RmIe_E=0.0;				# [mV]
$RmIe_I=0.0;				# [mV]

# Simulation parameters #
$dt=0.05;     				# Step size. [ms]
$tEnd=1000;					# [ms]

$is_turn_on_pert=0;
$pert_size=0.0;
$show_each_second=1.0;	# [ms]
$is_rec_spkp=1;
$is_rec_traj=0;
$is_rec_LEs=1;
$is_rec_LEs_End=1;
$is_rec_LEs_End_Pert=1;
$is_rec_LEs_Pert=1;

$is_hetero_Ie_E = 1;
$is_hetero_Ie_I = 1;

$sigma_Ie_E = 0.0;
$sigma_Ie_I = 0.0;

$rnd_init_perts_seed = 1;
$is_rand_V_init_X = 1;
$rnd_init_V = 1;
$rnd_init_RmIe_seed = 1;
$rnd_determineConnections = 1;
$rnd_duringExecution = 1;

$sigma_white_noise_E = 0.0;	# [mV]
$sigma_white_noise_I = 0.0;	# [mV]
$rnd_white_noise_E = 1;
$rnd_white_noise_I = 1;

$is_use_GSL_2_solve_ODE = 0;		# dt = 0.05 when is_use_GSL_2_solve_ODE = 1 and dt = 0.01 when is_use_GSL_2_solve_ODE = 0

$type_of_external_inputs = $POISSON_EXTERNAL_INPUTS;

if ($host eq "patron")
{
$mjob_file = "./mjob_RAP_patron_type2int.sh";
$store_dir = "/local/";      
}
else
{
$mjob_file = "./mjob_RAP_type2int.sh";
$store_dir = "/vol/neuroinf/viriyopa/scratch/";	
}

$isGJExactNumberOfConns = 0;

##################################################################################################################################################################################
####################### Parameters #######################
$which_reg = $REG_CA1;

$N_pyr = 4000;
$N_basket = 1000;

$FN="ING";
$N_e=1;
$N_i=1;

####################### Parameters #######################
if ($which_reg == $REG_CA1)
{
# Realistic parameters for CA1
$N_pyr 		=300000;
$N_basket 	=5000;
$pyr_cell 	= $PYR_CELL_CA1; 	
}
else
{
# Realistic parameters for CA3
$N_pyr 		=200000;
$N_basket 	=1000;
$pyr_cell 	= $PYR_CELL_CA3;	
}

if ($which_reg == $REG_CA1)
{
# Probability for each pair groups of neurons for CA1 #
$p_EE=0.0067;   	# E to E
$p_EI=0.3000;   	# E to I
$p_IE=0.6700;   	# I to E
$p_II=0.0600;   	# I to I
}
else
{
# Probability for each pair groups of neurons for CA3 #
$p_EE=0.0200;   	# E to E
$p_EI=0.1000;   	# E to I
$p_IE=0.3000;   	# I to E
$p_II=0.0080;   	# I to I	
}		
##################################################################################################################################################################################
################################################## Dynamic of interneurons ##################################################
######## wb_g_L, [wb_g_L] = mS/cm^2
$begin_Val41		=0.1;$end_Val41		=$begin_Val41;$len_Val41	=1;
#$begin_Val41		=0.05;$end_Val41	=0.35;$len_Val41	=11;

################################################## External constant inputs ##################################################
######## Iapp_E, [Iapp_E] = pA
$is_Iapp_E_linear = 1;
#$begin_Val13		=0.0;$end_Val13		=$begin_Val13;$len_Val13	=1;
#$begin_Val13		=1258.9;$end_Val13=$begin_Val13;$len_Val13=1;					
#$begin_Val13		=3548.1;$end_Val13=$begin_Val13;$len_Val13=1;					
#$begin_Val13		=750.0;$end_Val13=$begin_Val13;$len_Val13=1;

#$begin_Val13		=60000.0;$end_Val13=$begin_Val13;$len_Val13=1;
$begin_Val13		=19000.0;$end_Val13=$begin_Val13;$len_Val13=1;

#$begin_Val13		=0.0;$end_Val13		=300000.0;$len_Val13			=11;
#$begin_Val13		=2400.0;$end_Val13		=19000.0;$len_Val13			=41;
#$begin_Val13		=300.0;$end_Val13	=300000.0;$len_Val13			=11;$is_Iapp_E_linear = 0;		
#$begin_Val13		=5650.0;$end_Val13	=18928.0;$len_Val13			=41;$is_Iapp_E_linear = 0;		


######## sigma_Ie_E, in term of percent of the mean current
$begin_Val16		=0.0;$end_Val16		=$begin_Val16;$len_Val16	=1;$is_hetero_Ie_E = 0;
#$begin_Val16		=100.0;$end_Val16		=$begin_Val16;$len_Val16	=1;$is_hetero_Ie_E = 1;
#$begin_Val16		=0.0;$end_Val16		=100.0;$len_Val16	=41;$is_hetero_Ie_E = 1;

######## Iapp_I, [Iapp_I] = pA
$is_Iapp_I_linear = 1;
$begin_Val1			=0.0;$end_Val1		=$begin_Val1;$len_Val1		=1;
#$begin_Val1			=90.0;$end_Val1=$begin_Val1;$len_Val1=1;					# Threshold current to excite a neuron is around 30.0 pA.
#$begin_Val1			=1280.0;$end_Val1=$begin_Val1;$len_Val1=1;					# Threshold current to excite a neuron is around 30.0 pA.
#$begin_Val1			=1385.0;$end_Val1=$begin_Val1;$len_Val1=1;					# Threshold current to excite a neuron is around 30.0 pA.
#$begin_Val1			=1520.0;$end_Val1=$begin_Val1;$len_Val1=1;					# Threshold current to excite a neuron is around 30.0 pA.
#$begin_Val1			=1175.0;$end_Val1		=2170.0;$len_Val1			=41;
#$begin_Val1			=1175.0;$end_Val1		=1770.0;$len_Val1			=41;
#$begin_Val1			=1806.0;$end_Val1		=2170.0;$len_Val1			=41;

#$begin_Val1			=0.0;$end_Val1		=2000.0;$len_Val1			=11;
#$begin_Val1			=1200.0;$end_Val1		=2000.0;$len_Val1			=11;
#$begin_Val1			=1200.0;$end_Val1		=1580.0;$len_Val1			=41;
#$begin_Val1			=1200.0;$end_Val1		=1580.0;$len_Val1			=41;
#$begin_Val1			=1200.0;$end_Val1		=1201.0;$len_Val1			=11;

#$begin_Val1			=0.0;$end_Val1		=20000.0;$len_Val1			=2001;
#$begin_Val1			=90.0;$end_Val1		=300.0;$len_Val1			=41;
#$begin_Val1			=1175.0;$end_Val1		=2170.0;$len_Val1		=11;$is_Iapp_I_linear = 0;

######## sigma_Ie_I, in term of percent of the mean current
$begin_Val4			=0.0;$end_Val4		=$begin_Val4;$len_Val4		=1;$is_hetero_Ie_I = 0;
#$begin_Val4			=200.0;$end_Val4		=$begin_Val4;$len_Val4		=1;$is_hetero_Ie_I = 1;
#$begin_Val4			=0.0;$end_Val4		=1000.0;$len_Val4	=11;$is_hetero_Ie_I = 1;

################################################## Poissson inputs ##################################################
######## N_ext
$begin_Val36		=1.0;$end_Val36		=$begin_Val36;$len_Val36	=1;

######## lambda_E_on_E, $val5*1000 Hz
$begin_Val28		=0.0;$end_Val28		=$begin_Val28;$len_Val28	=1;
#$begin_Val28		=3.0;$end_Val28=$begin_Val28;$len_Val28=1;					# Threshold
#$begin_Val28		=0.003;$end_Val28=$begin_Val28;$len_Val28=1;
#$begin_Val28		=0.3;$end_Val28=$begin_Val28;$len_Val28=1;
#$begin_Val28		=0.0;$end_Val28		=5.0;$len_Val28	=6;

######## lambda_E_on_I, $val5*1000 Hz
$begin_Val5			=0.0;$end_Val5		=$begin_Val5;$len_Val5		=1;
#$begin_Val5			=1.0;$end_Val5=$begin_Val5;$len_Val5=1;						# ING with ~50Hz.
#$begin_Val5			=0.05;$end_Val5=$begin_Val5;$len_Val5=1;
#$begin_Val5			=0.0;$end_Val5		=10.0;$len_Val5		=11;

######## lambda_I_on_E, $val5*1000 Hz
$begin_Val29		=0.0;$end_Val29		=$begin_Val29;$len_Val29	=1;
#$begin_Val29		=3.0;$end_Val29=$begin_Val29;$len_Val29=1;
#$begin_Val29		=0.0;$end_Val29		=5.0;$len_Val29=6;

######## lambda_I_on_I, $val5*1000 Hz
$begin_Val30		=0.0;$end_Val30		=$begin_Val30;$len_Val30	=1;
#$begin_Val30		=3.0;$end_Val30=$begin_Val30;$len_Val30=1;
#$begin_Val30		=0.0;$end_Val30		=10.0;$len_Val30	=11;

################################################## Current noise (White noise), dV/dt = f(V)/C_m + sigma/sqrt(tau_m)*w(t) ##################################################
$is_use_sigma_white_noise = 0;
######## sigma_white_noise_E
$begin_Val17		=0.0;$end_Val17		=$begin_Val17;$len_Val17	=1;
#$begin_Val17		=60.0;$end_Val17	=$begin_Val17;$len_Val17	=1;$is_use_sigma_white_noise=1;
#$begin_Val17		=0.0;$end_Val17	=150.0;$len_Val17			=11;$is_use_sigma_white_noise=1;

######## sigma_white_noise_I
$begin_Val12		=0.0;$end_Val12		=$begin_Val12;$len_Val12	=1;
#$begin_Val12		=0.5;$end_Val12	=$begin_Val12;$len_Val12	=1;$is_use_sigma_white_noise=1;
#$begin_Val12		=0.0;$end_Val12		=10.0;$len_Val12			=11;$is_use_sigma_white_noise=1;

################################################## tau_l of Synpases ##################################################
######## External Poisson Input also uses this values
######## tau_l_AMPA_on_E
$begin_Val38		=2.5;$end_Val38		=$begin_Val38;$len_Val38	=1;
#$begin_Val38		=3.0;$end_Val38		=$begin_Val38;$len_Val38		=1;
#$begin_Val38		=0.0;$end_Val38		=5.0;$len_Val38		=21;

######## tau_l_AMPA_on_I
$begin_Val22		=1.3;$end_Val22		=$begin_Val22;$len_Val22	=1;
#$begin_Val22		=0.8;$end_Val22		=1.8;$len_Val22		=5;

######## tau_l_GABA_on_E
$begin_Val25		=0.95;$end_Val25	=$begin_Val25;$len_Val25	=1;
#$begin_Val25		=0.5;$end_Val25	=$begin_Val25;$len_Val25		=1;
#$begin_Val25		=0.45;$end_Val25	=0.45;$len_Val25	=1;
#$begin_Val25		=0.45;$end_Val25	=1.45;$len_Val25	=5;

######## tau_l_GABA_on_I
$begin_Val7			=0.6;$end_Val7		=$begin_Val7;$len_Val7		=1;
#$begin_Val7			=1.1;$end_Val7		=$begin_Val7;$len_Val7		=1;
#$begin_Val7			=0.5;$end_Val7		=0.5;$len_Val7		=1;
#$begin_Val7		=0.5;$end_Val7		=1.5;$len_Val7		=5;

################################################## tau_r of Synpases ##################################################
######## External Poisson Input also uses this values
######## tau_r_AMPA_on_E
$begin_Val39		=0.5;$end_Val39		=$begin_Val39;$len_Val39	=1;
#$begin_Val39		=1.0;$end_Val39		=$begin_Val39;$len_Val39	=1;
#$begin_Val39		=0.1;$end_Val39		=1.1;$len_Val39		=21;

######## tau_r_AMPA_on_I
#$begin_Val23		=0.45;$end_Val23	=$begin_Val23;$len_Val23	=1;
$begin_Val23		=1.0;$end_Val23		=$begin_Val23;$len_Val23	=1;
#$begin_Val23		=0.1;$end_Val23		=0.1;$len_Val23	=1;
#$begin_Val23		=0.1;$end_Val23		=1.0;$len_Val23		=11;

######## tau_r_GABA_on_E
$begin_Val26		=0.25;$end_Val26	=$begin_Val26;$len_Val26	=1;
#$begin_Val26		=0.1;$end_Val26		=$begin_Val26;$len_Val26	=1;
#$begin_Val26		=0.1;$end_Val26	=0.1;$len_Val26	=1;
#$begin_Val26		=0.1;$end_Val26		=0.5;$len_Val26				=5;

######## tau_r_GABA_on_I
$begin_Val8			=0.3;$end_Val8		=$begin_Val8;$len_Val8		=1;
#$begin_Val8			=0.1;$end_Val8		=$begin_Val8;$len_Val8		=1;
#$begin_Val8			=0.1;$end_Val8		=0.1;$len_Val8		=1;
#$begin_Val8			=0.1;$end_Val8		=0.5;$len_Val8				=5;

################################################## tau_d of Synpases ##################################################
######## External Poisson Input also uses this values
######## tau_d_AMPA_on_E
$begin_Val40		=2.5;$end_Val40		=$begin_Val40;$len_Val40	=1;
#$begin_Val40		=3.0;$end_Val40		=$begin_Val40;$len_Val40	=1;
#$begin_Val40		=2.0;$end_Val40		=3.0;$len_Val40		=21;

######## tau_d_AMPA_on_I
$begin_Val24		=1.00;$end_Val24	=$begin_Val24;$len_Val24	=1;
#$begin_Val24		=0.5;$end_Val24		=$begin_Val24;$len_Val24	=1;
#$begin_Val24		=0.7;$end_Val24		=0.7;$len_Val24	=1;
#$begin_Val24		=0.7;$end_Val24		=1.7;$len_Val24		=5;

######## tau_d_GABA_on_E
$begin_Val27		=4.00;$end_Val27	=$begin_Val27;$len_Val27	=1;
#$begin_Val27		=3.5;$end_Val27		=$begin_Val27;$len_Val27	=1;
#$begin_Val27		=3.5;$end_Val27		=3.5;$len_Val27	=1;
#$begin_Val27		=1.5;$end_Val27		=4.0;$len_Val27				=5;

######## tau_d_GABA_on_I
$begin_Val9			=2.0;$end_Val9		=$begin_Val9;$len_Val9		=1;
#$begin_Val9			=2.5;$end_Val9		=$begin_Val9;$len_Val9		=1;
#$begin_Val9			=1.5;$end_Val9		=4.0;$len_Val9				=5;

################################################## Condutance of Chemical Synpases ##################################################
######## External Poisson Input also uses this values
######## g_syn_AMPA_on_E, gE->E
$begin_Val21		=8.598;$end_Val21	=$begin_Val21;$len_Val21	=1;	# Default CA1
#$begin_Val21		=9.408;$end_Val21=$begin_Val21;$len_Val21=1;	# Bal.
#$begin_Val21		=8.598*$p_EE*4000;$end_Val21	=$begin_Val21;$len_Val21	=1;	# Default CA1
#$begin_Val21		=9.598;$end_Val21=$begin_Val21;$len_Val21=1;	# old.		

######## g_syn_AMPA_on_I, gE->I
$begin_Val14		=6.150;$end_Val14	=$begin_Val14;$len_Val14	=1;	# Default gEI
#$begin_Val14		=4.430;$end_Val14=$begin_Val14;$len_Val14=1;	# Bal.
#$begin_Val14		=6.150*$p_EI*4000;$end_Val14	=$begin_Val14;$len_Val14	=1;	# Default gEI
#$begin_Val14		=0.0;$end_Val14=$begin_Val14;$len_Val14=1;
#$begin_Val14		=0.0;$end_Val14			=6.150;$len_Val14		=11;	

######## g_syn_GABA_on_E, gI->E
$begin_Val15		=24.06;$end_Val15	=$begin_Val15;$len_Val15	=1;	# Default gIE
#$begin_Val15		=23.739;$end_Val15=$begin_Val15;$len_Val15=1;	# Bal.
#$begin_Val15		=24.06*$p_IE*1000;$end_Val15	=$begin_Val15;$len_Val15	=1;	# Default gIE
#$begin_Val15		=23.06;$end_Val15=$begin_Val15;$len_Val15=1;	
#$begin_Val15		=0.0;$end_Val15		=50.0;$len_Val15		=6;	

######## g_syn_GABA_on_I, gI->I
$is_gII_linear = 1;
$begin_Val2			=11.18;$end_Val2	=$begin_Val2;$len_Val2		=1;	# Default gII
#$begin_Val2			=0.0;$end_Val2		=$begin_Val2;$len_Val2		=1;	
#$begin_Val2			=11.18*0.3*1000;$end_Val2		=$begin_Val2;$len_Val2		=1;	# Default gII
#$begin_Val2			=5.0;$end_Val2			=15.0;$len_Val2		=11;
#$begin_Val2			=0.0;$end_Val2		=20.0;$len_Val2				=11;
#$begin_Val2			=0.01;$end_Val2		=20.0;$len_Val2				=41;$is_gII_linear = 0;

################################################## Connection probability of Chemical Synpases ##################################################
######################## p_EE
if ($which_reg == $REG_CA1)
{
$begin_Val18		=0.0067;$end_Val18	=$begin_Val18;$len_Val18	=1; # default CA1
$tmp_p_EE = $begin_Val18;
}
else
{
$begin_Val18		=0.02;$end_Val18	=$begin_Val18;$len_Val18	=1; # default CA3
}

#$begin_Val18		=0.0;$end_Val18		=$begin_Val18;$len_Val18	=1;
$begin_Val18		=1.0;$end_Val18		=$begin_Val18;$len_Val18	=1;

######################## p_EI
if ($which_reg == $REG_CA1)
{
$begin_Val19		=0.3;$end_Val19		=$begin_Val19;$len_Val19	=1; # default CA1
$tmp_p_EI = $begin_Val19;
}
else
{
$begin_Val19		=0.1;$end_Val19		=$begin_Val19;$len_Val19	=1; # default CA3
}

#$begin_Val19		=0.0;$end_Val19		=$begin_Val19;$len_Val19	=1;
$begin_Val19		=1.0;$end_Val19		=$begin_Val19;$len_Val19	=1;
#$begin_Val19		=0.0;$end_Val19		=0.4;$len_Val19		=6;

######################## p_IE
if ($which_reg == $REG_CA1)
{
$begin_Val20		=0.67;$end_Val20	=$begin_Val20;$len_Val20	=1; # default CA1
$tmp_p_IE = $begin_Val20;
}
else
{
$begin_Val20		=0.3;$end_Val20		=$begin_Val20;$len_Val20	=1; # default CA3
}

##$begin_Val20		=0.0;$end_Val20		=$begin_Val20;$len_Val20	=1;
$begin_Val20		=1.0;$end_Val20		=$begin_Val20;$len_Val20	=1;
#$begin_Val20		=0.0;$end_Val20		=0.7;$len_Val20		=6;

######################## p_II
$is_pII_linear = 1;
if ($which_reg == $REG_CA1)
{
$begin_Val3			=0.3;$end_Val3		=$begin_Val3;$len_Val3		=1; # default CA1
$tmp_p_II = $begin_Val3;
}
else
{
$begin_Val3			=0.3;$end_Val3		=$begin_Val3;$len_Val3		=1; # default CA3. This value is actually 0.008.
}

#$begin_Val3			=0.0;$end_Val3		=$begin_Val3;$len_Val3		=1;
$begin_Val3			=1.0;$end_Val3		=$begin_Val3;$len_Val3		=1;
#$begin_Val3			=0.008;$end_Val3=$begin_Val3;$len_Val3=1;			# default CA3
#$begin_Val3			=0.06;$end_Val3=$begin_Val3;$len_Val3=1;			# default CA1
#$begin_Val3			=0.3;$end_Val3=$begin_Val3;$len_Val3=1;			
#$begin_Val3			=0.017;$end_Val3		=0.05;$len_Val3		=11;
#$begin_Val3			=0.001;$end_Val3		=0.4;$len_Val3		=11;$is_pII_linear = 0;
#$begin_Val3			=0.001;$end_Val3		=1.0;$len_Val3		=11;$is_pII_linear = 0;

################################################## Condutance of Synpases, Gap junction ##################################################
######## g_syn_II_gap_junction
$begin_Val35		=1.8069;$end_Val35	=$begin_Val35;$len_Val35	=1;	# CA1, Bartos 2002, PNAS
#$begin_Val35		=0.0;$end_Val35		=$begin_Val35;$len_Val35	=1;
#$begin_Val35		=0.0;$end_Val35		=4.0;$len_Val35				=41;


################################################## Connection probability, Gap Junction ##################################################
######## p_EE_gap_junction
$begin_Val32		=0.0;$end_Val32		=$begin_Val32;$len_Val32	=1;			# default, Bartos 2007

######## p_EI_gap_junction
$begin_Val33		=0.0;$end_Val33		=$begin_Val33;$len_Val33	=1;			# default, Bartos 2007

######## p_IE_gap_junction
$begin_Val34		=0.0;$end_Val34		=$begin_Val34;$len_Val34	=1;			# default, Bartos 2007

######## p_II_gap_junction
$is_pGJII_linear = 1;
$begin_Val31		=0.0;$end_Val31		=$begin_Val31;$len_Val31	=1;
#$begin_Val31		=0.004;$end_Val31	=$begin_Val31;$len_Val31	=1;
#$begin_Val31		=0.001;$end_Val31		=0.01;$len_Val31	=11;$is_pGJII_linear = 0;

################################################## Simulation ##################################################
################ Display 
$show_each_second=100.0;	# [ms]
#$show_each_second=10000.0;	# [ms]

######## Run multiple jobs?
#$blast_job = 1;$wait2see = 0;
$blast_job = 0;$wait2see = 1;
#$blast_job = 0;$wait2see = 0;

######## Record trajectories?
#$is_rec_traj=0;
$is_rec_traj=1;

######## Do we randomize the initial voltages?
$is_rand_V_init_X = 1;
#$is_rand_V_init_X = 0;

######## rnd_init_V
$begin_Val11		=0;$end_Val11		=$begin_Val11;$len_Val11	=1;
#$begin_Val11	=2;$end_Val11	=3;$len_Val11	=2;

######## $rnd_init_RmIe_seed
$rnd_init_RmIe_seed = 0;

######## ODE solver 
if ($is_use_sigma_white_noise == 0)
{
#$is_use_GSL_2_solve_ODE = 1;$begin_Val10=0.05;$end_Val10=$begin_Val10;$len_Val10=1;			# For GSL, dt = 0.05 is not so different from dt = 0.01.
$is_use_GSL_2_solve_ODE = 1;$begin_Val10=0.01;$end_Val10=$begin_Val10;$len_Val10=1;			# For GSL, dt = 0.05 is not so different from dt = 0.01.
}
else
{
$is_use_GSL_2_solve_ODE = 0;$begin_Val10=0.01;$end_Val10=$begin_Val10;$len_Val10=1;		# dt = 0.001 to get close to GSL results		
}

######## Type of External input
$type_of_external_inputs = $POISSON_EXTERNAL_INPUTS;
#$type_of_external_inputs = $PERIODIC_EXTERNAL_INPUTS;

######## tEnd
#$begin_Val37	=2000.0;$end_Val37=$begin_Val37;$len_Val37=1;
$begin_Val37	=500.0;$end_Val37=$begin_Val37;$len_Val37=1;
#$begin_Val37	=10000.0;$end_Val37=$begin_Val37;$len_Val37=1;
#$begin_Val37	=2000.0;$end_Val37	=20000.0;$len_Val37=11;

######## Self connection of a cell
$allow_EE_connections = $IS_SELF_EE_ALLOW;
$allow_II_connections = $IS_SELF_II_ALLOW;

#$allow_EE_connections = $IS_SELF_EE_NOTALLOW;
#$allow_II_connections = $IS_SELF_II_NOTALLOW;

######## Exact number of connections of the gap junctions
$isGJExactNumberOfConns = 0;
#$isGJExactNumberOfConns = 1;

########################################################################################################################

if ($len_Val1 > 1){$tmp = $end_Val1 - $begin_Val1;$tmp1 = $len_Val1 - 1;$dx_Val1 = $tmp/$tmp1;}else{$dx_Val1 = 0;}
if ($len_Val2 > 1){$tmp = $end_Val2 - $begin_Val2;$tmp1 = $len_Val2 - 1;$dx_Val2 = $tmp/$tmp1;}else{$dx_Val2 = 0;}
if ($len_Val3 > 1){$tmp = $end_Val3 - $begin_Val3;$tmp1 = $len_Val3 - 1;$dx_Val3 = $tmp/$tmp1;}else{$dx_Val3 = 0;}
if ($len_Val4 > 1){$tmp = $end_Val4 - $begin_Val4;$tmp1 = $len_Val4 - 1;$dx_Val4 = $tmp/$tmp1;}else{$dx_Val4 = 0;}
if ($len_Val5 > 1){$tmp = $end_Val5 - $begin_Val5;$tmp1 = $len_Val5 - 1;$dx_Val5 = $tmp/$tmp1;}else{$dx_Val5 = 0;}
if ($len_Val7 > 1){$tmp = $end_Val7 - $begin_Val7;$tmp1 = $len_Val7 - 1;$dx_Val7 = $tmp/$tmp1;}else{$dx_Val7 = 0;}
if ($len_Val8 > 1){$tmp = $end_Val8 - $begin_Val8;$tmp1 = $len_Val8 - 1;$dx_Val8 = $tmp/$tmp1;}else{$dx_Val8 = 0;}
if ($len_Val9 > 1){$tmp = $end_Val9 - $begin_Val9;$tmp1 = $len_Val9 - 1;$dx_Val9 = $tmp/$tmp1;}else{$dx_Val9 = 0;}
if ($len_Val10 > 1){$tmp = $end_Val10 - $begin_Val10;$tmp1 = $len_Val10 - 1;$dx_Val10 = $tmp/$tmp1;}else{$dx_Val10 = 0;}
if ($len_Val11 > 1){$tmp = $end_Val11 - $begin_Val11;$tmp1 = $len_Val11 - 1;$dx_Val11 = $tmp/$tmp1;}else{$dx_Val11 = 0;}
if ($len_Val12 > 1){$tmp = $end_Val12 - $begin_Val12;$tmp1 = $len_Val12 - 1;$dx_Val12 = $tmp/$tmp1;}else{$dx_Val12 = 0;}
if ($len_Val13 > 1){$tmp = $end_Val13 - $begin_Val13;$tmp1 = $len_Val13 - 1;$dx_Val13 = $tmp/$tmp1;}else{$dx_Val13 = 0;}
if ($len_Val14 > 1){$tmp = $end_Val14 - $begin_Val14;$tmp1 = $len_Val14 - 1;$dx_Val14 = $tmp/$tmp1;}else{$dx_Val14 = 0;}
if ($len_Val15 > 1){$tmp = $end_Val15 - $begin_Val15;$tmp1 = $len_Val15 - 1;$dx_Val15 = $tmp/$tmp1;}else{$dx_Val15 = 0;}
if ($len_Val16 > 1){$tmp = $end_Val16 - $begin_Val16;$tmp1 = $len_Val16 - 1;$dx_Val16 = $tmp/$tmp1;}else{$dx_Val16 = 0;}
if ($len_Val17 > 1){$tmp = $end_Val17 - $begin_Val17;$tmp1 = $len_Val17 - 1;$dx_Val17 = $tmp/$tmp1;}else{$dx_Val17 = 0;}
if ($len_Val18 > 1){$tmp = $end_Val18 - $begin_Val18;$tmp1 = $len_Val18 - 1;$dx_Val18 = $tmp/$tmp1;}else{$dx_Val18 = 0;}
if ($len_Val19 > 1){$tmp = $end_Val19 - $begin_Val19;$tmp1 = $len_Val19 - 1;$dx_Val19 = $tmp/$tmp1;}else{$dx_Val19 = 0;}
if ($len_Val20 > 1){$tmp = $end_Val20 - $begin_Val20;$tmp1 = $len_Val20 - 1;$dx_Val20 = $tmp/$tmp1;}else{$dx_Val20 = 0;}
if ($len_Val21 > 1){$tmp = $end_Val21 - $begin_Val21;$tmp1 = $len_Val21 - 1;$dx_Val21 = $tmp/$tmp1;}else{$dx_Val21 = 0;}
if ($len_Val22 > 1){$tmp = $end_Val22 - $begin_Val22;$tmp1 = $len_Val22 - 1;$dx_Val22 = $tmp/$tmp1;}else{$dx_Val22 = 0;}
if ($len_Val23 > 1){$tmp = $end_Val23 - $begin_Val23;$tmp1 = $len_Val23 - 1;$dx_Val23 = $tmp/$tmp1;}else{$dx_Val23 = 0;}
if ($len_Val24 > 1){$tmp = $end_Val24 - $begin_Val24;$tmp1 = $len_Val24 - 1;$dx_Val24 = $tmp/$tmp1;}else{$dx_Val24 = 0;}
if ($len_Val25 > 1){$tmp = $end_Val25 - $begin_Val25;$tmp1 = $len_Val25 - 1;$dx_Val25 = $tmp/$tmp1;}else{$dx_Val25 = 0;}
if ($len_Val26 > 1){$tmp = $end_Val26 - $begin_Val26;$tmp1 = $len_Val26 - 1;$dx_Val26 = $tmp/$tmp1;}else{$dx_Val26 = 0;}
if ($len_Val27 > 1){$tmp = $end_Val27 - $begin_Val27;$tmp1 = $len_Val27 - 1;$dx_Val27 = $tmp/$tmp1;}else{$dx_Val27 = 0;}
if ($len_Val28 > 1){$tmp = $end_Val28 - $begin_Val28;$tmp1 = $len_Val28 - 1;$dx_Val28 = $tmp/$tmp1;}else{$dx_Val28 = 0;}
if ($len_Val29 > 1){$tmp = $end_Val29 - $begin_Val29;$tmp1 = $len_Val29 - 1;$dx_Val29 = $tmp/$tmp1;}else{$dx_Val29 = 0;}
if ($len_Val30 > 1){$tmp = $end_Val30 - $begin_Val30;$tmp1 = $len_Val30 - 1;$dx_Val30 = $tmp/$tmp1;}else{$dx_Val30 = 0;}
if ($len_Val31 > 1){$tmp = $end_Val31 - $begin_Val31;$tmp1 = $len_Val31 - 1;$dx_Val31 = $tmp/$tmp1;}else{$dx_Val31 = 0;}
if ($len_Val32 > 1){$tmp = $end_Val32 - $begin_Val32;$tmp1 = $len_Val32 - 1;$dx_Val32 = $tmp/$tmp1;}else{$dx_Val32 = 0;}
if ($len_Val33 > 1){$tmp = $end_Val33 - $begin_Val33;$tmp1 = $len_Val33 - 1;$dx_Val33 = $tmp/$tmp1;}else{$dx_Val33 = 0;}
if ($len_Val34 > 1){$tmp = $end_Val34 - $begin_Val34;$tmp1 = $len_Val34 - 1;$dx_Val34 = $tmp/$tmp1;}else{$dx_Val34 = 0;}
if ($len_Val35 > 1){$tmp = $end_Val35 - $begin_Val35;$tmp1 = $len_Val35 - 1;$dx_Val35 = $tmp/$tmp1;}else{$dx_Val35 = 0;}
if ($len_Val36 > 1){$tmp = $end_Val36 - $begin_Val36;$tmp1 = $len_Val36 - 1;$dx_Val36 = $tmp/$tmp1;}else{$dx_Val36 = 0;}
if ($len_Val37 > 1){$tmp = $end_Val37 - $begin_Val37;$tmp1 = $len_Val37 - 1;$dx_Val37 = $tmp/$tmp1;}else{$dx_Val37 = 0;}
if ($len_Val38 > 1){$tmp = $end_Val38 - $begin_Val38;$tmp1 = $len_Val38 - 1;$dx_Val38 = $tmp/$tmp1;}else{$dx_Val38 = 0;}
if ($len_Val39 > 1){$tmp = $end_Val39 - $begin_Val39;$tmp1 = $len_Val39 - 1;$dx_Val39 = $tmp/$tmp1;}else{$dx_Val39 = 0;}
if ($len_Val40 > 1){$tmp = $end_Val40 - $begin_Val40;$tmp1 = $len_Val40 - 1;$dx_Val40 = $tmp/$tmp1;}else{$dx_Val40 = 0;}
if ($len_Val41 > 1){$tmp = $end_Val41 - $begin_Val41;$tmp1 = $len_Val41 - 1;$dx_Val41 = $tmp/$tmp1;}else{$dx_Val41 = 0;}

for ($i = 0; $i < $len_Val1; $i++){
for ($j = 0; $j < $len_Val2; $j++){		
for ($k = 0; $k < $len_Val3; $k++){	
for ($l = 0; $l < $len_Val4; $l++){
for ($m = 0; $m < $len_Val5; $m++){										
for ($o = 0; $o < $len_Val7; $o++){			
for ($p = 0; $p < $len_Val8; $p++){			
for ($q = 0; $q < $len_Val9; $q++){
for ($r = 0; $r < $len_Val10; $r++){
for ($s = 0; $s < $len_Val11; $s++){
for ($t = 0; $t < $len_Val12; $t++){
for ($u = 0; $u < $len_Val13; $u++){
for ($v = 0; $v < $len_Val14; $v++){
for ($w = 0; $w < $len_Val15; $w++){
for ($x = 0; $x < $len_Val16; $x++){																																																																																																																																					
for ($y = 0; $y < $len_Val17; $y++){
for ($z = 0; $z < $len_Val18; $z++){
for ($a = 0; $a < $len_Val19; $a++){
for ($b = 0; $b < $len_Val20; $b++){
for ($c = 0; $c < $len_Val21; $c++){
for ($d = 0; $d < $len_Val22; $d++){
for ($e = 0; $e < $len_Val23; $e++){
for ($f = 0; $f < $len_Val24; $f++){			
for ($g = 0; $g < $len_Val25; $g++){
for ($h = 0; $h < $len_Val26; $h++){
for ($h1 = 0; $h1 < $len_Val27; $h1++){
for ($h2 = 0; $h2 < $len_Val28; $h2++){
for ($h3 = 0; $h3 < $len_Val29; $h3++){
for ($h4 = 0; $h4 < $len_Val30; $h4++){
for ($h5 = 0; $h5 < $len_Val31; $h5++){
for ($h6 = 0; $h6 < $len_Val32; $h6++){
for ($h7 = 0; $h7 < $len_Val33; $h7++){
for ($h8 = 0; $h8 < $len_Val34; $h8++){
for ($h9 = 0; $h9 < $len_Val35; $h9++){
for ($h10 = 0; $h10 < $len_Val36; $h10++){
for ($h11 = 0; $h11 < $len_Val37; $h11++){
for ($h12 = 0; $h12 < $len_Val38; $h12++){
for ($h13 = 0; $h13 < $len_Val39; $h13++){
for ($h14 = 0; $h14 < $len_Val40; $h14++){
for ($h15 = 0; $h15 < $len_Val41; $h15++){	
																														
if ($is_Iapp_I_linear == 0)
{
$val1 = $begin_Val1*(($end_Val1/$begin_Val1)**($i/($len_Val1 - 1.0)));
}
else
{
$val1 = $begin_Val1 + $i*$dx_Val1;	
}

if ($is_gII_linear == 0)
{
$val2 = $begin_Val2*(($end_Val2/$begin_Val2)**($j/($len_Val2 - 1.0)));
}
else
{
$val2 = $begin_Val2 + $j*$dx_Val2;	
}


if ($is_pII_linear == 0)
{
$val3 = $begin_Val3*(($end_Val3/$begin_Val3)**($k/($len_Val3 - 1.0)));
}
else
{
$val3 = $begin_Val3 + $k*$dx_Val3;
}

$val4 = $begin_Val4 + $l*$dx_Val4;
$val5 = $begin_Val5 + $m*$dx_Val5;
$val7 = $begin_Val7 + $o*$dx_Val7;
$val8 = $begin_Val8 + $p*$dx_Val8;
$val9 = $begin_Val9 + $q*$dx_Val9;
$val10 = $begin_Val10 + $r*$dx_Val10;
$val11 = $begin_Val11 + $s*$dx_Val11;
$val12 = $begin_Val12 + $t*$dx_Val12;

if ($is_Iapp_E_linear == 0)
{
$val13 = $begin_Val13*(($end_Val13/$begin_Val13)**($u/($len_Val13 - 1.0)));
}
else
{
$val13 = $begin_Val13 + $u*$dx_Val13;	
}

$val14 = $begin_Val14 + $v*$dx_Val14;
$val15 = $begin_Val15 + $w*$dx_Val15;
$val16 = $begin_Val16 + $x*$dx_Val16;
$val17 = $begin_Val17 + $y*$dx_Val17;
$val18 = $begin_Val18 + $z*$dx_Val18;
$val19 = $begin_Val19 + $a*$dx_Val19;
$val20 = $begin_Val20 + $b*$dx_Val20;
$val21 = $begin_Val21 + $c*$dx_Val21;
$val22 = $begin_Val22 + $d*$dx_Val22;
$val23 = $begin_Val23 + $e*$dx_Val23;
$val24 = $begin_Val24 + $f*$dx_Val24;
$val25 = $begin_Val25 + $g*$dx_Val25;
$val26 = $begin_Val26 + $h*$dx_Val26;
$val27 = $begin_Val27 + $h1*$dx_Val27;
$val28 = $begin_Val28 + $h2*$dx_Val28;
$val29 = $begin_Val29 + $h3*$dx_Val29;
$val30 = $begin_Val30 + $h4*$dx_Val30;

if ($is_pGJII_linear == 0)
{
$val31 = $begin_Val31*(($end_Val31/$begin_Val31)**($h5/($len_Val31 - 1.0)));
}
else
{
$val31 = $begin_Val31 + $h5*$dx_Val31;
}

$val32 = $begin_Val32 + $h6*$dx_Val32;
$val33 = $begin_Val33 + $h7*$dx_Val33;
$val34 = $begin_Val34 + $h8*$dx_Val34;
$val35 = $begin_Val35 + $h9*$dx_Val35;
$val36 = $begin_Val36 + $h10*$dx_Val36;
$val37 = $begin_Val37 + $h11*$dx_Val37;
$val38 = $begin_Val38 + $h12*$dx_Val38;
$val39 = $begin_Val39 + $h13*$dx_Val39;
$val40 = $begin_Val40 + $h14*$dx_Val40;
$val41 = $begin_Val41 + $h15*$dx_Val41;
	
$wb_g_L = $val41;

$Iapp_E = $val13;
$Iapp_I = $val1; 

$RmE = $tau_m_e/$C_m_e; 
$RmI = $tau_m_i/$C_m_i;
		
$RmIe_E=$RmE*$Iapp_E;	
$RmIe_I=$RmI*$Iapp_I;
		
$dt=$val10;
$tEnd=$val37;
	
$is_rec_spkp=1;
$is_rec_LEs=0;
$is_rec_LEs_End=0;
$is_rec_LEs_End_Pert=0;
$is_rec_LEs_Pert=0;
		
# Conductance #
$g_syn_AMPA_on_E=$val21;
$g_syn_AMPA_on_I=$val14;
$g_syn_GABA_on_E=$val15;
$g_syn_GABA_on_I=$val2;

#$g_syn_GABA_on_I=$val2*1000.0/$N_i;

#if ($N_e > 0){$g_syn_AMPA_on_E=$val21*sqrt($N_pyr/$N_e);}else{$g_syn_AMPA_on_E=$val21;}
#if ($N_e > 0){$g_syn_AMPA_on_I=$val14*sqrt($N_pyr/$N_e);}else{$g_syn_AMPA_on_I=$val14;}
#if ($N_i > 0){$g_syn_GABA_on_E=$val15*$N_basket/$N_i;}else{$g_syn_GABA_on_E=$val15;}
#if ($N_i > 0){$g_syn_GABA_on_I=$val2*$N_basket/$N_i;}else{$g_syn_GABA_on_I=$val2;}

#$E_firing_rate = 40.0; #Hz
$E_firing_rate = 3; #Hz
$Gamma_freq = 40.0; #Hz

#$g_syn_AMPA_on_E=$g_syn_AMPA_on_E*$N_pyr*$p_EE;
#$g_syn_AMPA_on_I=$g_syn_AMPA_on_I*$N_pyr*$p_EI;

$g_syn_AMPA_on_E=($g_syn_AMPA_on_E*$N_pyr*$tmp_p_EE*3.0)/$Gamma_freq;
$g_syn_AMPA_on_I=$g_syn_AMPA_on_I*$N_pyr*$tmp_p_EI*$E_firing_rate/$Gamma_freq;

$g_syn_GABA_on_E=$g_syn_GABA_on_E*$N_basket*$tmp_p_IE;
$g_syn_GABA_on_I=$g_syn_GABA_on_I*$N_basket*$tmp_p_II;

$g_syn_II_gap_junction=$val35;

# Connection probability #
$p_EE=$val18;
$p_EI=$val19;
$p_IE=$val20;
$p_II=$val3;

$p_EE_gap_junction=$val32;
$p_EI_gap_junction=$val33;
$p_IE_gap_junction=$val34;
$p_II_gap_junction=$val31;

# Variation of external currents #
$sigma_Ie_E=$val16;
$sigma_Ie_I=$val4;

# White noise
$sigma_white_noise_E = $val17;
$sigma_white_noise_I = $val12;

# Poisson noise #
$N_ext = $val36;
$lambda_AMPA_on_E=$val28;
$lambda_AMPA_on_I=$val5;
$lambda_GABA_on_E=$val29;
$lambda_GABA_on_I=$val30;

# Kinetic of connection
$tau_l_AMPA_on_E=$val38;
$tau_r_AMPA_on_E=$val39;
$tau_d_AMPA_on_E=$val40;

$tau_l_AMPA_on_I=$val22;
$tau_r_AMPA_on_I=$val23;
$tau_d_AMPA_on_I=$val24;

$tau_d_AMPA_on_I=$tau_r_AMPA_on_I/0.45;

$tau_l_GABA_on_E=$val25;
$tau_r_GABA_on_E=$val26;
$tau_d_GABA_on_E=$val27;

$tau_l_GABA_on_I=$val7;
$tau_r_GABA_on_I=$val8;
$tau_d_GABA_on_I=$val9;

# Random numbers
$rnd_init_V = $val11;

if ($N_e == 0)
{
$REG_TAG = "BW_";
}
else
{
if ($which_reg == $REG_CA3)
{
$REG_TAG = "NWCA3_";
}
else
{
$REG_TAG = "NWCA1_";
}
}

#$FN=$REG_TAG."EEEIIEII_sigmaIeEx${sigma_Ie_E}_sigmaWNEx${sigma_white_noise_E}_IappE${Iapp_E}_rndV${rnd_init_V}_rndRmIe${rnd_init_RmIe_seed}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_gIIGJx${g_syn_II_gap_junction}_sigmaWNIx${sigma_white_noise_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_IappIx${Iapp_I}_sigmaWNE${sigma_white_noise_E}_pEI${p_EI}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_IappIx${Iapp_I}_sigmaWNE${sigma_white_noise_E}_sigmaWNI${sigma_white_noise_I}_pEI${p_EI}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaIeEx${sigma_Ie_E}_IappI${Iapp_I}_sigmaWNE${sigma_white_noise_E}_pEI${p_EI}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNIx${sigma_white_noise_I}_IappI${Iapp_I}_sigmaWNE${sigma_white_noise_E}_pEI${p_EI}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_lEEx${tau_l_AMPA_on_E}_rEEx${tau_r_AMPA_on_E}_dEEx${tau_d_AMPA_on_E}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}";

$FN=$REG_TAG."test_biexp_slow";
#$FN=$REG_TAG."test_biexp_fast2";
#$FN=$REG_TAG."test_normal";
#$FN=$REG_TAG."IappEx${Iapp_E}";
#$FN=$REG_TAG."IappIx${Iapp_I}";
#$FN="test_IappEx${Iapp_E}_dt${dt}_gEEx${g_syn_AMPA_on_E}_gEIx${g_syn_AMPA_on_I}_gIEx${g_syn_GABA_on_E}_gIIx${g_syn_GABA_on_I}";
#$FN="test_dt${dt}_is_use_GSL_2_solve_ODE${is_use_GSL_2_solve_ODE}";
#$FN="test_pIIGJ${p_II_gap_junction}_pII${p_II}_Iapp_I${Iapp_I}_gII${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}";
#$FN="test_gIIx${g_syn_GABA_on_I}_IappIx${Iapp_I}";
#$FN="test_lambda_AMPA_on_E${lambda_AMPA_on_E}_dt${dt}";
#$FN="test_lambda_AMPA_on_E${lambda_AMPA_on_E}_dt${dt}";
#$FN=$REG_TAG."EEEIIEII_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";

#$FN=$REG_TAG."EEEIIEII_tauREEx${tau_r_AMPA_on_E}_tauDEEx${tau_d_GABA_on_E}_tauRIIx${tau_r_GABA_on_I}_tauDIIx${tau_d_GABA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaIex${sigma_Ie_E}_sigmaWNE${sigma_white_noise_E}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_sigmaWNEx${sigma_white_noise_E}_sigmaIex${sigma_Ie_E}_IappE${Iapp_E}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_sigmaWNEx${sigma_white_noise_E}_sigmaWNIx${sigma_white_noise_I}_IappE${Iapp_E}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_lamEEx${lambda_AMPA_on_E}_lamEIx${lambda_AMPA_on_I}_lamIEx${lambda_GABA_on_E}_lamIIx${lambda_GABA_on_I}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_pIIx${p_II}_pIIGJx${p_II_gap_junction}_lamEE${lambda_AMPA_on_E}_lamEI${lambda_AMPA_on_I}_lamIE${lambda_GABA_on_E}_lamII${lambda_GABA_on_I}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappIx${Iapp_I}_pIIGJxi${h5}_${begin_Val31}_${end_Val31}_${len_Val31}_pII${p_II}_lamEE${lambda_AMPA_on_E}_lamEI${lambda_AMPA_on_I}_lamIE${lambda_GABA_on_E}_lamII${lambda_GABA_on_I}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_gIE${g_syn_GABA_on_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_gEI${g_syn_AMPA_on_I}_gIE${g_syn_GABA_on_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_gEIx${g_syn_AMPA_on_I}_gIEx${g_syn_GABA_on_E}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_pEIx${p_EI}_pIEx${p_IE}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauRIE${tau_r_GABA_on_E}_tauDIE${tau_d_GABA_on_E}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_gEIx${g_syn_AMPA_on_I}_gIEx${g_syn_GABA_on_E}_pEIx${p_EI}_pIEx${p_IE}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_pEE${p_EE}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_sigmaWNIx${sigma_white_noise_I}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauLII${tau_l_GABA_on_I}_tauRII${tau_r_GABA_on_I}_tauDII${tau_d_GABA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauLIE${tau_l_GABA_on_E}_tauRIE${tau_r_GABA_on_E}_tauDIE${tau_d_GABA_on_E}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauLEI${tau_l_AMPA_on_I}_tauREI${tau_r_AMPA_on_I}_tauDEI${tau_d_AMPA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauLEE${tau_l_AMPA_on_E}_tauREE${tau_r_AMPA_on_E}_tauDEE${tau_d_AMPA_on_E}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_wb_g_L${wb_g_L}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauLIEx${tau_l_GABA_on_E}_tauRIEx${tau_r_GABA_on_E}_tauDIEx${tau_d_GABA_on_E}_tauLIIx${tau_l_GABA_on_I}_tauRIIx${tau_r_GABA_on_I}_tauDIIx${tau_d_GABA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tLEEx${tau_l_AMPA_on_E}_tREEx${tau_r_AMPA_on_E}_tDEEx${tau_d_AMPA_on_E}_tLEIx${tau_l_AMPA_on_I}_tREIx${tau_r_AMPA_on_I}_tDEIx${tau_d_AMPA_on_I}_tLIEx${tau_l_GABA_on_E}_tRIEx${tau_r_GABA_on_E}_tDIEx${tau_d_GABA_on_E}_tLIIx${tau_l_GABA_on_I}_tRIIx${tau_r_GABA_on_I}_tDIIx${tau_d_GABA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauRIIx${tau_r_GABA_on_I}_tauDIIx${tau_d_GABA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";
#$FN=$REG_TAG."EEEIIEII_tauRIEx${tau_r_GABA_on_E}_tauDIEx${tau_d_GABA_on_E}_tauRIIx${tau_r_GABA_on_I}_tauDIIx${tau_d_GABA_on_I}_IappE${Iapp_E}_sigmaWNE${sigma_white_noise_E}_pII${p_II}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ne${N_e}_Ni${N_i}";

#$FN=$REG_TAG."BW_EEEIIEII_IappE${Iapp_E}_IappIx${Iapp_I}_sigmaWNE${sigma_white_noise_E}_sigmaWNI${sigma_white_noise_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEEIIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappIx${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappIx${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEEIIEII_tauRx${tau_r_AMPA_on_I}_IappI${Iapp_I}_IappE${Iapp_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_tauRx${tau_r_AMPA_on_I}_IappIx${Iapp_I}_IappE${Iapp_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEEIIEII_gEIx${g_syn_AMPA_on_I}_IappIx${Iapp_I}_IappE${Iapp_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEEIIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappI${Iapp_I}_supraE_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappIx${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEIEII_IappExi${u}_${begin_Val13}_${end_Val13}_${len_Val13}_IappIx${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEEIIEII_IappEx${Iapp_E}_IappI${Iapp_I}_supraE_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_IappEx${Iapp_E}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_IappEx${Iapp_E}_IappIx${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEIEII_IappEx${Iapp_E}_IappIx${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEIEII_IappE${Iapp_E}_IappIxi${i}_${begin_Val1}_${end_Val1}_${len_Val1}_sigmaWNE${sigma_white_noise_E}_sigmaWNI${sigma_white_noise_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";
#$FN=$REG_TAG."BW_EEEIIEII_IappE${Iapp_E}_IappIxi${i}_${begin_Val1}_${end_Val1}_${len_Val1}_sigmaWNE${sigma_white_noise_E}_sigmaWNI${sigma_white_noise_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEEIIEII_IappEx${Iapp_E}_sigmaWNEx${sigma_white_noise_E}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."BW_EEIEII_gIIxi${j}_gIIGJx${g_syn_II_gap_junction}_${begin_Val2}_${end_Val2}_${len_Val2}_IappI${Iapp_I}_IappE${Iapp_E}_sigmaWNI${sigma_white_noise_I}_sigmaWNE${sigma_white_noise_E}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}_Ne${N_e}";

#$FN=$REG_TAG."II_sigmaIeIx${sigma_Ie_I}_IappI${Iapp_I}_gII${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_sigmaWNIx${sigma_white_noise_I}_IappI${Iapp_I}_gII${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIx${Iapp_I}_gIIx${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_gIIxi${j}_gIIGJx${g_syn_II_gap_junction}_${begin_Val2}_${end_Val2}_${len_Val2}_pIIGJ${p_II_gap_junction}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_gIIx${g_syn_GABA_on_I}_pIIx${p_II}_IappIx${Iapp_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_gIIx${g_syn_GABA_on_I}_gIIGJx${g_syn_II_gap_junction}_pII${p_II}_pIIGJ${p_II_gap_junction}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_tEnd${tEnd}_gII${g_syn_GABA_on_I}_pII${p_II}_IappI${Iapp_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_Ni${N_i}";
#$FN=$REG_TAG."II_tEnd${tEnd}_gII${g_syn_GABA_on_I}_pII${p_II}_IappI${Iapp_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_Ni${N_i}";
#$FN=$REG_TAG."II_pIIxi${k}_gIIGJx${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_gII${g_syn_GABA_on_I}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_pIIxi${k}_pIIGJxi${h5}_gIIGJ${g_syn_II_gap_junction}_gII${g_syn_GABA_on_I}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIx${Iapp_I}_pII${p_II}_gII${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_pIIxi${k}_pIIGJxi${h5}_${begin_Val3}_${end_Val3}_${len_Val3}_${begin_Val31}_${end_Val31}_${len_Val31}_gIIGJ${g_syn_II_gap_junction}_gII${g_syn_GABA_on_I}_IappI${Iapp_I}_rndV${rnd_init_V}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_lamEIx${lambda_AMPA_on_I}_lamIIx${lambda_GABA_on_I}_pIIGJxi${h5}_${begin_Val31}_${end_Val31}_${len_Val31}_pII${p_II}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIx${Iapp_I}_pIIGJxi${h5}_${begin_Val31}_${end_Val31}_${len_Val31}_pII${p_II}_lamEI${lambda_AMPA_on_I}_lamII${lambda_GABA_on_I}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIx${Iapp_I}_pIIGJxi${h5}_${begin_Val31}_${end_Val31}_${len_Val31}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIxi${i}_pIIGJxi${h5}_${begin_Val1}_${end_Val1}_${len_Val1}_${begin_Val31}_${end_Val31}_${len_Val31}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIxi${i}_${begin_Val1}_${end_Val1}_${len_Val1}_pII${p_II}_gII${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ni${N_i}";
#$FN=$REG_TAG."II_IappIx${Iapp_I}_pII${p_II}_gII${g_syn_GABA_on_I}_gIIGJ${g_syn_II_gap_junction}_pIIGJ${p_II_gap_junction}_dt${dt}_tEnd${tEnd}_Ni${N_i}";

######### Reserve

####################### Derived Parameters #######################
	$Nt=int($tEnd/$dt);			# The number of rounds of simulation, i.e. tEnd=dt*Nt

	$abs_file_name_Str = $store_dir.$FN;
	$abs_file_name_Str1 = $store_dir.$FN."_Traj";
	$abs_file_name_Str2 = $store_dir.$FN."_LEs";
	$abs_file_name_Str3 = $store_dir.$FN."_LEs_End_Pert";
	$abs_file_name_Str4 = $store_dir.$FN."_LEs_Pert";
	$abs_file_name_Str5 = $store_dir.$FN."_LEs_End";

	$done="DONE".$FN;

	# Typical synaptic properties #
	$tau_m_on_E=$tau_m_e;     	# [ms]
	$tau_m_on_I=$tau_m_i;     	# [ms]

	# Typical external input properties #
	$v_ext_syn_AMPA=$v_syn_AMPA;
	$v_ext_syn_GABA=$v_syn_GABA;
	
	$tau_m_ext_on_E=$tau_m_e;   # [ms]
	$tau_m_ext_on_I=$tau_m_i;   # [ms]

	$g_ext_AMPA_on_E=$g_syn_AMPA_on_E;   	# [nS]
	$g_ext_AMPA_on_I=$g_syn_AMPA_on_I;   	# [nS]
	$g_ext_GABA_on_E=$g_syn_GABA_on_E;   	# [nS]
	$g_ext_GABA_on_I=$g_syn_GABA_on_I;   	# [nS]
	
	$tau_l_ext_AMPA_on_E=$tau_l_AMPA_on_E;
	$tau_l_ext_AMPA_on_I=$tau_l_AMPA_on_I;
	$tau_l_ext_GABA_on_E=$tau_l_GABA_on_E;
	$tau_l_ext_GABA_on_I=$tau_l_GABA_on_I;	
	
	$tau_r_ext_AMPA_on_E=$tau_r_AMPA_on_E;
	$tau_r_ext_AMPA_on_I=$tau_r_AMPA_on_I;
	$tau_r_ext_GABA_on_E=$tau_r_GABA_on_E;
	$tau_r_ext_GABA_on_I=$tau_r_GABA_on_I;	
	
	$tau_d_ext_AMPA_on_E=$tau_d_AMPA_on_E;
	$tau_d_ext_AMPA_on_I=$tau_d_AMPA_on_I;
	$tau_d_ext_GABA_on_E=$tau_d_GABA_on_E;
	$tau_d_ext_GABA_on_I=$tau_d_GABA_on_I;
		
####################### Send commands #######################
	
	$MYFLAG_STR=$abs_file_name_Str." ".$abs_file_name_Str1." ".$abs_file_name_Str2." ".$abs_file_name_Str3." ".$abs_file_name_Str4." ".$abs_file_name_Str5." ".$seed." ".$isUsingSeed." ".$N_i." ".$N_e." ".$N_ext." ".$p_EE." ".$p_EI." ".$p_IE." ".$p_II." ".$V_e_rest." ".$V_e_th." ".$V_e_reset." ".$V_e_refrac." ".$tau_m_e." ".$C_m_e." ".$V_i_rest." ".$V_i_th." ".$V_i_reset." ".$V_i_refrac." ".$tau_m_i." ".$C_m_i." ".$g_syn_AMPA_on_E." ".$g_syn_AMPA_on_I." ".$g_syn_GABA_on_E." ".$g_syn_GABA_on_I." ".$v_syn_AMPA." ".$v_syn_GABA." ".$tau_m_on_E." ".$tau_m_on_I." ".$tau_l_AMPA_on_E." ".$tau_l_AMPA_on_I." ".$tau_l_GABA_on_E." ".$tau_l_GABA_on_I." ".$tau_r_AMPA_on_E." ".$tau_r_AMPA_on_I." ".$tau_r_GABA_on_E." ".$tau_r_GABA_on_I." ".$tau_d_AMPA_on_E." ".$tau_d_AMPA_on_I." ".$tau_d_GABA_on_E." ".$tau_d_GABA_on_I." ".$lambda_AMPA_on_E." ".$lambda_AMPA_on_I." ".$lambda_GABA_on_E." ".$lambda_GABA_on_I." ".$v_ext_syn_AMPA." ".$v_ext_syn_GABA." ".$g_ext_AMPA_on_E." ".$g_ext_AMPA_on_I." ".$g_ext_GABA_on_E." ".$g_ext_GABA_on_I." ".$tau_m_ext_on_E." ".$tau_m_ext_on_I." ".$tau_l_ext_AMPA_on_E." ".$tau_l_ext_AMPA_on_I." ".$tau_l_ext_GABA_on_E." ".$tau_l_ext_GABA_on_I." ".$tau_r_ext_AMPA_on_E." ".$tau_r_ext_AMPA_on_I." ".$tau_r_ext_GABA_on_E." ".$tau_r_ext_GABA_on_I." ".$tau_d_ext_AMPA_on_E." ".$tau_d_ext_AMPA_on_I." ".$tau_d_ext_GABA_on_E." ".$tau_d_ext_GABA_on_I." ".$RmIe_E." ".$RmIe_I." ".$dt." ".$Nt." ".$done." ".$is_turn_on_pert." ".$pert_size." ".$show_each_second." ".$is_rec_spkp." ".$is_rec_traj." ".$is_rec_LEs." ".$is_rec_LEs_End_Pert." ".$is_rec_LEs_Pert." ".$is_rec_LEs_End." ".$is_hetero_Ie_E." ".$is_hetero_Ie_I." ".$sigma_Ie_E." ".$sigma_Ie_I." ".$rnd_init_perts_seed." ".$rnd_init_V." ".$rnd_init_RmIe_seed." ".$rnd_determineConnections." ".$rnd_duringExecution." ".$sigma_white_noise_E." ".$sigma_white_noise_I." ".$rnd_white_noise_E." ".$rnd_white_noise_I." ".$is_use_GSL_2_solve_ODE." ".$type_of_external_inputs." ".$p_EE_gap_junction." ".$p_EI_gap_junction." ".$p_IE_gap_junction." ".$p_II_gap_junction." ".$g_syn_II_gap_junction." ".$is_rand_V_init_X." ".$wb_g_L." ".$pyr_cell." ".$allow_EE_connections." ".$allow_II_connections." ".$isGJExactNumberOfConns;

	if ($blast_job == 1)
	{		
		system("echo \"[`date`] Submit $FN\"");
		$cmd="qsub -v MYFLAG=\"$MYFLAG_STR\" $mjob_file";
		system($cmd);
	}
	else
	{
		system("echo \"[`date`]  run  $FN\"");
		
		if ($wait2see == 1)
		{
			$cmd="./hh_network_PING_ING_BW \"$MYFLAG_STR\" ";
			system($cmd);
		}
		else
		{
			$cmd="nohup nice -n 19 ./hh_network_PING_ING_BW \"$MYFLAG_STR\" >nohup\"$done\".out 2>&1 &";
			system($cmd);
		}
				
	}
	

	sleep(0.1);
}
}
}
}
}
}
}
}
}				
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}

system("echo \"[`date`] Stop script: Success\"");