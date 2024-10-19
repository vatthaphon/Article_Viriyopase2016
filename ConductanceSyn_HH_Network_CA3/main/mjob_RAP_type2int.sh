#!/bin/sh
# Script to submit a job

# Set the interpreting shell
#$ -S /bin/bash

# Execute the job from the current working directory
#$ -cwd

# Set the executable file
LD_LIBRARY_PATH=/home/viriyopa/utils/nlopt/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/utils/glpk/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/utils/gsl/lib:$LD_LIBRARY_PATH

LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Designing_Network/main:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Designing_Network/misc:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Designing_Network_Sim_Network/main:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Designing_Network_Sim_Network/misc:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Lyap_Exp/main:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Lyap_Exp/misc:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Lyap_Exp_Ext_Synapse/main:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Lyap_Exp_Ext_Synapse/misc:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Lyap_Exp_Ext_Synapse_PC/main:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/viriyopa/source_code/paper2_raoul/Lyap_Exp_Ext_Synapse_PC/misc:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH

#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_CA3
#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_CA3_NWB
export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_BW
#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_BW_PulseSyn
#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_BW_incdec_I
#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_WB

# Set the name of the job
#$ -N z_tmp_Output

# Merge the standard error and standard output streams
#$ -j y

# Mail to user at beginning/end/abort/on suspension
# -m beas

# Mail to this email, rather than submitting user
#$ -M a.viriyopase@donders.ru.nl

# Set priority
#$ -p -500

# -q '*@@snnhosts-4'
# -q '*@@snnhosts-12'
# -q '*@@neuroinfhosts-8'
#$ -q '*@@neuroinfhosts-64'
# -q '*@@neuroinfhosts'

$MYAPP "$MYFLAG"

