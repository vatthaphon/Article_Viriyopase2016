#!/bin/sh
# Script to submit an MPI job

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
#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_BW
#export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_WB
export MYAPP=/home/viriyopa/source_code/paper2_raoul/ConductanceSyn_HH_Network_CA3/main/hh_network_PING_ING_WHITE_WB

# Set the name of the job
#$ -N z_tmp_Output

# Merge the standard error and standard output streams
#$ -j y

# Mail to user at beginning/end/abort/on suspension
# -m beas

# Mail to this email, rather than submitting user
#$ -M a.viriyopase@donders.ru.nl

# Add necessary modules
. /etc/profile.d/modules.sh
module add default-environment sge/6.2 mvapich2/gcc/64/1.2

#########################################################################################
#priority              20
#$ -q verylong.q
#$ -p -1000

#priority              15
# -q long.q
# -p -800

#priority              10
# -q short.q
# -p -600

#priority              5
# -q veryshort.q
# -p -400

#priority              20
# -q all.q@@hosts_001_010

$MYAPP "$MYFLAG"
#########################################################################################

