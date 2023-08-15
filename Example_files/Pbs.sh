#!/bin/bash

#PBS -N my_job             # Job name
#PBS -l nodes=1:ppn=12     # Number of nodes and CPUs
#PBS -l walltime=999:00:00 # Job walltime
#PBS -q queue_name
#PBS -V
#PBS -p 1023
#PBS -S /bin/bash

# Change to the working directory
cd /your_path/work_directory

# Execute the command
<command_or_script>
