#!/bin/bash

#PBS -N job_name             # Set the job name
#PBS -l nodes=1:ppn=12       # Set the required number of nodes and CPU cores per node
#PBS -l walltime=999:00:00   # Set the maximum runtime for the job
# Optional directives:
#PBS -q queue_name   # Set the queue name
#PBS -o output_file  # Set the output file path
#PBS -e error_file   # Set the error file path

# Change to the work_directory
cd your_path/outputDir1_S17/

# Run your commands or executable
Macrosyn1 -i1 inputDir1_S17/ -o1 outputDir1_S17/ -a diamond
