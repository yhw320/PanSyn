#!/bin/bash
#SBATCH --job-name=job_name  # Set the job name
#SBATCH --nodes=1  # Set the required number of nodes
#SBATCH --ntasks-per-node=4  # Set the required number of CPU cores per node
#SBATCH --time=01:00:00  # Set the maximum runtime for the job

# Optional directives:
#SBATCH --partition=queue_name  # Set the queue name
#SBATCH --output=output_file  # Set the output file path
#SBATCH --error=error_file  # Set the error file path

# Change to the work_directory
cd your_path/outputDir1_S17/

# Run your commands or executable
Macrosyn1 -i1 inputDir1_S17/ -o1 outputDir1_S17/ -a diamond
