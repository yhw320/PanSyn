#!/bin/bash
#SBATCH --job-name=my_job           # Job name
#SBATCH --output=output.log         # Output log file
#SBATCH --error=error.log           # Error log file
#SBATCH --partition=partition_name  # Partition name
#SBATCH --nodes=1                   # Number of nodes to use
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=1           # Number of CPUs per task
#SBATCH --mem=1G                    # Memory limit per node
#SBATCH --time=999:00:00             # Job run time limit

# Execute command or script
<command_or_script>