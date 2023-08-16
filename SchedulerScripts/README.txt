
PanSyn Pipeline
The PanSyn pipeline can work on both a stand-alone server and a cluster.

Job Scheduling
We recommend using a job scheduler to initiate the execution of commands on computation nodes in a cluster environment, including PBS, SLURM, and others. We have provided example scripts, "Pbs.sh" and "Slurm.sh", on the website (https://github.com/yhw320/PanSyn/tree/main/SchedulerScripts). 

Please note that the above examples provide a basic framework for PBS/SLURM scripts and require appropriate modifications based on your cluster environment and job requirements. For more detailed information on PBS/SLURM scripts, refer to your cluster documentation or consult with your system administrator.

To submit a script using PBS (Portable Batch System), you can follow these steps:
1. Prepare a PBS script:
   Modify the "Pbs.sh" script file that contains the commands and parameters needed for your job.
   In the PBS script, you need to set some PBS directives to define the characteristics of your job, such as the 
   job name, the number of nodes, CPU cores, runtime, etc. 

   Here are some common examples of PBS directives:
   #PBS -N job_name             # Set the job name
   #PBS -l nodes=1:ppn=12       # Set the required number of nodes and CPU cores per node
   #PBS -l walltime=999:00:00   # Set the maximum runtime for the job
  
2. Add task commands: 
   Add the actual task commands you want to execute on the compute nodes in your script. For example, if you want 
   to run PanSyn, you can add it to the PBS script:
   >Macrosyn1 -i1 inputDir1_S17/ -o1 outputDir1_S17/ -a diamond

3. Submit the PBS job: 
   Use the qsub command to submit the PBS script. Execute the following command in the terminal:
   >qsub Pbs.sh

4. Check job status: 
   You can use the qstat command to check the status and progress of your jobs.
   >qstat


To submit a script using SLURM (Simple Linux Utility for Resource Management), you can follow these steps:
1. Prepare a SLURM script:
   Modify the "Slurm.sh" file that contains the commands and parameters needed for your job.

2. Set SLURM directives: 
   In the SLURM script, you need to set some SLURM directives to define the characteristics of your job, such as the job name, the number of nodes, CPU cores, runtime, etc. 
   
   Here are some common examples of SLURM directives:
   #SBATCH --job-name=job_name    # Set the job name
   #SBATCH --nodes=1              # Set the required number of nodes
   #SBATCH --ntasks-per-node=4    # Set the required number of CPU cores per node
   #SBATCH --time=01:00:00        # Set the maximum runtime for the job

3. Add task commands: 
   Add the actual task commands you want to execute on the compute nodes in your script. For example, if you want to run PanSyn, you can add it to the SLURM script:
   >Macrosyn1 -i1 inputDir1_S17/ -o1 outputDir1_S17/ -a diamond

4. Submit the SLURM job: 
   Use the sbatch command to submit the SLURM script. Execute the following command in the terminal:
   >sbatch Slurm.sh

5. Check job status: 
   You can use the squeue command to check the status and progress of your jobs.
   >squeue
