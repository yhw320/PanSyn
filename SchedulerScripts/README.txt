
PanSyn Pipeline
The PanSyn pipeline can work on both a stand-alone server and a cluster.

Job Scheduling
We recommend using a job scheduler to initiate the execution of commands on computation nodes in a cluster environment. Popular job schedulers include PBS, SLURM, and others. We have provided example scripts, "Pbs.sh" and "Slurm.sh", on the website (https://github.com/yhw320/PanSyn/tree/main/SchedulerScripts). 

Please note that the above examples provide a basic framework for PBS/SLURM scripts and may require appropriate modifications based on your cluster environment and job requirements. For more detailed information on PBS/SLURM scripts, refer to your cluster documentation or consult with your system administrator.

Submitting Commands with PBS and SLURM
PBS
1. Open a terminal and log in to your cluster node.
2. Retrieve the "Pbs.sh" file and make necessary modifications based on your needs.
3. Run the following command to submit the job to the PBS scheduling system:
   
   >qsub Pbs.sh

SLURM
1. Open a terminal and log in to your cluster node.
2. Retrieve the "Slurm.sh" file and make necessary modifications based on your needs.
3. Run the following command to submit the job to the SLURM scheduling system:
   
   >sbatch Slurm.sh
