# PanSyn

--------------------------
Table of Contents
--------------------------
> 1. Introduction
> 2. System requirements
> 3. Installation guide
> 4. Demo datasets
> 5. Instructions for use


## 1. Introduction
Based on the published algorithms or tools developed by our and other groups, we introduce a detailed protocol for the most comprehensive and up-to-date genome synteny pipeline (called PanSyn) and provide step-by-step instructions as well as application examples for demonstrating how to use it. PanSyn pipeline includes three major modules (microsynteny analysis, macrosynteny analysis, and integrated micro & macro synteny analysis). 


## 2. System requirements
### Hardware：
A computer cluster running Linux (e.g., CentOS, Ubuntu) with the demands for random-access memory (RAM) and free disk space heavily depend on the number of species for analysis. We conducted the demo data (containing about 20 species) under 32 GB of RAM and ~100 GB of free disk space. Internet connection is required for downloading and installing PanSyn from Conda or Docker

### Software prerequisites：
We have packaged PanSyn with all its dependencies as one Conda package and made a Docker image of PanSyn with all needed programs and dependencies.
Three ways to install PanSyn. 
**To install PanSyn from Conda or Docker**, make sure that you have preinstall conda or Docker. The installation was tested on Conda version 23.7.1 and the Docker version 1.13.1, build 7d71120/1.13.1. 
**To install PanSyn from GitHub**, make sure that you have preinstall the following dependencies and add to the PATH environment variable. Refer to the list of dependencies in the protocol [Software prerequisites], that indicate which module(s) use which packages.


## 3. Installation guide
### Three ways to install PanSyn:
#### Install from Conda:
    (1) To add channels in conda, you can use the commands:
	`conda config --add channels bioconda`
	`conda config --add channels conda-forge`
	'conda config --add channels seqera'
	'conda config --add channels dnachun'
       
    (2) Verify that the channels have been added successfully.
	'conda config --show channels'
       
    (3) Create an environment named pansyn and active it.
	'conda create --name pansyn'
	'conda activate pansyn'
       
    (4) Install PanSyn
	'conda install -c micromacro pansyn -y'

#### Install from Docker:
    (1) Pull image from Dockerhub.
	'docker pull micromacro/pansyn:last'

    (2) Mount local files into docker container. Replace <local_dir_path> with your local dir path.
	'docker run -it -v <local_dir_path>:/root/workspace/ micromacro/pansyn:last /bin/bash'

    (3) Activate the environment and script.
	'source activate Pansyn'
      	'source /opt/conda/envs/Pansyn/cns_solve_1.3/cns_solve_env.sh'

#### Install from GitHub:
    (1) Download and unpack https://github.com/yhw320/PanSyn/archive/refs/heads/main.zip. Or using the following command.
	'git clone https://github.com/yhw320/PanSyn/archive/refs/heads/main.zip'

    (2) PanSyn software package includes scripts located in the directory "scripts" that users can run directly without compilation.
	'cd scripts'
	'perl *.pl'
       
### Timing: 
Installing PanSyn via Conda or Docker can be completed within approximately 30 minutes, but it may take longer depending your internet speed.

## 4. Demo datasets
We provide the simplified demo datasets and full demo datasets in the website (https://doi.org/10.5281/zenodo.7588262), respectively.
Simplified demo datasets: including input files and main result files.
Full demo datasets: including input files, all processed data and result files.

Check out our user protocol [timing] for more information about the expected run time.

## 5. Instructions for use
### Prepare input files
Please refer to the [inputDir] folder in demo datasets(https://doi.org/10.5281/zenodo.7588262) to prepare input files.

### How to run PanSyn
PanSyn has multiple subroutines (Microsyteny, Macrosyteny and Integrated Micro & Macro synteny analysis). Users only need to prepare input files and corresponding command parameters to execute them. 
Detailed function description and parameter settings are described in the protocol [Procedure]. 

### Contact
If you have any questions, please feel free to contact: liyuli@ouc.edu.cn or hongweiyu@stu.ouc.edu.cn
