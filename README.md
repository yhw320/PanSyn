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
A computer cluster running Linux (e.g., CentOS, Ubuntu) with the demands for random-access memory (RAM) and free disk space depend on the number of species for analysis. The demo data (containing about 20 species) were successfully run under 32 GB of RAM and ~100 GB of free disk space. Internet connection is required for downloading and installing PanSyn from Conda or Docker.

### Software prerequisites：
We have packaged PanSyn with all its dependencies as one Conda package and made a Docker image of PanSyn with all needed programs and dependencies.<br>Three ways to install PanSyn.<br>**To install PanSyn from Conda or Docker**, make sure that you have preinstall conda or Docker. The installation was tested on Conda version 23.7.1 and the Docker version 1.13.1, build 7d71120/1.13.1.<br>**To install PanSyn from GitHub**, make sure that you have preinstall the dependencies and add to the PATH environment variable. Refer to the list of dependencies in the protocol [Software prerequisites], that indicate which module(s) use which packages.


## 3. Installation guide
### Three ways to install PanSyn:
#### A. Install from Conda:
(1) To add channels in conda, you can use the commands:<br>
```
conda config --add channels seqera
conda config --add channels dnachun
conda config --add channels bioconda
conda config --add channels r
conda config --add channels conda-forge
```
 
(2) Verify that the channels have been added successfully.<br>
```
conda config --show channels  
```
(3) Create an environment named pansyn and active it.<br>
```
conda create --name pansyn  
conda activate pansyn  
```
(4) Install PanSyn.<br>
```
conda install -c micromacro pansyn -y
```
**CRITICAL**: If the Step 3A(4) takes a long time (>30 min), users can accelerate the installation of PanSyn using Conda’s mamba. Replace Steps 3A(3-4) with the following command lines:<br>
 ```
	conda create -n pansyn python=3.10.12 -y
	conda activate pansyn
	conda install mamba -c conda-forge -y
	mamba install -c micromacro pansyn -y
 ```
#### B. Install from Docker:
(1) Pull image from Dockerhub.<br>
```
docker pull micromacro/pansyn:latest
```
(2) Mount local files into docker container.<br>
```
docker run -it -v <your_host_path>:<your_container_path> micromacro/pansyn:latest /bin/bash
```
(3) Activate the environment and script.<br>
```
source activate Pansyn
source /opt/conda/envs/Pansyn/cns_solve_1.3/cns_solve_env.sh
```
#### C. Install from GitHub:
(1) Download and unpack https://github.com/yhw320/PanSyn/archive/refs/heads/main.zip. Or using the following command.<br>
```
git clone https://github.com/yhw320/PanSyn.git  
```
(2) PanSyn package includes scripts located in the directory "scripts" which users can run directly without compilation.<br>
```
cd PanSyn/scripts  
perl *.pl
```
### Timing: 
~30 min (using Conda), ~15 min (using Conda’s mamba or Docker), or ~2 h (using GitHub). But it may take longer depending your internet speed.

## 4. Demo datasets
We provide the full demo datasets in the website (https://doi.org/10.5281/zenodo.10115240), which includes input files, all processed data and result files.

Check out our user protocol [timing] for more information about the expected run time.

## 5. Instructions for use
### Prepare input files
Please refer to the [inputDir] folder in demo datasets(https://doi.org/10.5281/zenodo.10115240) to prepare input files.

### How to run PanSyn
PanSyn has multiple subroutines (Microsyteny, Macrosyteny and Integrated Micro & Macro synteny analysis). Users only need to prepare input files and corresponding command parameters to execute them. 
Detailed function description and parameter settings are described in the protocol [Procedure]. 

Note: The PanSyn pipeline can work on both a stand-alone server and a cluster. We have provided instructions and scripts for running PanSyn in a cluster environment (e.g., with PBS or SLURM job schedulers; see https://github.com/yhw320/PanSyn/tree/main/SchedulerScripts).

### Contact
If you have any questions, please feel free to contact: liyuli@ouc.edu.cn or hongweiyu@stu.ouc.edu.cn
