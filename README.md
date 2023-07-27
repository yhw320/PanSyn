# PanSyn

--------------------------
Table of Contents
--------------------------
> 1. Introduction
> 2. System requirements
> 3. Installation guide
> 4. Demo
> 5. Instructions for use


## 1. Introduction
Based on the published algorithms or tools developed by our and other groups, we introduce a detailed protocol for the most comprehensive and up-to-date genome synteny pipeline (called PanSyn) and provides step-by-step instructions as well as application examples for demonstrating how to use it. PanSyn pipeline includes three major modules (microsynteny analysis, macrosynteny analysis, and integrated micro & macro analysis). 


## 2. System requirements
A computer cluster with the Unix-based operating system (e.g., Linux, Ubuntu, CentOS) is needed, and the demands for memory and disk space heavily depend on the number of species for analysis. We conducted the analysis of demo data (containing about 20 species) under 32 GB of random-access memory and ~100 GB of free disk space. 

Basic Dependencies: Perl 5.16.3.

Packages | Version used in Research|
---------| --------|
SVG      | 2.86    |

To run PanSyn, users need to install the following third-party software and add them to the PATH environment variables. We advise installing particular software when users run the program, because various analysis modules need different software. All the software versions description are described in the protocol [Software prerequisites]. 


## 3. Installation guide
Source code is freely available at https://github.com/yhw320/PanSyn.
Download the software package and save it to the specified path. PanSyn is an executable file written in Perl, and users can run it directly.
The time from download to installation is less than 10 minutes.


## 4. Demo
We provide the simplified demo datasets and full demo datasets in the website (https://doi.org/10.5281/zenodo.7588262), respectively.
Simplified demo datasets: including input files and main result files.
Full demo datasets: including input files, all processed data and result files.
Check out our user protocol [timing] for more information about the expected run time.

## 5. Instructions for use
### Prepare input files
Please refer to the [inputDir] folder in demo datasets(https://doi.org/10.5281/zenodo.7588262) to generate the input files.
### PanSyn run
PanSyn has multiple subroutines (Microsyteny, Macrosyteny and Integrated Micro & Macro analysis). Users only need to prepare input files, simply modify the configuration file and corresponding command parameters to execute them. 
Detailed function description and parameter settings are described in the protocol [Procedure]. 
### Quick start
```
perl PanSyn/Syn_scripts/*.pl -parameter
```
## Contact
If you have any questions, please feel free to contact: liyuli@ouc.edu.cn
