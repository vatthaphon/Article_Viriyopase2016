# Cooperation and competition of gamma oscillation mechanisms

<p align="justify">When ING and PING co-exist, how do they interact?</p>

<p align="justify">This repository contains the source code accompanying the paper titled "Cooperation and competition of gamma oscillation mechanisms", publised at Journal of Neurophysiology. 
The goal of this study is to investigate how ING and PING interact in a network that both co-exist.</p>

## Table of Contents  
[1. Introduction](#Introduction)  
[2. Installation](#Installation)  
[3. Usage](#Usage)  
[4. File Structure](#FileStructure)  
[5. License](#License)  
          
## Introduction<a name="Introduction"/>
<p align="justify">This repository provides the implementation of the methods and algorithms described in the paper "Cooperation and competition of gamma oscillation mechanisms". 
The main goal of this study is to numerically investigate how ING and PING interact when both co-exist in a network. We generally find that the network assumes either ING or PING depending
on which mechanisms have higher intrinsic frequency.</p>

You can find the full paper [here] (To be filled).
  
## Installation<a name="Installation"/>
- Matlab 2015.
- GNU C compiler.

## Usage<a name="Usage"/>
- Run the files in /ConductanceSyn_HH_Network_CA3 to simulate network dynamics. Store simulation results in corresponding /Data directory.
- Run the Matlab files to generate each figure in the paper.

## File Structure<a name="FileStructure"/>
├── ConductanceSyn_HH_Network_CA3/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Contains files for cleaning the tabular data.  
├── Fig1/&nbsp;&nbsp;&nbsp;&nbsp;# Contains files for transmitting EEG data from participants' PCs  to the server.  
├── EEG data analysis/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Contains files for cleaning and analyzing EEG data.  
├── EEG sanity checks/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Contains files for sanity checks.  
├── README.md  
└── LICENSE


## License<a name="License"/>
This project is licensed under the MIT License - see the LICENSE file for details.

