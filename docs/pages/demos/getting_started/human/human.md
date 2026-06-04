---
layout: default
title: Human
parent: Getting started
grand_parent: Demos
has_children: false
nav_order: 1
---

# Human

## Overview

This page explains how to simulate 20 seconds of human cardiovascular function.

Important notes
{: .label .label-blue}

* FiberVent will <i>only</i> run on Windows PCs. 
* This simulation will take a few minutes to run on a modern PC.
* The manuscript summarizes results from hundreds of simulations. Generating the complete set of simulations could take a week on a laptop.
* The simulations can be generated using files from this repository and a Python environment that can be created via [Conda](https://conda-forge.org/download/).

## Instructions

### Setup to run simulations

+ Clone the [FiberVent](https://github.com/Campbell-Muscle-Lab/FiberVent) to a local directory that will be referred to as `<FiberVent_repo>`
+ Install [Conda](https://conda-forge.org/download/) if not already available
+ Install the FiberVent environment by
  + opening a terminal and changing the directory to `<FiberVent_repo>/code/FiberVentPy/environment`
  + run `conda env create -f environment.yml`

### Run simulations

Warning - this simulation takes several minutes to run.
{: .label .label-blue}

+ Open a terminal
+ Activate the conda environment with the command `conda activate FiberVent`
+ Change the working directory to `<FiberVent_repo>/code/FiberVentPy/FiberVentPy`
+ Run the command
  + `python fiberventcpy.py characterize ../../../../demos/getting_started/human/setup/setup_human_1.json`


## Output

### Figure 1
<img src="figures/fig_example_beats.svg">

