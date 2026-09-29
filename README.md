# Noise_Plasticity

#### Manuscript Title: Environmental impacts on gene expression noise and its relationship with fitness 

#### Author : Taslima Haque

#### Last modified: 09/29/2026

#### Clone this github repository 

```
$ git clone https://github.com/tahia/Noise_Plasticity.git
```

#### Directory structure:
```
ParentDirectory:
           Data
           Results
           Plots
           Scripts
```

Data directory has all the raw and initial files required to reproduce the analysis

```
$ ls Data/*
```
#### Scipt Descriptions:
```
1. NoisePlasticitySummary.R    
Script to test for noise plasticity across environments and get summary statistics

2. NoiseExpFitnessRelationship.R
Script to model non-linear expression-fitness relationships and get residual noise

3. PopsimExpFitDeapConstantNoise.py
Run different simulations using DEAP package 

4. SimparseFunctions.R
Script with all the summary functions to parse simulations

5. SimSummary.R
Generate and compare simulation summaries 

6. MakeSimJobs.R
Generate a list of jobs for a given combination of expression, noise level, and fitness function 

```

#### Detailed description of running simulations with PopsimExpFitDeapConstantNoise.py

##### Overview
PopsimExpFitDeapConstantNoise.py simulates competition between two populations, A and B, whose individuals differ in a quantitative expression phenotype.
Each individual has:

An expression value
A population label (A or B)
A fitness determined by the reciprocal of its expression-dependent doubling time

At each simulation step, individuals are selected using fitness-proportionate roulette selection. Offspring expression is determined by the parent’s expression, population-specific baseline parameters, heritability, and a constant environmental-noise component.

The script supports Gaussian, log-normal, and two-component mixture fitness functions. Summary statistics are written to a CSV file for each time point and replicate simulation.

To run multiple parameters (expression, noise, and different fitness functions) simultatnously, please use MakeSimJobs.R. 

##### Requirements

Python 3\
NumPy\
pandas\
SciPy\
DEAP

Install the required packages with:

```
$ pip install numpy pandas scipy deap
```

Basic usage:

```
python PopsimExpFitDeapConstantNoise.py \
    --EXPR_MEAN_A 0.5 \
    --EXPR_SD_A 0.1 \
    --EXPR_MEAN_B 1.0 \
    --EXPR_SD_B 0.1 \
    --FIT_var1 0.8 \
    --FIT_var2 0.2 \
    --h 0.5 \
    --fitness_function gaussian \
    --outfile simulation_results.csv
```
By default, the output file is written to the current directory.

Simulation model
Initial populations
Each replicate begins with two populations:

Population A contains pop_size individuals.\
Population B contains pop_size individuals.

Therefore, the initial combined population size is:
2 × pop_size
\
\
Initial expression values are independently sampled from normal distributions:\
\
Population A: Normal(EXPR_MEAN_A, EXPR_SD_A)\
Population B: Normal(EXPR_MEAN_B, EXPR_SD_B)
\
Offspring expression values are clipped to:
0.01 ≤ expression ≤ 2.0

Fitness
An individual’s fitness is the reciprocal of its predicted doubling time:\
\
fitness = 1 / doubling_time(expression)\
Individuals with shorter doubling times therefore have greater fitness and are more likely to be selected as parents.
Selection is performed using DEAP’s roulette-wheel selection:
tools.selRoulette

Offspring expression
For a mother with expression xm, offspring expression is calculated as:\
xd=μ+ sqrt(h^2)(xm−μ)+ sqrt((1−h^2))ϵ

where:

xd is the daughter’s expression.\
xm is the mother’s expression.\
μ is the initial expression mean for the mother’s population.\
h is the heritability parameter.\
ϵ∼N(0,σ).\
σ is the initial expression SD for the mother’s population.

The environmental-noise SD remains constant throughout the simulation. Population labels are inherited directly from the mother.
