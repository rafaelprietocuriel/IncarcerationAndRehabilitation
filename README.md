# Incarceration And Rehabilitation
Code created in November 2025 by Rafael Prieto-Curiel.

This repository contains the code to reproduce the analysis of incarceration and rehabilitation. 
The software required to compile the project is R (RStudio). 
The code creates and saves in the corresponding repository the results of the execution. The process takes roughly one hour for the full analysis. 

# The code is divided into three sections
S1 - Time series analysis of cartels
S2 - Age assignment of agents
S3 - Analysis of the model

# S1 - Time series
The code creates a time series that corresponds to weekly steps of the model. Five time series are created: Number of active members, recruited, incapacitated, killed and retired each time series indicates the number of people for that weekly step

# S2 - Age assignment of agents
Using the time series produced in S1, an age is assigned to each person. The age is updated each week for all active members of the cartel.

# S3 - Analysis of the model
Based on a population estimate of 1.35 million males born in 1990, first get the whole cartels with all population, then distribute those according to some peak 

### considerations
The code is based on functions that can be reused for other cases. Some parameters can be easily adjusted to fit different attributes of the population.

# relevant parameters
19300 recruits in a year
6500 deaths
5700 incapacitations
on date 44562 (01/01/2022) we get 175000
with the current parameters we get 2012 and 2022

