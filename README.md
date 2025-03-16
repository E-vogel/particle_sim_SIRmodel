# Contagion Simulation Using Particle-Based Model and SIR Model

## Overview
This project simulates the spread of an infection using a particle-based model and compares it with the classical SIR (Susceptible-Infected-Recovered) model. The simulation is implemented in MATLAB and visualizes how an infection spreads through a population over time.

## Purpose of the Simulation
The goal of this simulation is to:
- Model the transmission of an infectious disease in a population using a particle-based approach.
- Compare the results with the SIR model to observe similarities and differences.
- Provide a visual representation of infection dynamics.

## SIR Model
The SIR model is a mathematical model used to describe the spread of infectious diseases. It consists of three compartments:
- **S (Susceptible):** Individuals who can be infected.
- **I (Infected):** Individuals who are currently infected and can transmit the disease.
- **R (Recovered):** Individuals who have recovered from the infection and are now immune.

The model is governed by the following differential equations:
```math
\frac{dS}{dt} = -\beta S I
```
```math
\frac{dI}{dt} = \beta S I - \gamma I
```
```math
\frac{dR}{dt} = \gamma I
```
where:
- `\beta` is the infection rate.
- `\gamma` is the recovery rate.

## How to Run the Simulation
1. Ensure that you have MATLAB installed.
2. Run the provided MATLAB script.
3. The simulation will generate a visualization comparing the particle-based approach with the SIR model.
4. The results will be saved as a video file (`contagion_sim.avi`).

## Output
![Image](https://github.com/user-attachments/assets/dc8686fd-c62c-4958-b3a8-ca1c9ae7589e)