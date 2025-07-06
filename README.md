# massTransferVoF
OpenFoam 11 VoF solver for mass transfer of H2S and O2 (and other species)

massTransferVoF is an OpenFOAM solver describing mass transfer processes in volume of fluid simulations. It is compatible with OpenFOAM 11 and based on imcompressibleVoF.

## Authors
This code has been created by Katharina Teuber at Jade Hochschule (University of Applied Sciences)

## Installation
Place the solver into your OpenFOAM environment (usually USERNAME11/run/) and compile via ./Allwmake. For further instructions on how to compile custom OpenFOAM solvers, please check out the OpenFOAM documentation and discussion forums on cfd-online.

## Usage/Examples
massTransferVOF can be used for different species and mass transfer simulations. It contains hard coded Henry coefficients of H2S and O2 mass transfer, custom Henry coefficients can be applied as well.
Test cases including reference solutions published here include:
- 2D Diffusion of a tracer on a plate
- Mass transfer across a horizontal interface at rest (equilibrium conditions)
- Mass transfer across a horizontal interface at rest (mass transfer over time)

## How to use?
After successful compilation, massTransferVoF can be used for simulations.
Modifications needed compared to simulations with incompressibleVoF:
- Boundary and initial conditions for tracer concentration C need to be defined (see test cases)
- A new file "transportProperties" needs to be placed in the constant-directory. This file should specify: Dwater, Dair, Schmidtnumber, Temp, HenryConstant or species
- Schmidtnumber needs to be defined but will not be considered for laminar simulations
- Temp needs to be defined but will not be considered if you specify the Henry coefficient via HenryConstant
- Examples for species specification (Hcp and Ctemp are already defined as default values, you do not need to specify them):

Custom Henry coefficient:
___
HenryConstant   0.1;
___
H2S Henry coefficient:
___
H2S
{
    // optional: override values
    Hcp   1e-3;     
    Ctemp   2200;
}
___
O2 Henry coefficient:
___
O2
{
    // optional: override values
    Hcp   1.3e-5;     
    Ctemp   1700;
}
___
Start simulations using the command "massTransferVoF".

## References and further reading
Please cite this code using the following publication:
Teuber, K., Dixit, A., Hinkelmann, R. (2025). Modular CFD Solver for Hydrogen Sulfide and Oxygen Transport - Implementation in OpenFOAM 11. Proceedings of the 13th Urban Drainage Modelling Conference, Innsbruck, Austria.

Explanations of each test cases are provided in this publication as well.
