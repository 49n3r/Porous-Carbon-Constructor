# Porous-Carbon-Constructor
Algorithm to construct porous carbon models (use the Jupyter Notebook script, the .py is currently being written to build this into a Streamlit App)

**POROUS CARBON CONSTRUCTOR ALGORITHM**

><span style="color:brown">`Associated Publication:` <br> Atomistic-to-Continuum Modeling of Carbon Foam: A New Approach to Finite Element Simulation <br>
C. Ugwumadu, W. Downs, R. Thapa, R. Olson III, M. Ali, J. Trembly, Y. Al-Majali and D. A. Drabolda </span>


This script creates a porous carbon model with desired porosity and pore distribution, which is then optimized using molecular dynamics simulation


**Authors Information** 

**Name**:: C. Ugwumadu, R. Thapa, and D. A. Drabold

**Affiliation** (at the time of creation):: Ohio University, Athens Ohio, 45701, USA.

**Date** (published Version):: May, 11 2024

**Contact**:: ugwumaduchinonso@gmail.com

><span style="color:red">*This algorithm is freely distributed and can be adapted as the user sees fit, all we ask is that you cite our work.* Thank You!!!</span>

## Parameters to compute Pore distribution from pre-defined porosity

As described in our paper, our method to determine the pore distribution from the pre-defined porosity $\xi_i$ is obtained as:

$$
\xi_i = \sum_{n=1}^N \frac{v_n}{V}
$$

Where $N$, $\nu_n$, and $V$ are the total number of pores, the pore volume sampled from a uniform distribution, and the bounding box volume, respectively.

## The Main Program

>Here we use the functions and the parameters to create the model. <font color=red> Please note that </font> the algorithm may take a long time (up to days) for huge models on a low-performance computer. This code was designed to run on a high-performance computer, as such, not much has been considered for code optimization.

At the end of running this script, 3 files will be created. 

1. A POSCAR (`VASP`) file containing the atoms without the pores
2. A POSCAR file containing the pores without the atoms
3. An extended *.xyz* file containing atoms and pores, for visualization in an atom visualizer like `OVITO`.

<font color=red> Please note that </font> the sizes of the pores in the *.xyz* file have been scaled by half to be visually presentable.

**The equation above** has been implemented in the `poreCreator` function below. Check and adjust the (default) parameter arguments accordingly

