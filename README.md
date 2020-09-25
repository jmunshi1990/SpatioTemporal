# Spatio Temporal Code for Covid Spread Dynamics
The code allows you to model epidemic spread such as Covid-19. The current implementation can be run on New York, Pennsylvania states and for individual counties as well. The code can be easily adopted for any input country, state or any demographic region. 

# Citation 
Developed by Joydeep Munshi and Indranil Roy at Lehigh University 

J. Munshi, I. Roy, G. Balasubramanian, Spatiotemporal dynamics in demography-sensitive disease transmission: COVID-19 spread in NY as a case study.
https://arxiv.org/ftp/arxiv/papers/2005/2005.01001.pdf

# Desctiprtion
This program implements a version of Cellular Automata framework
(CA) to investigate spatio-temporal dynamics of disease outbreak 
 sensitive to demographic features such as population density,
 mobility and employment status. The code implement here uses 
 Matlab's vectorization and mapping toolbox.

# Disclaimer                                                          
Different implementations may lead to slightly different
behavour and/or results, but there is nothing wrong with it, 
as this is the nature of random walks and all metaheuristics.

# Features
1. AutomataSpreadModel.m - This script contains all the required functions
to run the spatio-temporal dynamics

2. NYSpopulation.m       - This script contains preprocessing data from
NY state shape files using mapping toolbox. The function provides
normalized population density data for each county and initial cell
states inside a demography. Note: This script can be modiefied with
appropriate shape files for different regions/demography.

3. PASpopulation.m       - This script contains preprocessing data from
PA state shape files using mapping toolbox. The function provides
normalized population density data for each county and initial cell
states inside a demography.

4. script.m 	          - This script contains the sampling of AutomataSpreadModel
script and outputs the end states for a given simulation time.
