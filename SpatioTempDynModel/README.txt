% -----------------------------------------------------------------
% Cellular Automata (CA) algorithm by Munshi et. al.              %
% Programmed by Joydeep Munshi and Indranil Roy at Lehigh 	  %
% University               					  %
% Programming dates: Mar 2020 - May 2020                          %
% Last revised: May  2020    								      
% -----------------------------------------------------------------
% Papers -- Citation Details:
% J. Munshi, I. Roy, G. Balasubramanian, Spatiotemporal           %
% dynamics in demography-sensitive disease transmission:          %
% COVID-19 spread in NY as a case study.						  
% https://arxiv.org/ftp/arxiv/papers/2005/2005.01001.pdf
% ----------------------------------------------------------------%
% This program implements a version of Cellular Automata framework%
% (CA) to investigate spatio-temporal dynamics of disease outbreak% 
% sensitive to demographic features such as population density,   %
% mobility and employment status. The code implement here uses    %
% Matlab's vectorization and mapping toolbox.                     %

% --------------------------------------------------------------- %
% =============================================================== %
% Notes:                                                          %
% Different implementations may lead to slightly different        %
% behavour and/or results, but there is nothing wrong with it,    %
% as this is the nature of random walks and all metaheuristics.   %
% -----------------------------------------------------------------

Description of Matlab script are mentioned here:

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

---------------------------------------------------------------------

Data availibility:

1. Real-time data of deaths and infections can be found in this repository:
https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series

2. Shape files and map data are made availaible in the folder NYS_Civil_Boundaries.shp for NY state
and PA_Civil_Boundaries.shp for Pennsylvania state.
