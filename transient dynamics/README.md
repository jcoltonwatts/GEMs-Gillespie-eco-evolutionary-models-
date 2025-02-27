## README for "Evolution reduces the duration of transient coexistence in a simple two-species competition model"

##### Manuscript Authors: J. Colton Watts, Clayton E. Cressler, John P. DeLong README Author: J. Colton Watts Contact:  j.colton.watts\@gmail.com

#### Overview

The archive includes three subdirectories, "data", "RCode", and "MatlabCode". The "data" subdirectory holds all of the .csv files that are generated by running the scripts and functions in the "MatlabCode" subdirectory. The "RCode" subdirectory holds .R scripts that are used to analyze and plot the data contained in the .csv files. The complete contents of each subdirectory and their usage is described in detail below. The order of their description follows the order of the workflow necessary to reproduce the findings presented in the manuscript.

#### MatLab Code

The scripts and functions required to run the Gillespie Eco-evolutionary Model (GEM) are housed in the "MatlabCode" subdirectory. To reproduce the GEM results from our manuscript, all eight files must be present in the user's Current Folder in Matlab. The files and their usage are as follows:\
\
(1) Call_GEMv4_LT_20231025.m -\> This is the main GEM simulation file that sets the parameters and initial conditions for the GEM simulations, calls on the required additional scripts and functions to perform the GEM simulations, and summarizes their results by producing Figures 3 and 4 of the manuscript. It also creates .csv files containing the abundance, trait, and extinction data for all replicate simulations. Call_GEMv4_LT_20231025.m must be run for each initial ecological condition described in the manuscript (below, at, or above the unstable ecological equilibrium), manually changing the initial ecological conditions (y0 values, code lines 82-84) each time. The .csv files generated by are by default stored in the "MatlabCode" subdirectory. The original versions of these .csv files used to estimate extinction rates and create Figure 5 are stored in the "data" subdirectory (described below).\
\
(2) LTs_eqs_and_stability.m -\> This file contains a custom function for finding the location of the unstable equilibrium, should it exist (returned as [R1_hat1, R2_hat1]), and calculating the real part of the dominant eigenvalue of the Jacobian matrix evaluated at this equilibrium (eigen_v_max). Also returned are the isoclines for each competitor as a function of the other competitor's abundance (R1_at_R2, R2_at_R1), as well as the range of competitor abundances over which these isoclines are calculated (R1_Range, R2_Range). This function calls on an additional function, intersections.m, to find the unstable equilibrium using the isoclines.\
\
(3) intersection.m -\> This function is an algorithm for finding the location of the unstable equilibrium given the competitor isoclines.\
\
(4) bdcompmodel.m -\> This function is the deterministic, ordinary differential equation describing the birth-death parameterization of the two-species competition model. It is used to check the ecological dynamics predicted by the deterministic model given the chosen initial ecological states (y0), trait states (b_max1, b_max2), and model parameters (intraspecific competition coefficients b11, b22, d11, d22; intraspecific competition coefficients b12, d12, b21, d21).\
\
(5) GEMv4_bdCompModel_bdTradeoff.m -\> This function performs GEM simulations of the birth-death 2-species competition model given the model parameters (params), the mapping of these parameters to their respective populations (state_parameter_match), the coefficients of variation and heritabilities for the evolving traits (cv_vector and h2_vector, respectively), the desired number of replicates (num_replicates), initial ecological states (y0), and the time span over which the simulations should be run (t_max).\
\
(6) V4_pick_individuals.m -\> This function is called from within the GEM simulation (i.e., in GEMv4_bdCompModel_bdTradeoff) to pull new trait values (analogous to individuals in the GEM framework) from a lognormal distribution given the distribution of the trait in the parental generation.\
\
(7) V4_medians_and_CIs.m -\> A function for generating medians and upper/lower confidence intervals from the time series data generated by the GEM. This function only works for a sequence of standardized times, as generated by the Gillespie algorithm.\
\
(8) jbfill.m -\> This function fills in a space between two vectors with a colored region using MatLab's fill function.

#### Data Files

Running the Matlab code described above produces three .csv files-- one each for the trait values, abundances, and extinction events from the GEM simulations. Each file contains the corresponding values (i.e., traits, abundances, or extinctions) for all replicates across all time steps for a choice of initial ecological conditions and model parameters specified in Call_GEMv4_LT_20231025.m. The naming convention for the files is as follows:\
\
(1) Abundances_W_V_I.csv -\> These files contain the competitors' abundances for each replicate (row) across all time steps (columns) from the GEM simulation. The "W" takes on values of "i" and "j", denoting whether the abundance is of the (arbitrarily chosen) focal population i or competitor population j. The "V" takes on values "j1", "j2", and "j3" corresponding to the three GEM variants, where "j1" = no trait variation, "j2" = trait variation but no heritability (referred to as the "no evolution" scenario in the manuscript), and "j3" = heritable trait variation (referred to as "evolution" scenario). The "I" takes on values "below", "at", and "above", denoting the position of the initial ecological state of the community relative to the unstable coexistence equilibrium.\
\
(2) Extinctions_V_I.csv -\> These files contain the value of the community's "extinction state" (0 = both competitors present, 1 = one competitor extinct) for each replicate (row) across all time steps (columns) from the GEM simulation. The "V" and "I" notation correspond to evolutionary and ecological scenarios as described for Abundances above. Note that there is no "W" index in the file name, as the extinction state is recorded at the level of the community (i.e., whether either competitor is extinct).\
\
(3) Traits_W_V_I.csv -\> These files contain the competitors' median trait values for each replicate (row) across all time steps (columns) from the GEM simulation. The "W" takes on values of "i" and "j", denoting whether the median trait value is calculated for the (arbitrarily chosen) focal population i or competitor population j. The "V" takes on values "j1", "j2", and "j3" corresponding to the three GEM variants, where "j1" = no trait variation, "j2" = trait variation but no heritability (referred to as the "no evolution" scenario in the manuscript), and "j3" = heritable trait variation (referred to as "evolution" scenario). The "I" takes on values "below", "at", and "above", denoting the position of the initial ecological state of the community relative to the unstable coexistence equilibrium.

#### R Code

(1) ConceptualFigures_20231024.R -\> This R script contains code to produce Figures 1 and 2 from the manuscript. These are conceptual figures derived from analysis of the deterministic ordinary differential equation for the birth-death 2-species competition model. This code can be executed in the absence of any other scripts or functions in the repository.\
\
(2) ExtinctionRates_20231024.R -\> This R script contains code for fitting time-to-event (survival) models to the extinction data generated by the GEM simulations (see extinction data files above). This procedure is performed for each ecological scenario separately (i.e., simulations beginning with competitors below, at, or above the unstable ecological equilibrium). If you wish to reproduce the original analyses from the manuscript, be sure to run this code using relative file paths that point to the existing extinction data files in the data folder.\
\
(3) TraitsAtFinalTime_20231025.R -\> This R script contains code for plotting the final trait values of the 2 competitors against each other and coloring the points by their final non-zero abundances (Figure 5 in the manuscript). The organization and cleaning of the data is performed for each ecological scenario separately (i.e., simulations beginning with competitors below, at, or above the unstable ecological equilibrium), along with plots to check each. The results for each ecological scenario are then compiled in a single dataframe so they can be easily plotted together. If you wish to reproduce the original Figure 5 from the manuscript, be sure to run this code using relative file paths that point to the existing trait value data files in the data folder.
