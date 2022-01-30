# modelling-stellar-mass-evolution-of-galaxies

This module charts the stellar mass evolution of galaxies employing observational constraints & mass loss models. It starts
this evolutionary sequence at a given redshift and for seed galaxies of a given stellar mass. The mass evolution then
occurs as a result star formation periods at given redshift intervals. The code tracks the stellar mass loss of these populations for each epoch after their production, and their contribution to the overall stellar mass of the galaxy.

Below is the figure the module outputs for 6 galaxies, with seed masses ranging from 10^9.25 to 10^11 Solar masses. The star formation sequence for this case has been started at z = 2, with star formation happening at redshift steps of 0.05. The star formation rate (SFR) options in the module is defined as a function of the star formation rate main sequence relation by design, as this was the star formation history we wanted to compare our sample to in S. Deger, G. Rudnick et al. (2022, in prep.). For the figure below, the galaxies have SFR's that is always equal to the SFR of the star-forming main sequence, the functional form of which is defined in Whitaker et al. (2012). 


<br/><img src='sfr_mstar_evolution.png'>
