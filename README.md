# Modelling the Stellar Mass Evolution of Galaxies

This module charts the stellar mass evolution of galaxies employing observational constraints & mass loss models. It starts
this evolutionary sequence at a given redshift and for seed galaxies of a given stellar mass. The mass evolution then
occurs as a result star formation periods at given redshift intervals, and subsequent mass loss applied to these stellar populations. The code tracks the stellar mass loss of these populations for each epoch after their production, and their contribution to the overall stellar mass of the galaxy.

Below is the figure the module outputs for 6 galaxies, with seed masses ranging from 10^9.25 to 10^11 Solar masses, as shown in the figure legend. The star formation sequence for this case has been started at z = 2, with star formation happening at redshift steps of 0.05. The star formation rate (SFR) options in the module is defined as a function of the star formation rate main sequence relation by design, as this was the star formation history we wanted to compare our sample to in S. Deger, G. Rudnick et al. (2022, in prep.). For the figure below, the galaxies have SFR's that is always equal to the SFR of the star-forming main sequence, the functional form of which is defined in Whitaker et al. (2012). The three lines in the figure indicate the main sequence at three different redshift values, as identified within the legend. The models are evolved until a redshift of 0.4.


<br/><img src='sfr_mstar_evolution.png'>


This module is not intended as an exhaustive modeling of the stellar mass evolution of galaxies, yet it serves as a versatile tool that produces a realistic approximation. Some features that are intended to be added soon include an intrinsic quenching mechanism that shuts of the star formation of seed galaxies based on observational constraints, and the addition of more flexible star formation history mechanics.
