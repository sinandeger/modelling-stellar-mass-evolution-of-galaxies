"""
Author: Sinan Deger || Date: 15 January 2022

code origin: https://github.com/sinandeger/modelling-stellar-mass-evolution-of-galaxies.git
"""
import os
"""Base python libraries"""
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
"""Extra task-dependent libraries"""
from matplotlib.pyplot import cm
from astropy.cosmology import FlatLambdaCDM


"""This module charts the redshift evolution of of given """

"""Define the governing cosmology background here"""
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def sfr_mstar_relation(z, mstar):

    """
    This function returns the star formation rate a galaxy of mass mstar would have on the star-formation main sequence at redshift z.
    Adopted from Whitaker et al. (2012) (https://ui.adsabs.harvard.edu/abs/2012ApJ...754L..29W/abstract). Returns the SFR per year.

    Note that the z needs to be less than or equal to 2.5.

    :type z: float
    :type mstar: float, in solar masses
    """

    alpha_z = 0.70 - 0.13*z  # eq. (2) in Whitaker et al. (2012)
    beta_z = 0.38 + 1.14*z - 0.19*np.power(z, 2)   # eq. (3)

    log_sfr = alpha_z*(np.log10(mstar) - 10.5) + beta_z   # eq. (1)
    """Convert logarithm of SFR into SFR per year and return"""
    return np.power(10, log_sfr)


def time_interval_yr(higher_z, lower_z):

    """
    This function returns the time range in years between the two redshift values, lower_z and higher_z

    :type lower_z: float
    :type higher_z: float
    """

    age_upper_str = re.findall("[+-]?\d+\.\d+", str(cosmo.age(higher_z)))
    age_upper_nounit = float(next(iter(age_upper_str)))

    age_lower_str = re.findall("[+-]?\d+\.\d+", str(cosmo.age(lower_z)))
    age_lower_nounit = float(next(iter(age_lower_str)))

    time_elapsed = age_lower_nounit - age_upper_nounit

    return time_elapsed*np.power(10, 9)


def bc03_mstar_loss(redshift_upper_limit, redshift_lower_limit, bc03_df):

    """
    This function returns the expected stellar mass loss of a galaxy based on the evolution of BC03 simple stellar population (SSP) models

    This tool has been replaced by the mass loss mechanism described in poggianti_mass_loss in the current installation.

    """

    """First, turn the redshift range into an age range"""

    age_upper_str = re.findall("[+-]?\d+\.\d+", str(cosmo.age(redshift_upper_limit)))
    age_upper_nounit = float(next(iter(age_upper_str)))

    age_lower_str = re.findall("[+-]?\d+\.\d+", str(cosmo.age(redshift_lower_limit)))
    age_lower_nounit = float(next(iter(age_lower_str)))

    time_elapsed = age_lower_nounit - age_upper_nounit
    log_time_elapsed = 9.0 + np.log10(time_elapsed)
    print('Time elapsed in Gyr: ', round(time_elapsed, 5))
    print('log time elapsed', log_time_elapsed)

    """Take log(age) = 8 as the start of the mass loss cycle"""

    sfr_cycle_start = 8.0
    max_fmstar = bc03_df.loc[bc03_df['log-age-yr'] == sfr_cycle_start]['M*_liv'].item()/np.max(bc03_df['M*_liv'].values)

    nearest_bc03_log_age = find_nearest(bc03_df['log-age-yr'].values, log_time_elapsed)
    print('Nearest BC03 log age:', nearest_bc03_log_age)
    nearest_fmstar = bc03_df.loc[bc03_df['log-age-yr'] == nearest_bc03_log_age]['M*_liv'].item()/np.max(bc03_df['M*_liv'].values)

    fmstar_lost_in_time_elapsed = max_fmstar - nearest_fmstar
    print('Fraction of stellar mass lost in time elapsed: ', round(fmstar_lost_in_time_elapsed, 3))

    return fmstar_lost_in_time_elapsed


def find_nearest(array, value):

    """This function finds the array member nearest in absolute value to a given query value"""

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def poggianti_mass_loss(time_elapsed):

    """
    This function implements the mass loss approximation described in Poggianti et al. (2013)
    (https://iopscience.iop.org/article/10.1088/0004-637X/777/2/125/pdf)

    It needs the time elapsed in years for the interval the mass loss occurs.
    It then returns the fraction of the mass that remains from the beginning of the time interval.

    The underlying stellar population is assumed to have Solar metallicity, a Chabrier IMF, and obey the BC03 prescriptions.
    More detailed information can be found in Section 2.1, Poggianti et al. (2013).

    """

    if time_elapsed < 1.9 * np.power(10, 6):
        remaining_mass_fraction = 1.0  # Mass loss is set to be none in the first million years
    else:
        remaining_mass_fraction = 1.749 - 0.124 * np.log10(time_elapsed)

    return remaining_mass_fraction


def evolve_stellar_population(evolve_z, fiducial_mstar, fiducial_z, z_step, below_ms_sfr=False):

    """This function evolves a set of stellar populations, tracking their mass evolution and
     computing the mass at the formation of each new population"""

    """Note: below_ms_sfr sets the SFR to 0.3 dex below the main sequence by default"""

    """First, split the range between the fiducial_z and evolve_z into bins of size z_step, in descending order"""
    sf_epochs = np.flip(np.arange(start=evolve_z, stop=fiducial_z, step=z_step))
    """
    population_dict carries information about the stellar populations at the time of their formation.
    
    galaxy_dict carries information about the overall galaxy properties, such as stellar mass at redshift intervals.
    
    """
    population_dict = {}
    pop_evolve_df = pd.DataFrame()
    galaxy_dict = {}
    galaxy_dict.update({'mstar_post_sf': [fiducial_mstar]})
    if below_ms_sfr:
        galaxy_dict.update({'sfr_at_z': [sfr_mstar_relation(fiducial_z, fiducial_mstar)/2.0]})
    else:
        galaxy_dict.update({'sfr_at_z': [sfr_mstar_relation(fiducial_z, fiducial_mstar)]})

    for z_ind, z_sf in enumerate(sf_epochs):
        contributing_epochs = sf_epochs[:z_ind]
        pop_evolve_df.loc[z_ind, 'redshift'] = z_sf

        if len(contributing_epochs) == 0:
            mstar_at_epoch_z = fiducial_mstar
            """Take the average of SFR at the beginning and the end of the star formation epoch for sfr_at_z"""
            if below_ms_sfr:
                """Currently, the module only allows for SFR 0.3 dex below main sequence for the below_ms_sfr == True case.
                 This will be modified to allow for a user-defined SFR soon."""
                sfr_at_z = (sfr_mstar_relation(z_sf, mstar_at_epoch_z) + sfr_mstar_relation(z_sf - z_step,
                                                                                            mstar_at_epoch_z)) / 4.0
            else:
                sfr_at_z = (sfr_mstar_relation(z_sf, mstar_at_epoch_z) + sfr_mstar_relation(z_sf - z_step, mstar_at_epoch_z)) / 2.0

            """SFR duration is the time elapsed in years within the redshift step"""
            sfr_duration = time_interval_yr(z_sf, z_sf - z_step)
            mstar_z_sf = sfr_at_z * sfr_duration
            population_dict.update({'population_' + str(z_ind): [mstar_z_sf]})
            pop_evolve_df.loc[z_ind, 'time_spent_yr'] = sfr_duration
            pop_evolve_df.loc[z_ind, 'mstar_post_sf'] = mstar_at_epoch_z + mstar_z_sf
            pop_evolve_df.loc[z_ind, 'population_' + str(z_ind)] = mstar_z_sf
            galaxy_dict['mstar_post_sf'].append(mstar_at_epoch_z + mstar_z_sf)
            galaxy_dict['sfr_at_z'].append(sfr_at_z)

        mass_contrib = 0.0
        for contrib_ind, contrib_z in enumerate(contributing_epochs):
            individual_stellar_pop_contrib = population_dict['population_'+str(contrib_ind)][0] * poggianti_mass_loss(time_interval_yr(contrib_z, z_sf))
            mass_contrib += individual_stellar_pop_contrib
            pop_evolve_df.loc[z_ind, 'population_' + str(contrib_ind)] = individual_stellar_pop_contrib

        mstar_at_epoch_z = fiducial_mstar + mass_contrib
        """Take the average of SFR at the beginning and the end of the star formation epoch for sfr_at_z"""
        if below_ms_sfr:
            sfr_at_z = (sfr_mstar_relation(z_sf, mstar_at_epoch_z) + sfr_mstar_relation(z_sf - z_step,
                                                                                        mstar_at_epoch_z)) / 4.0
            sfr_indicator_tag = 'below_ms_sfr'

        else:
            sfr_at_z = (sfr_mstar_relation(z_sf, mstar_at_epoch_z) + sfr_mstar_relation(z_sf - z_step,
                                                                                        mstar_at_epoch_z)) / 2.0
            sfr_indicator_tag = 'ms_sfr'

        """SFR duration is the time elapsed in the redshift step"""
        sfr_duration = time_interval_yr(z_sf, z_sf - z_step)
        mstar_z_sf = sfr_at_z * sfr_duration
        galaxy_dict['mstar_post_sf'].append(mstar_at_epoch_z + mstar_z_sf)
        galaxy_dict['sfr_at_z'].append(sfr_at_z)

        population_dict.update({'population_' + str(z_ind): [mstar_z_sf]})
        pop_evolve_df.loc[z_ind, 'population_' + str(z_ind)] = mstar_z_sf
        pop_evolve_df.loc[z_ind, 'time_spent_yr'] = sfr_duration
        pop_evolve_df.loc[z_ind, 'mstar_post_sf'] = mstar_at_epoch_z + mstar_z_sf

    """Save information specific to the stellar populations and the galaxy as csv files"""
    population_df = pd.DataFrame.from_dict(population_dict)
    population_df.to_csv('_fid_logmass_'+str(np.log10(fiducial_mstar))+'_evolve_z_'+str(evolve_z)+'_population_dict.csv')

    galaxy_df = pd.DataFrame.from_dict(galaxy_dict)
    galaxy_df.to_csv('fid_logmass_'+str(np.log10(fiducial_mstar))+'_evolve_z_'+str(evolve_z)+'_galaxy_dict.csv')

    pop_evolve_df.to_csv('fid_logmass_'+str(np.log10(fiducial_mstar))+'_evolve_z_'+str(evolve_z)+'_population_evolve.csv')


def evolve_call(fiducial_masses, evolve_redshift, fiducial_z, z_step, below_ms_sfr=False, plot_sfr_mstar_sequence=False):

    """
    This function is the main call that triggers the evolution modelling sequence. This module was first designed to
    perform the evolution tracking the star formation main sequence only, but limited capability to allow for SFR's below
    that of the main sequence value has been added. Currently, only an SFR 0.3 dex below the main sequence value is allowed.

    Please note that currently no internal quenching mechanism is implemented, so seed galaxies would have ongoing star formation until
    the current evolve redshift.


    :param fiducial_masses: the seed galaxy stellar masses at redshift = fiducial_z  , type: list, array-like
    :param evolve_redshift: the redshift the seed galaxies will be evolved up to  , type: list, array-like
    :param fiducial_z: beginning of the evolution sequence, type: float
    :param z_step: the redshift step size at which a new generation of stars will be formed, type: float
    :param below_ms_sfr: whether the galaxy should have an SFR of 0.3 dex blow the star formation main sequence (SF MS),
     or at the SF MS.  type: boolean
    :param plot_sfr_mstar_sequence: whether a figure showing the evolution of seed masses on the SFR-M_star space should
     be plotted or not.   type: boolean

    """

    for fid_ind, fid_mass in enumerate(fiducial_masses):
        for fid_z_ind, fid_evolve_z in enumerate(evolve_redshift):
            evolve_stellar_population(fid_evolve_z, np.power(10, fid_mass), fiducial_z, z_step,
                                      below_ms_sfr=below_ms_sfr)

    """Plot the stellar mass evolution of the seed masses on the star formation rate versus stellar mass figure below"""
    if plot_sfr_mstar_sequence:

        """Define a color sequence the length of models to be evolved, to be drawn from a matplotlib colormap"""
        evolve_model_c = iter(cm.rainbow(np.linspace(0, 1, len(fiducial_masses))))

        for fid_mstar_ind, fid_mstar_ in enumerate(fiducial_masses):
            model_color = next(evolve_model_c)
            data_df = pd.read_csv('fid_logmass_' + str(fid_mstar_) + '_evolve_z_' + str(evolve_redshift[0]) + '_galaxy_dict.csv')

            plt.scatter(np.log10(data_df['mstar_post_sf']), np.log10(data_df['sfr_at_z']), marker='X', s=10,
                        color=model_color,
                        label='Fid. log(M_*): 10^{fid_m}'.format(fid_m=fid_mstar_))

        """To plot the main sequence relation at each evolve_redshift, we first find the MS SFR at each seed mass"""
        ms_fig_dict = {}
        """To have the figure extend past the largest value of fiducial masses, define a stellar mass range array"""
        mod_fiducial_masses = np.arange(np.min(fiducial_masses), np.max(fiducial_masses)+1.0, step=0.5)

        for z_ms_ind_, z_ms_ in enumerate(evolve_redshift):
            for mstar_ms_ind_, mstar_ms_ in enumerate(mod_fiducial_masses):
                if mstar_ms_ind_ == 0:
                    ms_fig_dict.update({'sfr_at_z'+str(z_ms_): [sfr_mstar_relation(z_ms_, np.power(10, mstar_ms_))]})
                    ms_fig_dict.update({'mstar_at_z' + str(z_ms_): [mstar_ms_]})
                else:
                    ms_fig_dict['sfr_at_z'+str(z_ms_)].append(sfr_mstar_relation(z_ms_, np.power(10, mstar_ms_)))
                    ms_fig_dict['mstar_at_z' + str(z_ms_)].append(mstar_ms_)

        """Convert ms_fig_dict to a pandas dataframe for ease of access"""
        ms_fig_df = pd.DataFrame.from_dict(ms_fig_dict)

        mseq_c = iter(cm.cividis(np.linspace(0, 1, len(evolve_redshift))))

        for z_ind, z_sfr_rel in enumerate(evolve_redshift):
            mseq_color = next(mseq_c)
            plt.plot(ms_fig_df['mstar_at_z' + str(z_sfr_rel)], np.log10(ms_fig_df['sfr_at_z'+str(z_sfr_rel)].values),
                     label='SFR MS Rel. at z = ' + str(z_sfr_rel), color=mseq_color)

        plt.xlabel('log($\mathrm{M_{*}}$) [$\mathrm{M_{\odot}}$]')
        plt.ylabel('log(SFR) [$\mathrm{M_{\odot}}$/yr]')
        plt.title('Models evolved until z = ' + str(fid_evolve_z_list[0]))
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        # plt.xlim(9.0, 12.2)
        plt.legend(by_label.values(), by_label.keys(), prop={'size': 6})
        plt.savefig('sfr_mstar_evolution.pdf', format='pdf')
        plt.close()


"""An example run is provided here. First, we set the lists of fiducial seed masses,
 and the redshift limits the models to which the models will be evolved."""
fiducial_mass_list = [9.25, 9.5, 10.0, 10.5, 10.75, 11.10]
fid_evolve_z_list = [0.4, 0.6, 0.8]  #  increasing order in z
"""Next, run the call that starts the evolution sequence"""
evolve_call(fiducial_mass_list, fid_evolve_z_list, fiducial_z=2.0, z_step=0.05, below_ms_sfr=False, plot_sfr_mstar_sequence=True)






