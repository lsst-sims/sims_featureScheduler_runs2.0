#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
import healpy as hp
from rubin_sim.scheduler.modelObservatory import Model_observatory
from rubin_sim.scheduler.schedulers import Core_scheduler, simple_filter_sched
from rubin_sim.scheduler.utils import (Sky_area_generator, Footprint, Footprints,
                                       make_rolling_footprints, Step_slopes, ra_dec_hp_map)
import rubin_sim.scheduler.basis_functions as bf
from rubin_sim.scheduler.surveys import (Greedy_survey, generate_dd_surveys,
                                         Blob_survey)
from rubin_sim.scheduler import sim_runner
import rubin_sim.scheduler.detailers as detailers
import sys
import subprocess
import os
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
from rubin_sim.scheduler.basis_functions import Base_basis_function


# Putting this here for now so we can keep running everything on the same git hash.
# If we like it, migrate to rubin_sim
class Area_timed_mask_basis_function(Base_basis_function):
    """Mask an area, with possible time restrictions

    Parameters
    ----------
    input_mask : np.array
        A healpix map with values np.nan where things should be masked. 
    """
    def __init__(self, input_mask, night_min=-1, night_max=np.inf, side=32):
        super().__init__(nside=nside)
        self.input_mask = input_mask
        self.night_min = night_min
        self.night_max = night_max
        # Should check that input mask is proper nside
        self.nomask = np.ones(input_mask.size)

    def _calc_value(self, conditions, indx=None):
        
        if (conditions.night > self.night_min) & (conditions.night < self.night_max):
            result = self.input_mask
        else:
            result = self.nomask + 0
        return result 


class Timed_feasability_basis_function(Base_basis_function):
    """set a survey to not execute for certiain times
    """
    def __init__(self, min_nights_allowed=[-1, -1], max_nights_allowed=[np.inf, np.inf]):
        super().__init__()
        self.min_nights_allowed = min_nights_allowed
        self.max_nights_allowed = max_nights_allowed

    def check_feasibility(self, conditions):
        for n1, n2 in zip(self.min_nights_allowed, self.max_nights_allowed):
            if (conditions.night > n1) & (conditions.night < n2):
                return True
        
        return False


def make_rolling_footprints_plain(fp_hp=None, mjd_start=60218., sun_RA_start=3.27717639,
                                  nslice=2, scale=0.8, nside=32, wfd_indx=None):
    """
    Generate rolling footprints

    Parameters
    ----------
    fp_hp : dict-like
        A dict with filtername keys and HEALpix map values
    mjd_start : float
        The starting date of the survey.
    sun_RA_start : float
        The RA of the sun at the start of the survey
    nslice : int (2)
        How much to slice the sky up. Can be 2 or 3. Value of 6 to be implemented.
    scale : float (0.8)
        The strength of the rolling, value of 1 is full power rolling, zero is no rolling.
    wfd_indx : array of ints (none)
        The indices of the HEALpix map that are to be included in the rolling.

    Returns
    -------
    Footprints object
    """

    hp_footprints = fp_hp
   
    down = 1.-scale
    up = nslice - down*(nslice-1)
    start = [1., 1., 1.]
    end = [1., 1., 1., 1., 1., 1.]
    if nslice == 2:
        rolling = [up, down, up, down, up, down]
    elif nslice == 3:
        rolling = [up, down, down, up, down, down]
    elif nslice == 6:
        rolling = [up, down, down, down, down, down]
    all_slopes = [start + np.roll(rolling, i).tolist()+end for i in range(nslice)]

    fp_non_wfd = Footprint(mjd_start, sun_RA_start=sun_RA_start)
    rolling_footprints = []
    for i in range(nslice):
        step_func = Step_slopes(rise=all_slopes[i])
        rolling_footprints.append(Footprint(mjd_start, sun_RA_start=sun_RA_start,
                                            step_func=step_func))

    wfd = hp_footprints['r'] * 0
    if wfd_indx is None:
        wfd_indx = np.where(hp_footprints['r'] == 1)[0]
        non_wfd_indx = np.where(hp_footprints['r'] != 1)[0]

    wfd[wfd_indx] = 1
    non_wfd_indx = np.where(wfd == 0)[0] 

    split_wfd_indices = slice_galactic_cut(hp_footprints, nslice=nslice,
                                       wfd_indx=wfd_indx)

    for key in hp_footprints:
        temp = hp_footprints[key] + 0
        temp[wfd_indx] = 0
        fp_non_wfd.set_footprint(key, temp)

        for i in range(nslice):
            # make a copy of the current filter
            temp = hp_footprints[key] + 0
            # Set the non-rolling area to zero
            temp[non_wfd_indx] = 0
            
            indx = split_wfd_indices[i]
            # invert the indices
            ze = temp * 0
            ze[indx] = 1
            temp = temp * ze
            rolling_footprints[i].set_footprint(key, temp)

    result = Footprints([fp_non_wfd] + rolling_footprints)
    return result


def slice_galactic_cut(target_map, nslice=2, wfd_indx=None):
    """
    Helper function for generating rolling footprints

    Parameters
    ----------
    target_map : dict of HEALpix maps
        The final desired footprint as HEALpix maps. Keys are filter names
    nslice : int (2)
        The number of slices to make, can be 2 or 3.
    wfd_indx : array of ints
        The indices of target_map that should be used for rolling. If None, assumes
        the rolling area should be where target_map['r'] == 1.
    """

    ra, dec = ra_dec_hp_map(nside=hp.npix2nside(target_map['r'].size))

    coord = SkyCoord(ra=ra*u.rad, dec=dec*u.rad)
    gal_lon, gal_lat = coord.galactic.l.deg, coord.galactic.b.deg

    indx_north = np.intersect1d(np.where(gal_lat >= 0)[0], wfd_indx)
    indx_south = np.intersect1d(np.where(gal_lat < 0)[0], wfd_indx)
    
    splits_north = slice_wfd_indx(target_map, nslice=nslice, wfd_indx=indx_north)
    splits_south = slice_wfd_indx(target_map, nslice=nslice, wfd_indx=indx_south)

    slice_indx = []
    for j in np.arange(nslice):
        indx_temp = []
        for i in np.arange(j+1, nslice+1, nslice):
            indx_temp += indx_north[splits_north[i-1]:splits_north[i]].tolist()
            indx_temp += indx_south[splits_south[i-1]:splits_south[i]].tolist()
        slice_indx.append(indx_temp)

    return slice_indx


def slice_wfd_indx(target_map, nslice=2, wfd_indx=None):
    """
    simple map split
    """

    wfd = target_map['r'] * 0
    if wfd_indx is None:
        wfd_indx = np.where(target_map['r'] == 1)[0]
    wfd[wfd_indx] = 1
    wfd_accum = np.cumsum(wfd)
    split_wfd_indices = np.floor(np.max(wfd_accum)/nslice*(np.arange(nslice)+1)).astype(int)
    split_wfd_indices = split_wfd_indices.tolist()
    split_wfd_indices = [0] + split_wfd_indices

    return split_wfd_indices



def gen_greedy_surveys(nside=32, nexp=2, exptime=30., filters=['r', 'i', 'z', 'y'],
                       camera_rot_limits=[-80., 80.],
                       shadow_minutes=60., max_alt=76., moon_distance=30., ignore_obs='DD',
                       m5_weight=3., footprint_weight=0.75, slewtime_weight=3.,
                       stayfilter_weight=3., footprints=None,
                       area_mask=None, night_min=None, night_max=None):
    """
    Make a quick set of greedy surveys

    This is a convienence function to generate a list of survey objects that can be used with
    rubin_sim.scheduler.schedulers.Core_scheduler.
    To ensure we are robust against changes in the sims_featureScheduler codebase, all kwargs are
    explicitly set.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filters : list of str (['r', 'i', 'z', 'y'])
        Which filters to generate surveys for.
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    """
    # Define the extra parameters that are used in the greedy survey. I
    # think these are fairly set, so no need to promote to utility func kwargs
    greed_survey_params = {'block_size': 1, 'smoothing_kernel': None,
                           'seed': 42, 'camera': 'LSST', 'dither': True,
                           'survey_name': 'greedy'}

    surveys = []
    detailer = detailers.Camera_rot_detailer(min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits))

    for filtername in filters:
        bfs = []
        bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight))
        bfs.append((bf.Footprint_basis_function(filtername=filtername,
                                                footprint=footprints,
                                                out_of_bounds_val=np.nan, nside=nside), footprint_weight))
        bfs.append((bf.Slewtime_basis_function(filtername=filtername, nside=nside), slewtime_weight))
        bfs.append((bf.Strict_filter_basis_function(filtername=filtername), stayfilter_weight))
        # Masks, give these 0 weight
        bfs.append((bf.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=shadow_minutes,
                                                         max_alt=max_alt), 0))
        bfs.append((bf.Moon_avoidance_basis_function(nside=nside, moon_distance=moon_distance), 0))

        bfs.append((bf.Filter_loaded_basis_function(filternames=filtername), 0))
        bfs.append((bf.Planet_mask_basis_function(nside=nside), 0))
        bfs.append((Area_timed_mask_basis_function(area_mask,
                                                   night_min=night_min,
                                                   night_max=night_max), 0))

        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        surveys.append(Greedy_survey(basis_functions, weights, exptime=exptime, filtername=filtername,
                                     nside=nside, ignore_obs=ignore_obs, nexp=nexp,
                                     detailers=[detailer], **greed_survey_params))

    return surveys


def generate_blobs(nside, nexp=2, exptime=30., filter1s=['u', 'u', 'g', 'r', 'i', 'z', 'y'],
                   filter2s=['g', 'r', 'r', 'i', 'z', 'y', 'y'], pair_time=33.,
                   camera_rot_limits=[-80., 80.], n_obs_template=3,
                   season=300., season_start_hour=-4., season_end_hour=2.,
                   shadow_minutes=60., max_alt=76., moon_distance=30., ignore_obs='DD',
                   m5_weight=6., footprint_weight=1.5, slewtime_weight=3.,
                   stayfilter_weight=3., template_weight=12., footprints=None, u_nexp1=True,
                   area_mask=None, night_min=None, night_max=None,
                   timed_feas_min=None, timed_feas_max=None):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (33)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : int (3)
        The number of observations to take every season in each filter
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_nexp1 : bool (True)
        Add a detailer to make sure the number of expossures in a visit is always 1 for u observations.
    """

    blob_survey_params = {'slew_approx': 7.5, 'filter_change_approx': 140.,
                          'read_approx': 2., 'min_pair_time': 15., 'search_radius': 30.,
                          'alt_max': 85., 'az_range': 90., 'flush_time': 30.,
                          'smoothing_kernel': None, 'nside': nside, 'seed': 42, 'dither': True,
                          'twilight_scale': True}

    surveys = []

    times_needed = [pair_time, pair_time*2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(detailers.Camera_rot_detailer(min_rot=np.min(camera_rot_limits),
                                                           max_rot=np.max(camera_rot_limits)))
        detailer_list.append(detailers.Close_alt_detailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        if filtername2 is not None:
            bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight/2.))
            bfs.append((bf.M5_diff_basis_function(filtername=filtername2, nside=nside), m5_weight/2.))

        else:
            bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight))

        if filtername2 is not None:
            bfs.append((bf.Footprint_basis_function(filtername=filtername,
                                                    footprint=footprints,
                                                    out_of_bounds_val=np.nan, nside=nside), footprint_weight/2.))
            bfs.append((bf.Footprint_basis_function(filtername=filtername2,
                                                    footprint=footprints,
                                                    out_of_bounds_val=np.nan, nside=nside), footprint_weight/2.))
        else:
            bfs.append((bf.Footprint_basis_function(filtername=filtername,
                                                    footprint=footprints,
                                                    out_of_bounds_val=np.nan, nside=nside), footprint_weight))

        bfs.append((bf.Slewtime_basis_function(filtername=filtername, nside=nside), slewtime_weight))
        bfs.append((bf.Strict_filter_basis_function(filtername=filtername), stayfilter_weight))

        if filtername2 is not None:
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername, nside=nside,
                                                         footprint=footprints.get_footprint(filtername),
                                                         n_obs=n_obs_template, season=season,
                                                         season_start_hour=season_start_hour,
                                                         season_end_hour=season_end_hour), template_weight/2.))
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername2, nside=nside,
                                                         footprint=footprints.get_footprint(filtername2),
                                                         n_obs=n_obs_template, season=season,
                                                         season_start_hour=season_start_hour,
                                                         season_end_hour=season_end_hour), template_weight/2.))
        else:
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername, nside=nside,
                                                         footprint=footprints.get_footprint(filtername),
                                                         n_obs=n_obs_template, season=season,
                                                         season_start_hour=season_start_hour,
                                                         season_end_hour=season_end_hour), template_weight))
        # Masks, give these 0 weight
        bfs.append((bf.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt,
                                                         penalty=np.nan, site='LSST'), 0.))
        bfs.append((bf.Moon_avoidance_basis_function(nside=nside, moon_distance=moon_distance), 0.))
        filternames = [fn for fn in [filtername, filtername2] if fn is not None]
        bfs.append((bf.Filter_loaded_basis_function(filternames=filternames), 0))
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.Time_to_twilight_basis_function(time_needed=time_needed), 0.))
        bfs.append((bf.Not_twilight_basis_function(), 0.))
        bfs.append((bf.Planet_mask_basis_function(nside=nside), 0.))
        bfs.append((Area_timed_mask_basis_function(area_mask,
                                                   night_min=night_min,
                                                   night_max=night_max), 0))
        if timed_feas_min is not None:
            bfs.append((Timed_feasability_basis_function(min_nights_allowed=timed_feas_min,
                                                         max_nights_allowed=timed_feas_max), 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = 'blob, %s' % filtername
        else:
            survey_name = 'blob, %s%s' % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.Take_as_pairs_detailer(filtername=filtername2))

        if u_nexp1:
            detailer_list.append(detailers.Filter_nexp(filtername='u', nexp=1))
        surveys.append(Blob_survey(basis_functions, weights, filtername1=filtername, filtername2=filtername2,
                                   exptime=exptime,
                                   ideal_pair_time=pair_time,
                                   survey_note=survey_name, ignore_obs=ignore_obs,
                                   nexp=nexp, detailers=detailer_list, **blob_survey_params))

    return surveys


def generate_twi_blobs(nside, nexp=2, exptime=30., filter1s=['r', 'i', 'z', 'y'],
                       filter2s=['i', 'z', 'y', 'y'], pair_time=15.,
                       camera_rot_limits=[-80., 80.], n_obs_template=3,
                       season=300., season_start_hour=-4., season_end_hour=2.,
                       shadow_minutes=60., max_alt=76., moon_distance=30., ignore_obs='DD',
                       m5_weight=6., footprint_weight=1.5, slewtime_weight=3.,
                       stayfilter_weight=3., template_weight=12., footprints=None, repeat_night_weight=None,
                       wfd_footprint=None,
                       area_mask=None, night_min=None, night_max=None):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (22)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : int (3)
        The number of observations to take every season in each filter
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    """

    blob_survey_params = {'slew_approx': 7.5, 'filter_change_approx': 140.,
                          'read_approx': 2., 'min_pair_time': 10., 'search_radius': 30.,
                          'alt_max': 85., 'az_range': 90., 'flush_time': 30.,
                          'smoothing_kernel': None, 'nside': nside, 'seed': 42, 'dither': True,
                          'twilight_scale': False, 'in_twilight': True}

    surveys = []

    times_needed = [pair_time, pair_time*2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(detailers.Camera_rot_detailer(min_rot=np.min(camera_rot_limits),
                                                           max_rot=np.max(camera_rot_limits)))
        detailer_list.append(detailers.Close_alt_detailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        if filtername2 is not None:
            bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight/2.))
            bfs.append((bf.M5_diff_basis_function(filtername=filtername2, nside=nside), m5_weight/2.))

        else:
            bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight))

        if filtername2 is not None:
            bfs.append((bf.Footprint_basis_function(filtername=filtername,
                                                    footprint=footprints,
                                                    out_of_bounds_val=np.nan, nside=nside), footprint_weight/2.))
            bfs.append((bf.Footprint_basis_function(filtername=filtername2,
                                                    footprint=footprints,
                                                    out_of_bounds_val=np.nan, nside=nside), footprint_weight/2.))
        else:
            bfs.append((bf.Footprint_basis_function(filtername=filtername,
                                                    footprint=footprints,
                                                    out_of_bounds_val=np.nan, nside=nside), footprint_weight))

        bfs.append((bf.Slewtime_basis_function(filtername=filtername, nside=nside), slewtime_weight))
        bfs.append((bf.Strict_filter_basis_function(filtername=filtername), stayfilter_weight))

        if filtername2 is not None:
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername, nside=nside,
                                                         footprint=footprints.get_footprint(filtername),
                                                         n_obs=n_obs_template, season=season,
                                                         season_start_hour=season_start_hour,
                                                         season_end_hour=season_end_hour), template_weight/2.))
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername2, nside=nside,
                                                         footprint=footprints.get_footprint(filtername2),
                                                         n_obs=n_obs_template, season=season,
                                                         season_start_hour=season_start_hour,
                                                         season_end_hour=season_end_hour), template_weight/2.))
        else:
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername, nside=nside,
                                                         footprint=footprints.get_footprint(filtername),
                                                         n_obs=n_obs_template, season=season,
                                                         season_start_hour=season_start_hour,
                                                         season_end_hour=season_end_hour), template_weight))
        if repeat_night_weight is not None:
            bfs.append((bf.Avoid_long_gaps_basis_function(nside=nside, filtername=None,
                                                          min_gap=0., max_gap=10./24., ha_limit=3.5,
                                                          footprint=wfd_footprint), repeat_night_weight))
        # Masks, give these 0 weight
        bfs.append((bf.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt,
                                                         penalty=np.nan, site='LSST'), 0.))
        bfs.append((bf.Moon_avoidance_basis_function(nside=nside, moon_distance=moon_distance), 0.))
        filternames = [fn for fn in [filtername, filtername2] if fn is not None]
        bfs.append((bf.Filter_loaded_basis_function(filternames=filternames), 0))
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.Time_to_twilight_basis_function(time_needed=time_needed, alt_limit=12), 0.))
        bfs.append((bf.Planet_mask_basis_function(nside=nside), 0.))
        bfs.append((Area_timed_mask_basis_function(area_mask,
                                                   night_min=night_min,
                                                   night_max=night_max), 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = 'blob_twi, %s' % filtername
        else:
            survey_name = 'blob_twi, %s%s' % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.Take_as_pairs_detailer(filtername=filtername2))
        surveys.append(Blob_survey(basis_functions, weights, filtername1=filtername, filtername2=filtername2,
                                   exptime=exptime,
                                   ideal_pair_time=pair_time,
                                   survey_note=survey_name, ignore_obs=ignore_obs,
                                   nexp=nexp, detailers=detailer_list, **blob_survey_params))

    return surveys


def run_sched(surveys, survey_length=365.25, nside=32, fileroot='baseline_', verbose=False,
              extra_info=None, illum_limit=40.):
    years = np.round(survey_length/365.25)
    scheduler = Core_scheduler(surveys, nside=nside)
    n_visit_limit = None
    filter_sched = simple_filter_sched(illum_limit=illum_limit)
    observatory = Model_observatory(nside=nside)
    observatory, scheduler, observations = sim_runner(observatory, scheduler,
                                                      survey_length=survey_length,
                                                      filename=fileroot+'%iyrs.db' % years,
                                                      delete_past=True, n_visit_limit=n_visit_limit,
                                                      verbose=verbose, extra_info=extra_info,
                                                      filter_scheduler=filter_sched)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", dest='verbose', action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument("--survey_length", type=float, default=365.25*10)
    parser.add_argument("--outDir", type=str, default="")
    parser.add_argument("--maxDither", type=float, default=0.7, help="Dither size for DDFs (deg)")
    parser.add_argument("--moon_illum_limit", type=float, default=40., help="illumination limit to remove u-band")
    parser.add_argument("--nexp", type=int, default=2)
    parser.add_argument("--rolling_nslice", type=int, default=2)
    parser.add_argument("--rolling_strength", type=float, default=0.9)
    parser.add_argument("--dbroot", type=str)
    parser.add_argument('--filters', help="filter distribution (default: u 0.07 g 0.09 r 0.22 i 0.22 z 0.20 y 0.20)")
    parser.add_argument("--same_pairs", action="store_true", default=False)

    args = parser.parse_args()
    survey_length = args.survey_length  # Days
    outDir = args.outDir
    verbose = args.verbose
    max_dither = args.maxDither
    illum_limit = args.moon_illum_limit
    nexp = args.nexp
    nslice = args.rolling_nslice
    scale = args.rolling_strength
    filters = args.filters
    dbroot = args.dbroot

    nside = 32
    per_night = True  # Dither DDF per night

    camera_ddf_rot_limit = 75.

    extra_info = {}
    exec_command = ''
    for arg in sys.argv:
        exec_command += ' ' + arg
    extra_info['exec command'] = exec_command
    try:
        extra_info['git hash'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
    except subprocess.CalledProcessError:
        extra_info['git hash'] = 'Not in git repo'

    extra_info['file executed'] = os.path.realpath(__file__)

    # Use the filename of the script to name the output database
    if dbroot is None:
        fileroot = os.path.basename(sys.argv[0]).replace('.py', '') + '_'
    else:
        fileroot = dbroot + '_'
    file_end = 'v2.0_'

    if filters is None:
        sm = Sky_area_generator(nside=nside)
    else:
        filter_split = filters.split(" ")
        filter_balance = {x: float(y) for x, y in zip(filter_split[::2], filter_split[1::2])}
        if not np.isclose(1.0, np.sum([v for v in filter_balance.values()])):
            raise ValueError('Make sure your filters sum to 1')
        sm = Sky_area_generator(nside=nside,
                                default_filter_balance=filter_balance)

    sm.set_maps()
    final_tot, footprints_hp = sm.return_maps()
    # Set the wfd, aka rolling, pixels
    wfd_footprint = sm.maps['dustfree']
    wfd_indx = np.where(wfd_footprint > 0)[0]

    normval = footprints_hp['r'][wfd_indx].min()
    for key in footprints_hp:
        footprints_hp[key] = footprints_hp[key]/normval

    # XXX--hacky hack. The sky_area_generator should be saving the actual bulge pixles for us.
    bulge_indx = np.where(footprints_hp['r'] > 1)[0]
    non_bulge_indx = np.where(footprints_hp['r'] <= 1)[0]

    repeat_night_weight = None

    observatory = Model_observatory(nside=nside)
    conditions = observatory.return_conditions()

    footprints = make_rolling_footprints(fp_hp=footprints_hp, mjd_start=conditions.mjd_start,
                                         sun_RA_start=conditions.sun_RA_start, nslice=nslice, scale=scale,
                                         nside=nside, wfd_indx=wfd_indx)

    bulge_footprints = make_rolling_footprints_plain(fp_hp=footprints_hp, mjd_start=conditions.mjd_start,
                                                     sun_RA_start=conditions.sun_RA_start, nslice=6, scale=scale,
                                                     nside=nside, wfd_indx=bulge_indx)

    #for fpl in bulge_footprints.footprint_list[1:]:
    #    for bfp, key in zip(fpl.footprints, footprints_hp):
    #        bfp[non_bulge_indx] = footprints_hp[key][non_bulge_indx] + 0

    # Because of the way the footprint basis function normalizes things, we 
    # need the regular footprints and the bulge footprints to have the same total values.
    # So we'll patch them up to match, then use a masking basis function to mask the bulge
    # in the regular footprints, and vice-versa in the bulge footprints. 

    # For the non-rolling years, we'll need to turn off the mask.

    # guess we need to do the same thing on the twilight blobs and greedy that we do to the main blobs.
    # this means we won't be doing bulge observations in twilight in rolling years, which is probably fine.

    mask_out_bulge = np.ones(hp.nside2npix(nside))
    mask_out_bulge[bulge_indx] = np.nan

    mask_out_wfd = np.zeros(hp.nside2npix(nside)) + np.nan
    mask_out_wfd[bulge_indx] = 1

    night_min = 365.25*1.5
    night_max = 365.25*7.5

    # Set up the DDF surveys to dither
    u_detailer = detailers.Filter_nexp(filtername='u', nexp=1)
    dither_detailer = detailers.Dither_detailer(per_night=per_night, max_dither=max_dither)
    details = [detailers.Camera_rot_detailer(min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit),
               dither_detailer, u_detailer]
    euclid_detailers = [detailers.Camera_rot_detailer(min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit),
                        detailers.Euclid_dither_detailer(), u_detailer]
    ddfs = generate_dd_surveys(nside=nside, nexp=nexp, detailers=details, euclid_detailers=euclid_detailers)

    greedy = gen_greedy_surveys(nside, nexp=nexp, footprints=footprints,
                                area_mask=mask_out_bulge,
                                night_min=night_min,
                                night_max=night_max)

    blobs = generate_blobs(nside, nexp=nexp, footprints=footprints,
                           area_mask=mask_out_bulge,
                           night_min=night_min,
                           night_max=night_max)

    bulge_blobs = generate_blobs(nside, nexp=nexp, footprints=bulge_footprints,
                                 filter1s=['u', 'g', 'r', 'i', 'z', 'y'],
                                 filter2s=[None, None, None, None, None, None], pair_time=5.,
                                 area_mask=mask_out_wfd,
                                 night_min=night_min,
                                 night_max=night_max,
                                 timed_feas_min=[night_min],
                                 timed_feas_max=[night_max])

    twi_blobs = generate_twi_blobs(nside, nexp=nexp,
                                   footprints=footprints,
                                   wfd_footprint=wfd_footprint,
                                   repeat_night_weight=repeat_night_weight,
                                   area_mask=mask_out_bulge,
                                   night_min=night_min,
                                   night_max=night_max)
   
    surveys = [ddfs, blobs+bulge_blobs, twi_blobs, greedy]
    run_sched(surveys, survey_length=survey_length, verbose=verbose,
              fileroot=os.path.join(outDir, fileroot+file_end), extra_info=extra_info,
              nside=nside, illum_limit=illum_limit)
