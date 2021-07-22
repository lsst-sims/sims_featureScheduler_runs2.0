import numpy as np
import matplotlib.pylab as plt
import healpy as hp
from rubin_sim.scheduler.modelObservatory import Model_observatory
from rubin_sim.scheduler.schedulers import Core_scheduler, simple_filter_sched
from rubin_sim.scheduler.utils import combo_dust_fp, Footprint, Footprints, Step_slopes
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
import astropy.units as u
from rubin_sim import data as rs_data
from rubin_sim.utils import _angularSeparation

class SurveyMap:
    def __init__(self, nside=64, default_filter_balance=None):
        self.nside = nside
        # healpix indexes
        self.hpid = np.arange(0, hp.nside2npix(nside))
        # Ra/dec in degrees and other coordinates
        self.ra, self.dec = hp.pix2ang(nside, self.hpid, lonlat=True)
        self.coord = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg, frame='icrs')
        self.eclip_lat = self.coord.barycentrictrueecliptic.lat.deg
        self.eclip_lon = self.coord.barycentrictrueecliptic.lon.deg
        self.gal_lon = self.coord.galactic.l.deg
        self.gal_lat = self.coord.galactic.b.deg
        # filterlist
        self.filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
        # SRD values
        self.nvis_min_srd = 750
        self.nvis_goal_srd = 825
        self.area_min_srd = 15000
        self.area_goal_srd = 18000
        if default_filter_balance is None:
            self.default_filter_balance = {'u': 0.07, 'g': 0.09, 'r': 0.22, 
                                           'i': 0.22, 'z': 0.20, 'y': 0.20}
        else:
            self.default_filter_balance = self._normalize_filter_balance(default_filter_balance)
        self.maps = {}
        self.maps_perfilter = {}
        self.nvis = {}
    
    def read_dustmap(self, dustmapFile=None):
        # Dustmap from rubin_sim_data  - this is basically just a data directory
        # The dustmap data is downloadable from 
        # https://lsst.ncsa.illinois.edu/sim-data/rubin_sim_data/maps_may_2021.tgz
        # (then just set RUBIN_SIM_DATA_DIR to where you downloaded it, after untarring the file)
        if dustmapFile is None:
            datadir = rs_data.get_data_dir()
            if datadir is None:
                raise Exception('Cannot find datadir, please set "RUBIN_SIM_DATA_DIR"')
            datadir = os.path.join(datadir, 'maps', 'DustMaps')
            filename = os.path.join(datadir, 'dust_nside_%i.npz' % self.nside)
        self.dustmap = np.load(filename)['ebvMap']
        
    def _normalize_filter_balance(self, filter_balance):
        filtersum = np.array(list(filter_balance.values())).sum()
        tmp = {k: round(v/filtersum, 2) for k, v in filter_balance.items()}
        for f in self.filterlist:
            if f not in tmp:
                tmp[f] = 0
        return tmp
    
    # The various regions take the approach that they should be independent
    # And after setting all of the regions, we take the max value (per filter?) in each part of the sky
    # The individual components are updated, so that we can still calculate survey fraction per part of the sky
    
    def _set_exwfd(self, dust_limit=0.199, dec_min=-67, dec_max=12, 
                   smoothing_cutoff=0.57, smoothing_beam=10,
                   nvis_exwfd=825*1.08,
                   exgal_filter_balance=None):
        # Define extragalactic WFD between dec_min and dec_max with low dust extinction 
        # These dec and dust limits are used to define the other survey areas as well.
        # This means setting each part of the footprint depends on having already set some previous pieces..
        self.dust_limit = dust_limit
        self.dec_min = dec_min
        self.dec_max = dec_max
        if exgal_filter_balance is None:
            self.exgal_filter_balance = self.default_filter_balance
        else:
            self.exgal_filter_balance = self._normalize_filter_balance(exgal_filter_balance)
        
        # Set the detailed dust boundary 
        self.dust_exwfd = np.where((self.dec > self.dec_min) & (self.dec < self.dec_max) 
                                   & (self.dustmap < self.dust_limit), 1, 0)
        # Set the smoothed dust boundary using the original dustmap and smoothing it with gaussian PSF
        self.exwfd = np.where((self.dustmap < self.dust_limit), 1, 0)
        self.exwfd = hp.smoothing(self.exwfd, fwhm=np.radians(smoothing_beam))
        self.exwfd = np.where((self.dec > self.dec_min) & (self.dec < self.dec_max) 
                              & (self.exwfd>smoothing_cutoff), 1, 0)

        # Make per-filter maps for the footprint
        self.exwfd_maps = {}
        for f in self.filterlist:
            self.exwfd_maps[f] = self.exwfd * self.exgal_filter_balance[f]
        # Make these easy to access
        self.maps['exwfd'] = self.exwfd
        self.maps_perfilter['exwfd'] = self.exwfd_maps
        self.nvis['exwfd'] = nvis_exwfd

    def _set_magellanic_clouds(self, lmc_radius=10, smc_radius=5, 
                              nvis_mcs = 825,
                              mcs_filter_balance=None):
        # Define the magellanic clouds region
        if mcs_filter_balance is None:
            self.mcs_filter_balance = self.default_filter_balance
        else:
            self.mcs_filter_balance = self._normalize_filter_balance(mcs_filter_balance)

        self.mcs = np.zeros(hp.nside2npix(self.nside))
        # Define the LMC center and size
        lmc_ra = np.radians(80.893860)
        lmc_dec = np.radians(-69.756126)
        self.lmc_radius = np.radians(lmc_radius)
        # Define the SMC center and size
        smc_ra = np.radians(13.186588)
        smc_dec = np.radians(-72.828599)
        self.smc_radius = np.radians(smc_radius)
        # Define the LMC pixels
        dist_to_lmc = _angularSeparation(lmc_ra, lmc_dec, np.radians(self.ra), np.radians(self.dec))
        lmc_pix = np.where(dist_to_lmc < self.lmc_radius)
        self.mcs[lmc_pix] = 1
        # Define the SMC pixels
        dist_to_smc = _angularSeparation(smc_ra, smc_dec, np.radians(self.ra), np.radians(self.dec))
        smc_pix = np.where(dist_to_smc < self.smc_radius)
        self.mcs[smc_pix] = 1
        
        # Make per-filter maps for the footprint
        self.mcs_maps = {}
        for f in self.filterlist:
            self.mcs_maps[f] = self.mcs * self.mcs_filter_balance[f]
        self.maps['mcs'] = self.mcs
        self.maps_perfilter['mcs'] = self.mcs_maps
        self.nvis['mcs'] = nvis_mcs
    
    def __set_bulge_diamond(self, center_width, end_width, gal_long1, gal_long2):
        """
        Define a Galactic Bulge diamond-ish region.

        Parameters
        ----------
        center_width : float
            Width at the center of the galactic plane region.
        end_width : float
            Width at the remainder of the galactic plane region.
        gal_long1 : float
            Longitude at which to start the GP region.
        gal_long2 : float
            Longitude at which to stop the GP region.
            Order matters for gal_long1 / gal_long2!
        Returns
        -------
        np.ndarray
        """
        # Reject anything beyond the central width.
        bulge = np.where(np.abs(self.gal_lat) < center_width, 1, 0)
        # Apply the galactic longitude cuts, so that plane goes between gal_long1 to gal_long2.
        # This is NOT the shortest distance between the angles.
        gp_length = (gal_long2 - gal_long1) % 360
        # If the length is greater than 0 then we can add additional cuts.
        if gp_length > 0:
            # First, remove anything outside the gal_long1/gal_long2 region.
            bulge = np.where(((self.gal_lon - gal_long1) % 360) < gp_length, bulge, 0)
            # Add the tapers.
            # These slope from the center (gp_center @ center_width)
            # to the edges (gp_center + gp_length/2 @ end_width).
            half_width = gp_length / 2.
            slope = (center_width - end_width) / half_width
            # The 'center' can have a wrap-around 0 problem
            gp_center = (gal_long1 + half_width) % 360
            # Calculate the longitude-distance between any point and the 'center'
            gp_dist = (self.gal_lon - gp_center) % 360
            gp_dist = np.abs(np.where((gp_dist > 180), (180 - gp_dist) % 180, gp_dist))
            lat_limit = np.abs(center_width - slope * gp_dist)
            bulge = np.where((np.abs(self.gal_lat)) < lat_limit, bulge, 0)
        return bulge
    
    def _set_galactic_plane(self, dec_max=12,
                            center_width_A=11, end_width_A=4, gal_long1_A=330, gal_long2_A=30, 
                           center_width_B=15, end_width_B=5, gal_long1_B=250, gal_long2_B=100,
                           gal_lat_width_max=23,
                           nvis_gal_A=825, nvis_gal_B=300, nvis_gal_min=250, gal_filter_balance=None):
        if gal_filter_balance is None:
            self.gal_filter_balance = {'u': 0.04, 'g': 0.22, 'r': 0.24,
                                      'i': 0.24, 'z': 0.22, 'y': 0.05}
        else:
            self.gal_filter_balance = self._normalize_filter_balance(gal_filter_balance)
        self.gal_dec_max = dec_max

        # Set up central bulge
        self.bulge_A = self.__set_bulge_diamond(center_width=center_width_A, end_width=end_width_A,
                                               gal_long1=gal_long1_A, gal_long2=gal_long2_A)
        self.bulge_A = np.where(self.dec > self.gal_dec_max, 0, self.bulge_A)
        # And a secondary bulge-ish region to follow further stars
        self.bulge_B = self.__set_bulge_diamond(center_width=center_width_B, end_width=end_width_B,
                                               gal_long1=gal_long1_B, gal_long2=gal_long2_B)
        self.bulge_B = np.where(self.dec > self.gal_dec_max, 0, self.bulge_B)
        # Remove regions of these bulges which go further north than dec_max
        # Set up 'background' galactic plane visits 
        self.gp_bkgnd = np.where((np.abs(self.gal_lat) < gal_lat_width_max) & (self.dec < self.dec_max), 1, 0)
        # Remove the areas that overlap
        self.bulge_B = self.bulge_B - self.bulge_A
        self.gp_bkgnd = self.gp_bkgnd - self.bulge_A - self.bulge_B
        
        # Add them together
        self.gal = (self.gp_bkgnd * nvis_gal_min / nvis_gal_A 
                   + self.bulge_B * nvis_gal_B / nvis_gal_A
                   + self.bulge_A)
        
        # Make per-filter maps for the footprint
        self.gal_maps = {}
        for f in self.filterlist:
            self.gal_maps[f] = self.gal * self.gal_filter_balance[f]
        self.maps['gal'] = self.gal
        self.maps_perfilter['gal'] = self.gal_maps
        self.nvis['gal'] = nvis_gal_A
            
    def _set_nes(self, eclat_min=-15, eclat_max=10, eclip_dec_min=-10, eclip_ra_max=180, 
                 nvis_nes=400, nes_filter_balance=None):
        if nes_filter_balance is None:
            self.nes_filter_balance =  {'u': 0.0, 'g': 0.2, 'r': 0.3, 'i': 0.3, 'z': 0.2, 'y': 0.0}
        else:
            self.nes_filter_balance = self._normalize_filter_balance(nes_filter_balance)
        # NES ecliptic latitude values tend to be assymetric because NES goes so far north
        self.eclat_min = eclat_min
        self.eclat_max = eclat_max
        self.eclip_dec_min = eclip_dec_min
        self.eclip_ra_max = eclip_ra_max

        self.nes = np.where(((self.eclip_lat > self.eclat_min) | 
                             ((self.dec > self.eclip_dec_min) & (self.ra < self.eclip_ra_max)))
                             & (self.eclip_lat < self.eclat_max), 1, 0)
        self.nes_maps = {}
        for f in self.filterlist:
            self.nes_maps[f] = self.nes * self.nes_filter_balance[f]
        self.maps['nes'] = self.nes
        self.maps_perfilter['nes'] = self.nes_maps
        self.nvis['nes'] = nvis_nes
            
    def _set_scp(self, nvis_scp=120, dec_max=12, scp_filter_balance=None):
        if scp_filter_balance is None:
            self.scp_filter_balance = {'u': 0.17, 'g': 0.17, 'r': 0.17, 'i': 0.17, 'z': 0.17, 'y': 0.17}
        else:
            self.scp_filter_balance = self._normalize_filter_balance(scp_filter_balance)
        # Basically this is a fill-in so that we don't have any gaps below the max dec limit for the survey
        # I would expect most of this to be ignored
        self.scp = np.where(self.dec < dec_max, 1, 0)
        self.scp_maps = {}
        for f in self.filterlist:
            self.scp_maps[f] = self.scp * self.scp_filter_balance[f]
        self.maps['scp'] = self.scp
        self.maps_perfilter['scp'] = self.scp_maps
        self.nvis['scp'] = nvis_scp
    
    def set_maps(self):
        self.read_dustmap()
        self._set_exwfd()
        self._set_magellanic_clouds()
        self._set_galactic_plane()
        self._set_nes()
        self._set_scp()
    
    def combine_maps(self):
        map_order = ['exwfd', 'gal', 'mcs', 'nes', 'scp']

        total_perfilter = {}
        for f in self.filterlist:
            total_perfilter[f] = np.zeros(len(self.hpid), float)
        for m in map_order:
            if m in self.maps:
                for f in self.filterlist:
                    total_perfilter[f] = np.maximum(total_perfilter[f], 
                                                    self.maps_perfilter[m][f] * self.nvis[m])
        total = np.zeros(len(self.hpid), float)
        for f in self.filterlist:
            total += total_perfilter[f]
        return total, total_perfilter


def slice_wfd_area_quad(target_map, nslice=2, wfd_indx=None):
    """
    Make a fancy double striped map
    """
    nslice2 = nslice * 2

    wfd = target_map['r'] * 0
    if wfd_indx is None:
        wfd_indices = np.where(target_map['r'] == 1)[0]
    else:
        wfd_indices = wfd_indx
    wfd[wfd_indices] = 1
    wfd_accum = np.cumsum(wfd)
    split_wfd_indices = np.floor(np.max(wfd_accum)/nslice2*(np.arange(nslice2)+1)).astype(int)
    split_wfd_indices = split_wfd_indices.tolist()
    split_wfd_indices = [0] + split_wfd_indices

    return split_wfd_indices


def make_rolling_footprints(fp_hp=None, mjd_start=60218., sun_RA_start=3.27717639,
                            nslice=2, scale=0.8, nside=32, wfd_indx=None):

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

    split_wfd_indices = slice_wfd_area_quad(hp_footprints, nslice=nslice, wfd_indx=wfd_indx)

    roll = np.zeros(nslice)
    roll[-1] = 1
    for key in hp_footprints:
        temp = hp_footprints[key] + 0
        temp[wfd_indx] = 0
        fp_non_wfd.set_footprint(key, temp)

        for i in range(nslice):
            temp = hp_footprints[key] + 0
            temp[non_wfd_indx] = 0
            for j in range(nslice*2):
                indx = wfd_indx[split_wfd_indices[j]:split_wfd_indices[j+1]]
                temp[indx] = temp[indx] * roll[(i+j) % nslice]
            rolling_footprints[i].set_footprint(key, temp)

    result = Footprints([fp_non_wfd] + rolling_footprints)
    return result


def gen_greedy_surveys(nside=32, nexp=2, exptime=30., filters=['r', 'i', 'z', 'y'],
                       camera_rot_limits=[-80., 80.],
                       shadow_minutes=60., max_alt=76., moon_distance=30., ignore_obs='DD',
                       m5_weight=3., footprint_weight=0.75, slewtime_weight=3.,
                       stayfilter_weight=3., footprints=None):
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
                   stayfilter_weight=3., template_weight=12., footprints=None, u_nexp1=True):
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
                       wfd_footprint=None):
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

    args = parser.parse_args()
    survey_length = args.survey_length  # Days
    outDir = args.outDir
    verbose = args.verbose
    max_dither = args.maxDither
    illum_limit = args.moon_illum_limit
    nexp = args.nexp
    nslice = args.rolling_nslice
    scale = args.rolling_strength

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
    fileroot = os.path.basename(sys.argv[0]).replace('.py', '') + '_'
    file_end = 'v2.0_'

    sm = SurveyMap(nside=nside)
    sm.set_maps()
    final_tot, footprints_hp = sm.combine_maps()
    wfd_footprint = footprints_hp['r']*0
    # Bad bad magic number
    wfd_footprint[np.where(np.round(final_tot) >= 891)] = 1
    wfd_indx = np.where(wfd_footprint == 1)[0]

    normval = footprints_hp['r'][np.where(np.round(final_tot) == 891)].max()
    for key in footprints_hp:
        footprints_hp[key] = footprints_hp[key]/normval

    repeat_night_weight = None

    observatory = Model_observatory(nside=nside)
    conditions = observatory.return_conditions()

    footprints = make_rolling_footprints(fp_hp=footprints_hp, mjd_start=conditions.mjd_start,
                                         sun_RA_start=conditions.sun_RA_start, nslice=nslice, scale=scale,
                                         nside=nside, wfd_indx=wfd_indx)

    # Set up the DDF surveys to dither
    u_detailer = detailers.Filter_nexp(filtername='u', nexp=1)
    dither_detailer = detailers.Dither_detailer(per_night=per_night, max_dither=max_dither)
    details = [detailers.Camera_rot_detailer(min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit),
               dither_detailer, u_detailer]
    euclid_detailers = [detailers.Camera_rot_detailer(min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit),
                        detailers.Euclid_dither_detailer(), u_detailer]
    ddfs = generate_dd_surveys(nside=nside, nexp=nexp, detailers=details, euclid_detailers=euclid_detailers)

    greedy = gen_greedy_surveys(nside, nexp=nexp, footprints=footprints)
    blobs = generate_blobs(nside, nexp=nexp, footprints=footprints)
    twi_blobs = generate_twi_blobs(nside, nexp=nexp, footprints=footprints, wfd_footprint=wfd_footprint,
                                   repeat_night_weight=repeat_night_weight)
    surveys = [ddfs, blobs, twi_blobs, greedy]
    run_sched(surveys, survey_length=survey_length, verbose=verbose,
              fileroot=os.path.join(outDir, fileroot+file_end), extra_info=extra_info,
              nside=nside, illum_limit=illum_limit)
