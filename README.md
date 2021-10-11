# sims_featureScheduler_runs2.0
v2.0 of the featureScheduler experiments


# Description of runs

# Baseline

* new footprint
* Pairs in non-twilight taken with ~33 minute separation (u+g, u+r, g+r, r+i, i+z, z+y, y+y)
* Pairs in twilight taken with ~15 min seperation (r+i, i+z, z+y, y+y)
* 1x30s exposures in u, 2x15s exposures in all other filters
* a half-sky rolling cadence in the exgal WFD area.

Filter distribution of:  {'u': 0.07, 'g': 0.09, 'r': 0.22, 'i': 0.22, 'z': 0.20, 'y': 0.20}


# Bluer

Testing two bluer filter distributions:

{'u': 0.07, 'g': 0.12, 'r': 0.21, 'i': 0.21, 'z': 0.19, 'y': 0.20}
{'u': 0.08, 'g': 0.11, 'r': 0.21, 'i': 0.21, 'z': 0.19, 'y': 0.20}

# DDF

Varying the amount of time given to observing the DDFs. Not running everything requested by the SCOC since we plan on doing DDF scheduling differently in the future.

# Long gaps

Takes observations in the east g+r, r+i, or i+z with 33 min separation. Then take a third observation 2-7 hours later.

Vary how often we attempt to observe triples (from every night, to every 7th night).

Adding an additional set where it doesn't start attempting long gaps until year 5 (tagged with "delayed1827"). Note, that means we are changing strategy in the middle of a rolling cadence season. There are non-rolling seasons at the end of the survey so the footprint should have plenty of time to recover, but could impact sciences that benefit from rolling.

# Long gaps no pair

Like long gaps, but the initial observation is unpaired. Again, varying the fraction of nights where long gaps are attempted and varying running for the whole survey and starting at year 5. 

# Long u

long_u1:  Increase u-band exposure time to 50s (leaving number of u band the same)
long_u2:  Increase u-band to 50s, decrease the number of exposures so total u-exposure time is similar to baseline.

# Microsurveys

## carina

Does intensive observations of carina nebula for 1 week each year.

## local_gals

Observing deeper on 10 nearby galaxies. Increases the g,r,i visits. Tried increasing the number of visits at three levels (50% more g, 100% more g, and 150% more g) with smaller increases in the other filters.  Note that the cadence note used minion_1016 as a reference point, so they were a bit optimistic on what the base coadded depths would be. May not be as feasible with the more realistic weather we have now. 

## multi_short

Takes sequences of 4 5-second exposures. Shoots for 12 total short exposures per year.

## north_stripe

Add a northern extension to the footprint. Probably just to have image templates to chance ToOs anywhere on the sky.

## short_exp

Takes 5 second exposures in all filters in year 1. 

## twilight_neo

Using twilight time for short exposure NEO search. Varying what fraction of nights are used for NEO searching. 

## virgo_cluster

Added the virgo galaxy cluster to the WFD area. This might be a no-brainer add, looks very low impact.


# noroll

Like baseline, but no rolling cadence.

# presto

taking pairs in the east (g+r, r+i, i+z) or (g+i, r+z, i+y) labeled as "mix", then taking a third observation 1.5-4.0 hours later. 


# retro

retro_baseline_v2.0_10yrs.db: Similar to previous baseline surveys. Mostly here to be used as a comparison run to check changes between previous baseline strategy and new runs.

baseline_retrofoot_v2.0_10yrs.db:  The v2.0 baseline settings, but with the classic footprint. Rolling might do some wacky things here because it's expecting a different footprint.

# rolling

Rolling cadences with half, or third of the sky. Trying 0.5 or 0.9 rolling strength. 


# rolling_all_sky

Doing rolling in the WFD, splitting the bulge in half and rolling, and rolling in the dusty plane. Probably needs a slightly more sophisticated way of splitting the dusty region (now it's splitting by area and not integrated etendue like it probably should be), but it's few enough visits that it's not a major impact.

# rolling_bulge_6

The standard baseline rolling in the WFD, and now the bulge is divided into 6 stripes which roll. During the non-rolling time, the bulge is observed like the baseline (I think). Once rolling starts, the bulge area is so small it is not feasible to do repeated contiguous regions, so the bulge observations are not paired. 

The rolling in the bulge is strong enough it may be pushing some observations to fairly high airmass compared to the baseline.  


# rolling bulge

Like the baseline, but split the MW bulge in half and roll there as well. 

# roll_early

Similar to the baseline, but begin rolling earlier so an additional season of rolling gets completed. 

# rolling_six

Rolling cadences where 1/6th of the sky is on for any season. 0.5 or 0.9 rolling strength.

# vary_gp

Vary the amount of time spent observing the galactic plane.

# vary_nes

Vary the amount of time spent observing the NES. Looks like some strange things can happen with the footprint if the NES is turned up very high.

