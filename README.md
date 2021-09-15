# sims_featureScheduler_runs2.0
v2.0 of the featureScheduler experiemnts


# Description of runs

# Baseline

* new footprint
* Pairs in non-twilight taken with ~33 minute seperation (u+g, u+r, g+r, r+i, i+z, z+y, y+y)
* Pairs in twilight taken with ~15 min seperation (r+i, i+z, z+y, y+y)
* 1x30s exposures in u, 2x15s exposures in all other filters
* a half-sky rolling cadence in the exgal WFD area.

Filter distribution of:  {'u': 0.07, 'g': 0.09, 'r': 0.22, 'i': 0.22, 'z': 0.20, 'y': 0.20}


# Bluer

Testing two bluer filter distributions:

{'u': 0.07, 'g': 0.12, 'r': 0.21, 'i': 0.21, 'z': 0.19, 'y': 0.20}
{'u': 0.08, 'g': 0.11, 'r': 0.21, 'i': 0.21, 'z': 0.19, 'y': 0.20}

# DDF

Varying the amount of time given to observing the DDFs.

# Long gaps

Takes observations in the east g+r, r+i, or i+z with 33 min separation. Then take a third observation 2-7 hours later.

Vary how often we attempt to observe triples (from every night, to every 7th night) 

# Long gaps no pair

Like long gaps, but the initial observation is unpaired. 

# Long u

long_u1:  Increase u-band exposure time to 50s (leaving number of u band the same)
long_u2:  Increase u-band to 50s, decrease the number of exposures so total u-exposure time is similar to baseline.

# Microsurveys

## carina

Does intensive observations of carina nebula for 1 week each year.

## short_exp

Takes 5 second exoposures in all filters in year 1.

## multi_short

Takes sequences of 4 5-second exposures. Shoots for 12 total short exposures per year.

# noroll

Like baseline, but no rolling cadence.

# presto

taking pairs in the east (g+r, r+i, i+z), then taking a third observation ~2.5 hours later. 


# retro

Similar to previous baseline surveys. Mostly here to be used as a comparison run to check changes between previous baseline strategy and new runs.

# rolling

Rolling cadences with half, third, or sixth of the sky. Trying 0.5 or 0.9 rolling strength.

# rolling bulge

Like the baseline, but split the MW bulge in half and roll there as well.

# vary_gp

Vary the amount of time spent observing the galactic plane.

# vary_nes

Vary the amount of time spent observing the NES.

