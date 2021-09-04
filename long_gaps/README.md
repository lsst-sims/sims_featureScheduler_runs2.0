Try to make some observations with gaps of 2-14 hours

Here's a plan:
* when it's not twilight, and haven't observed the start in the night, can draw a gap of between 2hours to (end of night-block time). Make a healpix map of things visible now, and visible at second time. Maybe even approximate the 5-sigma depth in the future. Generate a blob-pair, then add the long gap to the scheduled observations for later.

make a wrapper survey object, that holds some blob survey(s), and if the blob returns something, then make it scheduled for later too.

--------

How it's working

survey object has a blob and scripted component:

The blob can be g+r, r+i, or i+z. pairs separated by 33 minutes. The usual footprint and all the other standard blob parameters with additional constaints that it must be after evening twilight (18 degree) but not more than 30 minutes post-twilight. Also, the potential sky area with hour angle > -3.5 hours is masked off. 

When the blob meets it's observing criteria, it generates a list of around 100 observations (50 in filter 1, then those 40 repeated in filter 2). We take the first 50, reverse their order, and record them as scheduled observations for 2-7 hours in the future (picking the actual time from a uniform random distribution.).  


Potential issues:

* Not much look-ahead. So we may end up getting a g-band 3 hours later, after the moon has risen. Can add in logic to see if the moon will rise, then select the redder filter. Or only allow g+r if the moon is down and will stay down.
* Not paying attention to if things will track into the zenith exclusion zone.
* Forcing the blob to be rising in the east. I suppose we could soften this a bit if we wanted to. But I think this way the long gap observations will be near the meridian or in the west, so maybe it doesn't matter.
* Do I need to seed the random number generator? Should probably controll all the seeding in the conditions object.
* Is this going to smear out the rolling because there's not enough area after the HA mask?
* This could also over-expose the SCP, since it'll be available the most.
* Do we want to do this on only the WFD? Do we want the coverage of long-gap observations to be somewhat spatially even?