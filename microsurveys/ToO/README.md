Try to create a new ToO simulation.

Using https://arxiv.org/pdf/1812.04051.pdf as a point of figuring out a strategy.

Expect localizations on the order of 20-200 sq degrees, so 4 to 20 pointings needed to cover area. 

They want log-spacing of the observations at 1, 2, 4, 8 hours after merger.

They don't realize the filter wheel only holds 5 filters.

Maybe we can make a survey object that looks for new events, then generates scripted observations for each event. Give it a footprint that the center of the event has to be in. 



--------

Currently using circles with 74 sq deg size as the default. 
That takes ~22 pointings to cover. At 30s per exposure and filter swap, that's be ~15 min per filter. So, 75 minutes for the alert. And that is not taking additional dither positions. We could do something like 3 dither positions in r, then just 1 in the other filters. Then we're up to 105 minutes. This does make it tougher if we want log-spacing with the second observation 1 hour later. 

It's tough to automate the ToO followup. There can be filter thrashing if things are below the airmass limit. As is often the case, scheduling the sequence of followup would be better done with integer programming (or just manually).

This implementation of a followup code should only be considered a strawman for simulation purposes. Something more sophisticated is needed for production.

