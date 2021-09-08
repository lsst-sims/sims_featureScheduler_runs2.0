Setting up runs where we try and get tripples durring the night

Currently:

* Takes a regular blob pair, then schedules observations for 2.5 hours later
* Currently forcing first observations to be at HA < -2.5 hours. Could try to make this more sophisticated. 
* probably need more possible filter pairs so we don't take bad visits in bright time
* Could add some constraint on requiring enough area before starting a blob that will attempt triplets
