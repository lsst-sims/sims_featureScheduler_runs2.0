import numpy as np
import pandas as pd
from rubin_sim.maf import RunComparison


stats = RunComparison(baseDir='.')()
stats.to_pickle('combined_stats.pkl')
