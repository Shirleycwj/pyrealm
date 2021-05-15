import pyreadr
import numpy as np
import pandas as pd

def rda_read(fpath):
    rda_file = pyreadr.read_r(fpath)
    rda_data = np.array(list(rda_file.items()), dtype=object)[0, 1]
    return rda_data
