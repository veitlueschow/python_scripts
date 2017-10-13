import sys
sys.path.insert(0, '/home/mpim/m300522/mod_data')

import get_coastline as coasti
from netCDF4 import Dataset
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import decomp_plots as decplot

filename = "90er_mean_box"

[cline_s2, cline,coastline] = coasti.smooth_coastline(filename)

figure = plt.figure()
plt.plot(cline[:,1],cline[:,0])
plt.plot(cline_s2[:,1],cline_s2[:,0])
figure.savefig("coastline.png")


