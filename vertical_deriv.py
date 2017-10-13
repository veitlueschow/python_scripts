from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
import numpy.ma as ma

# This function computes a vertical derivative, e.g. to get the vertical velocity shear


# Get layer thickness!

filename='/work/mh0256/m300522/meta_storm/ddpo.nc'
fh_ddpo= Dataset(filename,mode='r')

ddpo_ = fh_ddpo-variables["ddpo"]
ddpo = ddpo_[:,:,:].copy()
