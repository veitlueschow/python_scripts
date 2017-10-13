from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
import xz_tools as xz_tools
import numpy.ma as ma
class Args: pass

#cdo -sellevel,2038.5 -selindexbox,2170,2500,700,1120 /work/mh0256/m300522/data_storm/tape/60-90/rhopoto_60-90_tm.nc rhopoto_box.nc
#cdo -sellevel,2038.5 -selindexbox,2170,2500,700,1120 ../div_U+rho+.nc div_U+rho+_box.nc 
#cdo -sellevel,2038.5 -selindexbox,2170,2500,700,1120 ../dx_rhopoto.nc dx_rhopoto_box.nc 

rho = tools.netread_data('rhopoto_box.nc','rhopoto') # density
lat,lon,depth = tools.netread_grid('rhopoto_box.nc','lat','lon','depth_2')

divUrho = tools.netread_data("div_U+rho+_box.nc","div_Urho_eddy") # eddy flux divergence

dx_rho = tools.netread_data('dx_rhopoto_box.nc','dx_rhopoto') # zonal density gradient


high = 1e-8
v = np.linspace(-high, high, 100, endpoint=True)
f1 = plt.figure()
plt.title('test')
frho = plt.contour(lon,lat, dx_rho[0,:,:], 15,colors='k', linewidths=1)
plt.clabel(frho, fontsize=7, inline=1)
fdiv = plt.contourf(lon,lat,divUrho[0,:,:],v,extend="both")
plt.colorbar(format='%.0e')
plt.show()
