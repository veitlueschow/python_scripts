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

#cdo selindexbox,2170,2500,700,1120 den1.nc den1_box.nc
#cdo selindexbox,2170,2500,700,1120 den2.nc den2_box.nc
#cdo selindexbox,2170,2500,700,1120 den3.nc den3_box.nc
#cdo selindexbox,2170,2500,700,1120 den4.nc den4_box.nc
#cdo selindexbox,2170,2500,700,1120 ../div_U+rho+.nc div_U+rho+_box.nc

divUrho = tools.netread_data('div_U+rho+_box.nc','div_Urho_eddy')
den1 = tools.netread_data('den1_box.nc','den1')
den2 = tools.netread_data('den2_box.nc','den2')
den3 = tools.netread_data('den3_box.nc','den3')
den4 = tools.netread_data('den4_box.nc','den4') 
lat,lon,depth = tools.netread_grid('den1_box.nc','lat_2','lon_2','depth')

eps = 2e-15
tmp = np.zeros((den1.shape))
for i in range(den1.shape[2]):
	print(i)
	for j in range(den1.shape[1]):
		#for k in range(den1.shape[0]):
		for k in range(56,57):	
			if np.abs(den1[k,j,i]) <= eps:
				den1[k,j,i] = 0.
			if np.abs(den2[k,j,i]) <= eps:
				den2[k,j,i] = 0.
			if np.abs(den3[k,j,i]) <= eps:
				den3[k,j,i] = 0.
			if np.abs(den4[k,j,i]) <= eps:
				den4[k,j,i] = 0.
			if np.abs(den1[k,j,i]) <= eps and np.abs(den2[k,j,i]) <= eps and np.abs(den3[k,j,i]) <= eps and np.abs(den4[k,j,i]) <= eps:
				tmp[k,j,i] = 1

kappa = np.zeros((den1.shape))
for i in range(den1.shape[2]):
	print(i)
	for j in range(den1.shape[1]):
		#for k in range(den1.shape[0]):
		for k in range(56,57):
			if tmp[k,j,i] == 0.:
				kappa[k,j,i] = -divUrho[k,j,i] / ( den1[k,j,i] + den2[k,j,i] - den3[k,j,i] - den4[k,j,i])

high=1e4
v = np.linspace(-high, high, 100, endpoint=True)
plt.contourf(lon,lat,kappa[56,:,:],v,extend="both")
plt.show()

