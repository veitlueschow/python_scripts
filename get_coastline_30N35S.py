from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
import tools as tools
import write_netCDF as write

# For this, a box with the command "cdo selindexbox,1700,3000,500,1350" was used ( I think). Like this, the values in clines() are to be read with reference to 1700

filename='box_vke.nc'
#filename='vke_box.nc'
fh_vke 		= Dataset(filename,mode='r')

vke_ 	= fh_vke.variables["vke"]
depth_	= fh_vke.variables["depth_4"]
vke	= vke_[0,:,:,:].copy()
depth 	= depth_[:].copy()

## FIRST get one coastline as reference ---------------------------------------
k=20 # starting level
var = vke[k,:,:].copy()

coastline_ = np.zeros((var.shape[0],var.shape[1]))
coastline = np.zeros((var.shape[0],var.shape[1]))
for y in range(var.shape[0]):
	for x in range(var.shape[1]):
		if var.mask[y,x] == True:
			coastline_[y,x] = 1

for y in range(0,150):
	for x in range(400,10,-1):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-3] == 1 and coastline_[y,x-2] == 1 and coastline_[y,x-1] == 1 and coastline_[y,x] == 1 and coastline_[y,x+1] == 0:
			coastline[y,x] = 1
			
for y in range(150,var.shape[0]):
	for x in range(800,10,-1):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-3] == 1 and coastline_[y,x-2] == 1 and coastline_[y,x-1] == 1 and coastline_[y,x] == 1 and coastline_[y,x+1] == 0:
			coastline[y,x] = 1
			
## FIND ENVELOPE for first level (k=20, depth=271 m)

cline_x = np.zeros((var.shape[0]))
for y in range(var.shape[0]):
	tmp = np.flatnonzero(coastline[y,:])
	cline_x[y] = tmp
	
cline_y = np.linspace(0,850,num=851)

q_u = np.zeros(cline_x.shape)
u_x = [cline_x[0],]
u_y = [cline_y[0],]

for y in range(1,99):
	u_x.append(cline_x[y])
	u_y.append(cline_y[y])

for y in range(100,300):
	tmp = 0.66*y + 257
	if cline_x[y] >= tmp:
		u_x.append(cline_x[y])
		u_y.append(cline_y[y])

for y in range(301,851):
	u_x.append(cline_x[y])
	u_y.append(cline_y[y])

u_p = interp1d(u_y,u_x, kind = 'cubic',bounds_error = False, fill_value=0.0)

for y in range(len(cline_x)):
    q_u[y] = u_p(cline_y[y])

window_len = 71

tmp = tools.smooth(q_u,window_len,"hanning")
smooth_x= tmp[window_len/2-1:cline_x.shape[0]+window_len/2-1]

smooth_x = np.around(smooth_x)
smooth_x = np.int64(smooth_x) # Now contains the x-positions of the coastline as integer

clines = np.zeros((60,851)) # 60 clines, because the first 20 are not considered
clines[0,:] = smooth_x

## NOW go on finding the other coastline, starting from smooth_x
coastlines = np.zeros((60,var.shape[0],var.shape[1]))
coastlines_ = np.zeros((60,var.shape[0],var.shape[1]))

for k in range(21,80):
#for k in range(20,80):
	var = vke[k,:,:].copy()
	ik=k-20
	smooth_x = clines[ik-1,:].copy()
	smooth_x = np.int64(smooth_x)

	coastline_ = np.zeros((var.shape[0],var.shape[1]))
	coastline = np.zeros((var.shape[0],var.shape[1]))
	
	for y in range(var.shape[0]):
		for x in range(var.shape[1]):
			if var.mask[y,x] == True:
				coastline_[y,x] = 1
	
	for y in range(var.shape[0]):
		for x in range(smooth_x[y]-5,smooth_x[y]+20):
			if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x] == 1 and coastline_[y,x+1] == 0 and coastline_[y,x+2] == 0:
				coastline[y,x] = 1
	
	coastlines[ik,:,:] = coastline[:,:]
	coastlines_[ik,:,:] = coastline_[:,:]
				
	u_x=[smooth_x[0],]
	u_y=[0.,]
		
	for y in range(1,var.shape[0]):
		if np.count_nonzero(coastline[y,:]) != 0:
			tmp = np.flatnonzero(coastline[y,:])
			u_x.append(tmp[0])
			u_y.append(y)
	
	if u_y[-1] <= 849:
		u_x.append(smooth_x[850])
		u_y.append(850.0)
			
	u_p = interp1d(u_y,u_x, kind = 'cubic',bounds_error = False, fill_value=0.0)
	
	for y in range(len(cline_y)):
	    q_u[y] = u_p(cline_y[y])
	
	window_len = 71
	
	tmp = tools.smooth(q_u,window_len,"hanning")
	smooth_x2= tmp[window_len/2-1:cline_y.shape[0]+window_len/2-1]
	
	smooth_x2 = np.around(smooth_x2)
	smooth_x2 = np.int64(smooth_x2) # Now contains the x-positions of the coastline as integer
	
	clines[ik,:] = smooth_x2
	
## CREATE coordinates for write_netCDF
depth2 = depth[20:80].copy()
y,z = np.meshgrid(cline_y,depth2)

write.write_depth_meridional("clines_30N35S.nc",clines,"coastlines",y,z)

## -------------------------- end smooth_coastline ------------
