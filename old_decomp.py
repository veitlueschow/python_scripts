import sys
sys.path.insert(0, '/home/mpim/m300522/mod_data')

import get_coastline_part as coasti
from netCDF4 import Dataset
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import decomp_plots as decplot
from classes import cart
from imp import reload
class Args: pass

filename = "90s_box_mean"

fh = Dataset(filename,mode='r')
rho_ 	= fh.variables["rhopoto"]
uko_ 	= fh.variables["uko"]
vke_ 	= fh.variables["vke"]
lat_2_ 	= fh.variables["lat_2"]
lon_2_ 	= fh.variables["lon_2"]
lat_3_ 	= fh.variables["lat_3"]
lon_3_ 	= fh.variables["lon_3"]

inn = Args()

inn.uko = uko_[0,0,100:1100,:].copy()
inn.rho = rho_[0,0,100:1100,:].copy()
inn.vke = vke_[0,0,100:1100,:].copy()
inn.lat_2 = lat_2_[100:1100,:].copy()
inn.lon_2 = lon_2_[100:1100,:].copy()
inn.lat_3 = lat_3_[100:1100,:].copy()
inn.lon_3 = lon_3_[100:1100,:].copy()

[cline,cline_s] = coasti.smooth_coastline(filename,111)
[tangent, normal] = coasti.compute_normals(cline_s)

np.ma.set_fill_value(inn.uko,0.)

[cline,cline_s] = coasti.smooth_coastline(filename,111)

gulf.uper = np.zeros((inn.uko.shape))
gulf.upar = np.zeros((inn.uko.shape))
global minimal
global minimaly
minimal = np.zeros((inn.uko.shape))
minimaly = np.zeros((inn.uko.shape))
for y in range(min(vke.shape[0],uko.shape[0])):
	for x in range(min(vke.shape[1],uko.shape[1])):
		if uko.mask[y,x] == False and vke.mask[y,x] == False and lon_2[y,x] >= (cline_s[y,1]): # positive distances
			dist = abs(lon_2[y,x] -cline_s[y,1])
			[minimaldist, ymin] = coasti.check_neighbours(x,y,cline_s,dist,lon_2,lat_2)
			minimal[y,x] = minimaldist
			minimaly[y,x] = ymin
			#print minimaldist
			if (ymin-1) >0 and (ymin+1) < min(vke.shape[0],uko.shape[0]):
				nx = normal[ymin-1,1]/4 + normal[ymin,1]/2 + normal[ymin+1,1]/4
				ny = normal[ymin-1,0]/4 + normal[ymin,0]/2 + normal[ymin+1,0]/4
				tx = tangent[ymin-1,1]/4 + tangent[ymin,1]/2 + tangent[ymin+1,1]/4
				ty = tangent[ymin-1,0]/4 + tangent[ymin,0]/2 + tangent[ymin+1,0]/4
			else:
				nx = normal[ymin,1]
				ny = normal[ymin,0]
				tx = tangent[ymin,1]
				ty = tangent[ymin,0]
			if np.sign(ny) == np.sign(vke[y,x]):
				uper[y,x] = nx*uko[y,x] + ny*vke[y,x]
			else:
				uper[y,x] = nx*uko[y,x] - ny*vke[y,x]
			if np.sign(tx) == np.sign(uko[y,x]):
				upar[y,x] = ty*vke[y,x] + tx*uko[y,x]
			else:
				upar[y,x] = ty*vke[y,x] - tx*uko[y,x]
			
			
		elif uko.mask[y,x] == False and vke.mask[y,x] == False and lon_2[y,x] < (cline_s[y,1]) and lon_2[y,x] >= (cline_s[y,1]-1): # negative distances
			dist = abs(lon_2[y,x] -cline_s[y,1])
			[minimaldist, ymin] = coasti.check_neighbours(x,y,cline_s,dist,lon_2,lat_2)
			minimal[y,x] = -minimaldist
			#print "negative", x,y, -minimaldist
			minimaly[y,x] = ymin
			#print minimaldist
			if (ymin-1) >0 and (ymin+1) < min(vke.shape[0],uko.shape[0]):
				nx = normal[ymin-1,1]/4 + normal[ymin,1]/2 + normal[ymin+1,1]/4
				ny = normal[ymin-1,0]/4 + normal[ymin,0]/2 + normal[ymin+1,0]/4
				tx = tangent[ymin-1,1]/4 + tangent[ymin,1]/2 + tangent[ymin+1,1]/4
				ty = tangent[ymin-1,0]/4 + tangent[ymin,0]/2 + tangent[ymin+1,0]/4
			else:
				nx = normal[ymin,1]
				ny = normal[ymin,0]
				tx = tangent[ymin,1]
				ty = tangent[ymin,0]
			if np.sign(ny) == np.sign(vke[y,x]):
				uper[y,x] = nx*uko[y,x] + ny*vke[y,x]
			else:
				uper[y,x] = nx*uko[y,x] - ny*vke[y,x]
			if np.sign(tx) == np.sign(uko[y,x]):
				upar[y,x] = ty*vke[y,x] + tx*uko[y,x]
			else:
				upar[y,x] = ty*vke[y,x] - tx*uko[y,x]

# Set all the missing values to zero
for x in range(min(vke.shape[1],uko.shape[1])):
	for y in range(min(vke.shape[0],uko.shape[0])):
		if abs(upar[y,x]) > 100.:
			upar[y,x] = 0.
		if abs(uper[y,x]) > 100.:
			uper[y,x] = 0.

upar = -upar
uper = -uper
[across,acr,along,data1,data2,data3,data4] = decplot.prepare_data(uko,vke,upar,uper,minimal,cline_s,minimaly)
[a,b1,b2,raw_upar,full_upar,raw_rho,full_rho] = decplot.preparation(uko,upar,minimal,cline_s,minimaly,rho)
#c = decplot.along_across(along,across,data1)
