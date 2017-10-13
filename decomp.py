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
po_ 	= fh.variables["po"]
uko_ 	= fh.variables["uko"]
vke_ 	= fh.variables["vke"]
lat_2_ 	= fh.variables["lat_2"]
lon_2_ 	= fh.variables["lon_2"]
lat_3_ 	= fh.variables["lat_3"]
lon_3_ 	= fh.variables["lon_3"]

inn = Args()
gulf = Args()

inn.uko = uko_[0,0,100:1100,:].copy()
inn.rho = rho_[0,0,100:1100,:].copy()
inn.po	= po_[0,0,100:1100,:].copy()
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
#global minimal
#global minimaly
gulf.minimal = np.zeros((inn.uko.shape))
gulf.minimaly = np.zeros((inn.uko.shape))
for y in range(min(inn.vke.shape[0],inn.uko.shape[0])):
	for x in range(min(inn.vke.shape[1],inn.uko.shape[1])):
		if inn.uko.mask[y,x] == False and inn.vke.mask[y,x] == False and inn.lon_2[y,x] >= (cline_s[y,1]): # positive distances
			dist = abs(inn.lon_2[y,x] -cline_s[y,1])
			[minimaldist, ymin] = coasti.check_neighbours(x,y,cline_s,dist,inn.lon_2,inn.lat_2)
			gulf.minimal[y,x] = minimaldist
			gulf.minimaly[y,x] = ymin
			#print minimaldist
			if (ymin-1) >0 and (ymin+1) < min(inn.vke.shape[0],inn.uko.shape[0]):
				nx = normal[ymin-1,1]/4 + normal[ymin,1]/2 + normal[ymin+1,1]/4
				ny = normal[ymin-1,0]/4 + normal[ymin,0]/2 + normal[ymin+1,0]/4
				tx = tangent[ymin-1,1]/4 + tangent[ymin,1]/2 + tangent[ymin+1,1]/4
				ty = tangent[ymin-1,0]/4 + tangent[ymin,0]/2 + tangent[ymin+1,0]/4
			else:
				nx = normal[ymin,1]
				ny = normal[ymin,0]
				tx = tangent[ymin,1]
				ty = tangent[ymin,0]
			if np.sign(ny) == np.sign(inn.vke[y,x]):
				gulf.uper[y,x] = nx*inn.uko[y,x] + ny*inn.vke[y,x]
			else:
				gulf.uper[y,x] = nx*inn.uko[y,x] - ny*inn.vke[y,x]
			if np.sign(tx) == np.sign(inn.uko[y,x]):
				gulf.upar[y,x] = ty*inn.vke[y,x] + tx*inn.uko[y,x]
			else:
				gulf.upar[y,x] = ty*inn.vke[y,x] - tx*inn.uko[y,x]
			
			
		elif inn.uko.mask[y,x] == False and inn.vke.mask[y,x] == False and inn.lon_2[y,x] < (cline_s[y,1]) and inn.lon_2[y,x] >= (cline_s[y,1]-1): # negative distances
			dist = abs(inn.lon_2[y,x] -cline_s[y,1])
			[minimaldist, ymin] = coasti.check_neighbours(x,y,cline_s,dist,inn.lon_2,inn.lat_2)
			gulf.minimal[y,x] = -minimaldist
			#print "negative", x,y, -minimaldist
			gulf.minimaly[y,x] = ymin
			#print minimaldist
			if (ymin-1) >0 and (ymin+1) < min(inn.vke.shape[0],inn.uko.shape[0]):
				nx = normal[ymin-1,1]/4 + normal[ymin,1]/2 + normal[ymin+1,1]/4
				ny = normal[ymin-1,0]/4 + normal[ymin,0]/2 + normal[ymin+1,0]/4
				tx = tangent[ymin-1,1]/4 + tangent[ymin,1]/2 + tangent[ymin+1,1]/4
				ty = tangent[ymin-1,0]/4 + tangent[ymin,0]/2 + tangent[ymin+1,0]/4
			else:
				nx = normal[ymin,1]
				ny = normal[ymin,0]
				tx = tangent[ymin,1]
				ty = tangent[ymin,0]
			if np.sign(ny) == np.sign(inn.vke[y,x]):
				gulf.uper[y,x] = nx*inn.uko[y,x] + ny*inn.vke[y,x]
			else:
				gulf.uper[y,x] = nx*inn.uko[y,x] - ny*inn.vke[y,x]
			if np.sign(tx) == np.sign(inn.uko[y,x]):
				gulf.upar[y,x] = ty*inn.vke[y,x] + tx*inn.uko[y,x]
			else:
				gulf.upar[y,x] = ty*inn.vke[y,x] - tx*inn.uko[y,x]

# Set all the missing values to zero
for x in range(min(inn.vke.shape[1],inn.uko.shape[1])):
	for y in range(min(inn.vke.shape[0],inn.uko.shape[0])):
		if abs(gulf.upar[y,x]) > 100.:
			gulf.upar[y,x] = 0.
		if abs(gulf.uper[y,x]) > 100.:
			gulf.uper[y,x] = 0.

gulf.upar = -gulf.upar
gulf.uper = -gulf.uper
#[across,acr,along,data1,data2,data3,data4] = decplot.prepare_data(inn,gulf,cline_s)
[along,across,acr,raw,full] = decplot.prepare_data(inn,gulf,cline_s)
decplot.rho_along_across(along,across,full.rho,raw.ind,350,450)
decplot.upar_along_across(along,across,full.upar,350,450)
decplot.po_along_across(along,across,full.po,raw.ind,350,450)

