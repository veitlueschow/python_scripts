from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
class Args: pass

  
# for "box" use : cdo selindexbox,1950,2600,500,1350

# Get the different data types...

rho = tools.netread_data('rhopoto_box.nc','rhopoto') # density
lat,lon,depth = tools.netread_grid('rhopoto_box.nc','lat','lon','depth_2')
divUrho = tools.netread_data("div_U+rho+_box.nc","div_Urho_eddy") # eddy flux divergence
#dx_rho = tools.netread_data('dx_rhopoto_box.nc','dx_rhopoto') # zonal density gradient
vke = tools.netread_data('vke_box.nc','vke') # vke
wrho = tools.netread_data('w+rho+_box.nc','wrho_eddy') # wrho

# Get coastlines!

filename='/work/mh0256/m300522/meta_storm/clines_30N35S.nc'
fh_clines= Dataset(filename,mode='r')

coastlines_ = fh_clines.variables["coastlines"]
clines = coastlines_[:,:].copy()
clines = clines - 250 # Correct for the fact the clines was computed from a different box. Like this, clines refers to the current setting!

# compute the tangent and normal component at each latitude and all depth levels from 20 to 70
#tangent,normal = decomp_all.get_normals(clines,lat,lon)


# Now choose a latitude at which you want to have a x-z-slice
latitude = 13.


tmp = np.where(np.around(lat[:,400],decimals=1)==latitude)[0]
ilat = tmp[0] # ilat is the latitude index that I want to look for

# Choose the depth levels that you are interested in
ztop 	= 1000.
zbot	= 3000.

tmp1 = min(depth[:], key=lambda x:abs(x-ztop))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
ktop = tmp[0]

tmp1 = min(depth[:], key=lambda x:abs(x-zbot))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
kbot = tmp[0]

# Take the zonal boundaries from clines
if ktop >= 20:
	ktop_ = ktop -20
else:
	ktop_ = 0
	
iw = np.int64(clines[ktop_,ilat]) -3
ie = iw + 40

# Now construct the cartesian grid
x,z = np.meshgrid(lon[ilat,iw:ie],depth[ktop:kbot])

# Extract data U want to plot
#data_vrho 	= vrho[ktop3:kbot3,ilat3,iw:ie].copy()
data_rho 	= rho[ktop:kbot,ilat,iw:ie].copy()
data_vke 	= vke[ktop:kbot,ilat,iw:ie].copy()
data_divUrho= divUrho[ktop:kbot,ilat,iw:ie].copy()
data_wrho= wrho[ktop:kbot,ilat,iw:ie].copy()
#data_us= us[ktop:kbot,ilat,iw:ie].copy()
#data_U_grad_rho= U_grad_rho[ktop:kbot,ilat,iw:ie].copy()
#data_upar=upar[ktop:kbot,ilat,iw:ie].copy()

#high=5.e-9
#a = aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_U_grad_rho,high)
#a.canvas.set_window_title('Mean advection of density ' + str(latitude))
high =2.e-8
b= aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_divUrho,high,"divUrho")
b.canvas.set_window_title('Eddy Fluxes of Density ' + str(latitude))

#high =3.e-7
#a= aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_wrho,high,"wrho")
#a.canvas.set_window_title('wrho ' + str(latitude))
