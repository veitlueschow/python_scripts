from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
import numpy.ma as ma
class Args: pass

  

# Get the different data types...

#cdo selindexbox,1950,2600,705,775 ../div_U+rho+.nc div_U+rho+_stripe.nc
#cdo selindexbox,1950,2600,705,775 ../v+rho+.nc v+rho+_stripe.nc
#cdo selindexbox,1950,2600,705,775 ../u+rho+.nc u+rho+_stripe.nc
#cdo selindexbox,1950,2600,705,775 ../U_grad_rho.nc U_grad_rho_stripe.nc
#cdo selindexbox,1950,2600,705,775 /work/mh0256/m300522/data_storm/tape/60-90/vke_60-90_tm.nc vke_stripe.nc
#cdo selindexbox,1950,2600,705,775 /work/mh0256/m300522/data_storm/tape/60-90/rhopoto_60-90_tm.nc rhopoto_stripe.nc

filename='v+rho+_stripe.nc'
fh_vrho= Dataset(filename,mode='r')

vrho_ 	= fh_vrho.variables["vrho_eddy"]
lat_3_ 	= fh_vrho.variables["lat_3"]
lon_3_ 	= fh_vrho.variables["lon_3"]
depth46_= fh_vrho.variables["depth_2"]

vrho	= vrho_[0,:,:,:].copy()
lon3	= lon_3_[:,:].copy()
lat3	= lat_3_[:,:].copy()
depth3	= depth46_[:].copy()

filename='vke_stripe.nc'
fh_vke= Dataset(filename,mode='r')
vke_	= fh_vke.variables["vke"]
vke	= vke_[0,:,:,:].copy()


#filename='uko_stripe.nc'
#fh_uko= Dataset(filename,mode='r')
#uko_	= fh_uko.variables["uko"]
#uko	= uko_[0,:,:,:].copy()

filename='u+rho+_stripe.nc'
fh_urho= Dataset(filename,mode='r')

urho_ 	= fh_urho.variables["urho_eddy"]
lat_2_ 	= fh_urho.variables["lat_2"]
lon_2_ 	= fh_urho.variables["lon_2"]
depth41_= fh_urho.variables["depth_2"]

urho	= urho_[0,:,:,:].copy()
lon2	= lon_2_[:,:].copy()
lat2	= lat_2_[:,:].copy()
depth2	= depth41_[:].copy()


filename='rhopoto_stripe.nc'
fh_rho= Dataset(filename,mode='r')

rho_ 	= fh_rho.variables["rhopoto"]
lat_ 	= fh_rho.variables["lat"]
lon_ 	= fh_rho.variables["lon"]
depth14_= fh_rho.variables["depth_2"]

rho	= rho_[0,:,:,:].copy()
lon	= lon_[:,:].copy()
lat	= lat_[:,:].copy()
depth	= depth14_[:].copy()

filename="div_U+rho+_stripe.nc"
fh_divUrho=Dataset(filename,mode='r')
divUrho_ 	= fh_divUrho.variables["div_Urho_eddy"]
divUrho		= divUrho_[0,:,:,:].copy()

filename="U_grad_rho_stripe.nc"
fh_mean_ad=Dataset(filename,mode='r')
U_grad_rho_ 	= fh_mean_ad.variables["mean_advection_rho"]
U_grad_rho		= U_grad_rho_[0,:,:,:].copy()



# Choose the depth levels that you are interested in
ztop 	= 1000.
zbot	= 3500.

tmp1 = min(depth3[:], key=lambda x:abs(x-ztop))
tmp = np.where(np.around(depth3[:],decimals=1)==tmp1)[0]
ktop3 = tmp[0]

tmp1 = min(depth3[:], key=lambda x:abs(x-zbot))
tmp = np.where(np.around(depth3[:],decimals=1)==tmp1)[0]
kbot3 = tmp[0]

tmp1 = min(depth[:], key=lambda x:abs(x-ztop))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
ktop = tmp[0]

tmp1 = min(depth[:], key=lambda x:abs(x-zbot))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
kbot = tmp[0]

tmp1 = min(depth2[:], key=lambda x:abs(x-ztop))
tmp = np.where(np.around(depth2[:],decimals=1)==tmp1)[0]
ktop2 = tmp[0]

tmp1 = min(depth2[:], key=lambda x:abs(x-zbot))
tmp = np.where(np.around(depth2[:],decimals=1)==tmp1)[0]
kbot2 = tmp[0]

tmp1 = min(depth[:], key=lambda x:abs(x-2000.))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
k2k = tmp[0]

# Now find maximum of DWBC in 2000 m (at k2k) for each meridional slice

maxpos = np.zeros((vke.shape[1]))
for j in range(vke.shape[1]):
		tmp1 = vke[k2k,j,:].min()
		tmp2 = np.where(vke[k2k,j,:]==tmp1)[0]
		maxpos[j] = tmp2[0]
	
# Now shift every merdional slice to a reference level, defined by j = 0

ref_i = np.int64(maxpos[0])
sum_vke = np.zeros((vke.shape[0],vke.shape[2]))
sum2_vke = np.zeros((vke.shape[0],vke.shape[2]))

sum_urho = np.zeros((vke.shape[0],vke.shape[2]))
sum2_urho = np.zeros((vke.shape[0],vke.shape[2]))

sum_rho = np.zeros((vke.shape[0],vke.shape[2]))
sum2_rho = np.zeros((vke.shape[0],vke.shape[2]))

sum_divUrho = np.zeros((vke.shape[0],vke.shape[2]))
sum2_divUrho = np.zeros((vke.shape[0],vke.shape[2]))

sum_U_grad_rho = np.zeros((vke.shape[0],vke.shape[2]))
sum2_U_grad_rho = np.zeros((vke.shape[0],vke.shape[2]))


	
for i in range(200,400):
	for k in range(ktop,kbot):
		count1 = 0
		count2 = 0
		count3 = 0
		count4 = 0
		count5 = 0
		for j in range(vke.shape[1]):
			diff_i = np.int64(maxpos[j] - ref_i)
			if vke.mask[k,j,i+diff_i] == False:
				count1 = count1+1
				sum_vke[k,i] = sum_vke[k,i] + vke[k,j,i+diff_i]
			if urho.mask[k,j,i+diff_i] == False:
				count2 = count2+1
				sum_urho[k,i] = sum_urho[k,i] + urho[k,j,i+diff_i]
			if rho.mask[k,j,i+diff_i] == False:
				count3 = count3+1
				sum_rho[k,i] = sum_rho[k,i] + rho[k,j,i+diff_i]
			if divUrho.mask[k,j,i+diff_i] == False:
				count4 = count4+1
				sum_divUrho[k,i] = sum_divUrho[k,i] + divUrho[k,j,i+diff_i]
			if U_grad_rho.mask[k,j,i+diff_i] == False:
				count5 = count5+1
				sum_U_grad_rho[k,i] = sum_U_grad_rho[k,i] + U_grad_rho[k,j,i+diff_i]
			sum2_rho[k,i] = sum_rho[k,i] / np.float64(count3)
			sum2_divUrho[k,i] = sum_divUrho[k,i] / np.float64(count4)
			sum2_U_grad_rho[k,i] = sum_U_grad_rho[k,i] / np.float64(count5)
			sum2_urho[k,i] = sum_urho[k,i] / np.float64(count2)
			sum2_vke[k,i] = sum_vke[k,i] / np.float64(count1)
mean_rho = ma.masked_array(data=sum2_rho, mask = rho.mask[:,0,:])
mean_vke = ma.masked_array(data=sum2_vke, mask = vke.mask[:,0,:])
mean_urho = ma.masked_array(data=sum2_urho, mask = urho.mask[:,0,:])
mean_divUrho = ma.masked_array(data=sum2_divUrho, mask = rho.mask[:,0,:])
mean_U_grad_rho = ma.masked_array(data=sum2_U_grad_rho, mask = rho.mask[:,0,:])


# Now construct the cartesian grid
rlon=240
llon=210
x,z = np.meshgrid(lon[0,llon:rlon],depth[ktop:kbot])
x3,z3 = np.meshgrid(lon3[0,llon:rlon],depth3[ktop:kbot])
x2,z2 = np.meshgrid(lon2[0,llon:rlon],depth3[ktop:kbot])
data_vke =mean_vke[ktop:kbot,llon:rlon]
data_vke = ma.masked_invalid(data_vke,copy=False)
data_urho =mean_urho[ktop:kbot,llon:rlon]
data_urho = ma.masked_invalid(data_urho,copy=False)
data_rho =mean_rho[ktop:kbot,llon:rlon]
data_rho = ma.masked_invalid(data_rho,copy=False)
data_divUrho =mean_divUrho[ktop:kbot,llon:rlon]
data_divUrho = ma.masked_invalid(data_divUrho,copy=False)
data_U_grad_rho =mean_U_grad_rho[ktop:kbot,llon:rlon]
data_U_grad_rho = ma.masked_invalid(data_U_grad_rho,copy=False)
# Extract data U want to plot
#data_vrho 	= vrho[ktop3:kbot3,ilat3,iw:ie].copy()

high=1e-8
b = aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_divUrho,high)
b.canvas.set_window_title('Eddy Fluxes of Density between 11N and 18N')

a = aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_U_grad_rho,high)
a.canvas.set_window_title('Mean advection of density between 11N and 18N')
