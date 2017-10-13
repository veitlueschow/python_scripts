from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import transport_meridional_slice as merslice
class Args: pass
inn = Args()
pas = Args()

def myround(x, prec=2, base=1.):
  return round(base * round(float(x)/base),prec)

filename="/work/mh0256/m300522/meta_storm/clines_30N35S.nc"
fh_clines	= Dataset(filename,mode='r')
filename='box_vke.nc'
fh_vke 		= Dataset(filename,mode='r')
filename='box_rbek.nc'
fh_rbek 	= Dataset(filename,mode='r')
filename='box_dlxv.nc'
fh_dlxv 	= Dataset(filename,mode='r')
filename='box_ddue.nc'
fh_ddue		= Dataset(filename,mode='r')

vke_ 	= fh_vke.variables["vke"]
lat_ 	= fh_vke.variables["lat_3"]
rbek_	= fh_rbek.variables["rbek"]
depth_	= fh_vke.variables["depth_4"]
dlxv_ 	= fh_dlxv.variables["dlxv"]
ddue_	= fh_ddue.variables["ddue"]
clines_ = fh_clines.variables["coastlines"]
clines_y_ = fh_clines.variables["lat"]
z_		= fh_clines.variables["depth"]


#pas.wmo	= wmo_[0,:,:,:].copy()
#pas.lat	= lat_[:,:].copy()
#pas.lon	= lon_[:,:].copy()

inn.vke	= vke_[0,:,:,:].copy()
inn.lat	= lat_[:,:].copy()
rbek	= rbek_[0,0,:,:].copy()
depth 	= depth_[:].copy()
dlxv 	= dlxv_[0,0,:,:].copy()
ddue 	= ddue_[0,:,:,:].copy()
clines	= clines_[:,:].copy()
clines_y= clines_y_[0,:].copy()
z		= z_[:,0].copy()

# Create coordinates

tmp = np.linspace(-35,30,num=66)
y,z = np.meshgrid(tmp,z)

full_tr = np.zeros((60,66))
dwbc_tr = np.zeros((60,66))

jbrei = 6
zbrei_inv = float(1) / float(2*jbrei+1)

for k in range(full_tr.shape[0]):
	print(k)
	ik=k+20
	cline_x = np.int64(clines[k,:])
	for i in range(inn.vke.shape[2]):
		for j in range(inn.vke.shape[1]):
			zlat=dlxv[j,i] / (np.float64(2*jbrei*111111))
			if inn.vke.mask[ik,j,i] == False:
				if rbek[j,i] == 4. or rbek[j,i] == 5. or rbek[j,i] == 3. or rbek[j,i] == 2. or rbek[j,i] == 1.:
					for l in range(-jbrei,jbrei+1):
						tmp = inn.lat[j,i]+l*zlat+35.
						tmp = np.float64(tmp)
						lbrei = myround(tmp)
						if lbrei >= 0 and lbrei <= 65:
							lbrei = np.int64(lbrei)
							full_tr[k,lbrei] = full_tr[k,lbrei] + inn.vke[ik,j,i]*dlxv[j,i]*ddue[ik,j,i]*zbrei_inv # Output in Sverdrup [m^3 s^-1]
							if i >= (cline_x[j] -5) and i <= (cline_x[j]+20):
								dwbc_tr[k,lbrei] = dwbc_tr[k,lbrei] + inn.vke[ik,j,i]*dlxv[j,i]*ddue[ik,j,i]*zbrei_inv

#write.write_depth_meridional("transport_meridional_full3.nc",full_tr,'transport_mer_full',y,z)
#write.write_depth_meridional("transport_meridional_dwbc3.nc",dwbc_tr,'transport_mer_dwbc',y,z)

plt.figure()
plt.contourf(y,z,full_tr)
plt.colorbar()
plt.show()
