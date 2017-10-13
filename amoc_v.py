from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
class Args: pass
inn = Args()
pas = Args()

def myround(x, prec=2, base=1.):
  return round(base * round(float(x)/base),prec)

filename='box2_vke.nc'
fh_vke 		= Dataset(filename,mode='r')
filename='box2_alat.nc'
fh_alat 	= Dataset(filename,mode='r')
filename='box2_rbek.nc'
fh_rbek 	= Dataset(filename,mode='r')
filename='box2_dlxv.nc'
fh_dlxv 	= Dataset(filename,mode='r')
filename='box2_ddue.nc'
fh_ddue		= Dataset(filename,mode='r')

vke_ 	= fh_vke.variables["vke"]
lat_ 	= fh_vke.variables["lat_3"]
rbek_	= fh_rbek.variables["rbek"]
depth_	= fh_vke.variables["depth_4"]
dlxv_ 	= fh_dlxv.variables["dlxv"]
ddue_	= fh_ddue.variables["ddue"]


#pas.wmo	= wmo_[0,:,:,:].copy()
#pas.lat	= lat_[:,:].copy()
#pas.lon	= lon_[:,:].copy()

inn.vke	= vke_[0,:,:,:].copy()
inn.lat	= lat_[:,:].copy()
rbek	= rbek_[0,0,:,:].copy()
depth 	= depth_[:].copy()
dlxv 	= dlxv_[0,0,:,:].copy()
ddue 	= ddue_[0,:,:,:].copy()

# Create coordinates
#tmp = np.linspace(0,35,num=36)
tmp = np.linspace(-40,40,num=81)


y,z = np.meshgrid(tmp,depth)

amoc = np.zeros((80,81))
#amoc_ = np.zeros((80,85))
amoc_count = np.zeros((80,81))

#amoc = np.zeros((81,36))
#amoc_ = np.zeros((81,36))
#amoc_count = np.zeros((81,36))


jbrei = 3
zbrei_inv = float(1) / float(2*jbrei+1)

for i in range(inn.vke.shape[2]): # approx atlantic
	print(i)
	for j in range(400,1500):
		lbrei = myround(inn.lat[j,i])
		if dlxv.mask[j,i] == False:
			#lbrei = np.int64(lbrei*2 + 98)
			if rbek[j,i] == 4. or rbek[j,i] == 5. or rbek[j,i] == 3. or rbek[j,i] == 2. or rbek[j,i] == 1.:
				#print(min(dlxp[j,i],dlxv[j,i]))
				zlat=dlxv[j,i] / (np.float64(2*jbrei*111111))
				for k in range(amoc.shape[0]):
					for l in range(-jbrei,jbrei+1):
						tmp = inn.lat[j,i]+l*zlat+40.
						tmp = np.float64(tmp)
						lbrei = myround(tmp)
						if lbrei >= 0 and lbrei <= 80 and inn.vke.mask[k,j,i] == False:
							lbrei = np.int64(lbrei)
							area = dlxv[j,i]*ddue[k,j,i]
							amoc[k,lbrei] = amoc[k,lbrei] + inn.vke[k,j,i]*area*1025.0* zbrei_inv

write.write_depth_meridional("amoc_untouched_v2.nc",amoc,'amoc',y,z)

amoc_old = amoc.copy()
amoc_untouched = amoc.copy()
for kb in range(1,amoc.shape[0]-2):
	#if abs(np.mean(amoc[:,kb])) > 1.:
	amoc[kb,:]=amoc[kb-1,:] + amoc[kb,:]
amoc_down = amoc.copy()

amoc = amoc_old.copy()
for kb in range(amoc.shape[0]-2,-1,-1):
	#if abs(np.mean(amoc[:,kb])) > 1.:
	amoc[kb,:]=amoc[kb+1,:] + amoc[kb,:]
amoc_up = amoc.copy()
							
f1 = plt.figure()
plt.set_cmap('coolwarm')
v = np.linspace(-3e10, 3e10, 100, endpoint=True)
plt.contourf(y,-z,amoc,v,extend="both")
plt.colorbar()
