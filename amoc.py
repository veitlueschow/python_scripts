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

filename='box2_wmo.nc'
fh_wmo 		= Dataset(filename,mode='r')
filename='box2_alat.nc'
fh_alat 	= Dataset(filename,mode='r')
filename='box2_rbek.nc'
fh_rbek 	= Dataset(filename,mode='r')
filename='box2_dlxp.nc'
fh_dlxp 	= Dataset(filename,mode='r')
filename='box2_dlyp.nc'
fh_dlyp 	= Dataset(filename,mode='r')

wmo_ 	= fh_wmo.variables["wmo"]
lat_ 	= fh_wmo.variables["lat"]
lon_ 	= fh_wmo.variables["lon"]
rbek_	= fh_rbek.variables["rbek"]
depth_	= fh_wmo.variables["depth_15"]
dlyp_ 	= fh_dlyp.variables["dlyp"]
dlxp_ 	= fh_dlxp.variables["dlxp"]


#pas.wmo	= wmo_[0,:,:,:].copy()
#pas.lat	= lat_[:,:].copy()
#pas.lon	= lon_[:,:].copy()

inn.wmo	= wmo_[0,:,:,:].copy()
inn.lat	= lat_[:,:].copy()
inn.lon	= lon_[:,:].copy()
rbek	= rbek_[0,0,:,:].copy()
depth 	= depth_[:].copy()
dlyp 	= dlyp_[0,0,:,:].copy()
dlxp 	= dlxp_[0,0,:,:].copy()

# Create coordinates
#tmp = np.linspace(0,35,num=36)
tmp = np.linspace(1,180,num=180)


y,z = np.meshgrid(tmp,depth)

amoc = np.zeros((81,180))
#amoc_ = np.zeros((80,85))
amoc_count = np.zeros((81,180))

#amoc = np.zeros((81,36))
#amoc_ = np.zeros((81,36))
#amoc_count = np.zeros((81,36))


jbrei = 3
zbrei_inv = float(1) / float(2*jbrei+1)

for i in range(inn.wmo.shape[2]): # approx atlantic
	print(i)
#for i in range(2500,2501):
	for j in range(2,inn.wmo.shape[1]-2):
	#for j in range(inn.wmo.shape[1]):
		lbrei = myround(inn.lat[j,i])
		if dlyp.mask[j,i] == False:
			#lbrei = np.int64(lbrei*2 + 98)
			if rbek[j,i] == 4. or rbek[j,i] == 5. or rbek[j,i] == 3. or rbek[j,i] == 2. or rbek[j,i] == 1.:
				#print(min(dlxp[j,i],dlyp[j,i]))
				zlat=min(dlxp[j,i],dlyp[j,i]) / (np.float64(2*jbrei*111111))
				for k in range(amoc.shape[0]):
					for l in range(-jbrei,jbrei+1):
						tmp = inn.lat[j,i]+l*zlat+90.
						tmp = np.float64(tmp)
						lbrei = myround(tmp)
						if lbrei >= 1 and lbrei <= 180 and inn.wmo.mask[k,j,i] == False:
							lbrei = np.int64(lbrei - 1)
							amoc[k,lbrei] = amoc[k,lbrei] - inn.wmo[k,j,i] * zbrei_inv

write.write_depth_meridional("amoc_untouched.nc",amoc,'amoc',y,z)

amoc_old = amoc.copy()
amoc_untouched = amoc.copy()
for lb in range(amoc.shape[1]-2,-1,-1):
	if abs(np.mean(amoc[:,lb])) > 1.:
		amoc[:,lb]=amoc[:,lb+1] + amoc[:,lb]
amoc_down = amoc.copy()
amoc = amoc_old.copy()
for lb in range(1,amoc.shape[1]-1):
	if abs(np.mean(amoc[:,lb])) > 1.:
		amoc[:,lb]=amoc[:,lb-1] + amoc[:,lb]
amoc_up = -amoc.copy()
							
f1 = plt.figure()
plt.set_cmap('coolwarm')
v = np.linspace(-3e10, 3e10, 100, endpoint=True)
plt.contourf(y,-z,amoc,v,extend="both")
plt.colorbar()
