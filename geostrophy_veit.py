from netCDF4 import Dataset
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from imp import reload
class Args: pass
inn = Args()
v = Args()
u = Args()

filename='box_po.nc'
fh 		= Dataset(filename,mode='r')
po_ 	= fh.variables["po"]
inn.po	= po_[0,0,:,:].copy()

filename='box_zo.nc'
fh 		= Dataset(filename,mode='r')
zo_ 	= fh.variables["zo"]
inn.zo	= zo_[0,0,:,:].copy()
lat_	= fh.variables["lat"]
lat	= lat_[:,:].copy()

filename='box_dlxp.nc'
fh 		= Dataset(filename,mode='r')
dlxp_ 	= fh.variables["dlxp"]
dlxp	= dlxp_[0,0,:,:].copy()

filename='box_dlyp.nc'
fh 		= Dataset(filename,mode='r')
dlyp_ 	= fh.variables["dlyp"]
dlyp	= dlyp_[0,0,:,:].copy()

filename='box_vke.nc'
fh 		= Dataset(filename,mode='r')
vke_ 	= fh.variables["vke"]
inn.vke	= vke_[0,0,:,:].copy()
lat3_	= fh.variables["lat_3"]
lat_3	= lat3_[:,:].copy()
lon3_ 	= fh.variables["lon_3"]
lon_3 	= lon3_[:,:].copy()

filename='box_uko.nc'
fh 		= Dataset(filename,mode='r')
uko_ 	= fh.variables["uko"]
inn.uko	= uko_[0,0,:,:].copy()
lat2_	= fh.variables["lat_2"]
lat_2	= lat2_[:,:].copy()
lon2_ 	= fh.variables["lon_2"]
lon_2 	= lon2_[:,:].copy()


inn.pzo = 9.81*inn.zo
ftwov = 0.0001458425*np.sin(lat*0.017453)
fdy=1./(dlyp*ftwov)
ftwou = 0.0001458425*np.sin(lat*0.017453)
fdx=1./(dlxp*ftwou)
w1 = np.zeros((inn.po.shape))
w2 = np.zeros((inn.po.shape))
w3 = np.zeros((inn.po.shape))
w4 = np.zeros((inn.po.shape))
v1 = np.zeros((inn.po.shape))
v2 = np.zeros((inn.po.shape))
v3 = np.zeros((inn.po.shape))
v4 = np.zeros((inn.po.shape))
l1 = np.zeros((inn.po.shape))
l2 = np.zeros((inn.po.shape))
l3 = np.zeros((inn.po.shape))
l4 = np.zeros((inn.po.shape))
m1 = np.zeros((inn.po.shape))
m2 = np.zeros((inn.po.shape))
m3 = np.zeros((inn.po.shape))
m4 = np.zeros((inn.po.shape))

#First deal with meridional velocity, therefore zonal pressure gradient. Seperate between contribution from pressure and elevation
for i in range(1,inn.po.shape[1]-2):
	for j in range(1,inn.po.shape[0]-2):
		if inn.po.mask[j,i] == False and inn.po.mask[j+1,i] == False \
		and inn.po.mask[j,i-1] == False and inn.po.mask[j+1,i-1] == False \
		and inn.po.mask[j,i+1] == False and inn.po.mask[j+1,i+1] == False:
			w1[j,i] = (inn.po[j,i+1]-inn.po[j,i])*fdx[j,i]
			w2[j,i] = (inn.po[j,i]-inn.po[j,i-1])*fdx[j,i]
			w3[j,i] = (inn.po[j+1,i+1]-inn.po[j+1,i])*fdx[j+1,i]
			w4[j,i] = (inn.po[j+1,i]-inn.po[j+1,i-1])*fdx[j+1,i]
			v1[j,i] = (inn.pzo[j,i+1]-inn.pzo[j,i])*fdx[j,i]
			v2[j,i] = (inn.pzo[j,i]-inn.pzo[j,i-1])*fdx[j,i]
			v3[j,i] = (inn.pzo[j+1,i+1]-inn.pzo[j+1,i])*fdx[j+1,i]
			v4[j,i] = (inn.pzo[j+1,i]-inn.pzo[j+1,i-1])*fdx[j+1,i]
		if inn.po.mask[j,i] == False and inn.po.mask[j,i+1] == False \
		and inn.po.mask[j-1,i] == False and inn.po.mask[j-1,i+1] == False \
		and inn.po.mask[j+1,i] == False and inn.po.mask[j+1,i+1] == False:
			l1[j,i] = (inn.po[j+1,i]-inn.po[j,i])*fdy[j,i]
			l1[j,i] = (inn.po[j,i]-inn.po[j-1,i])*fdy[j,i]
			l3[j,i] = (inn.po[j+1,i-1]-inn.po[j,i-1])*fdy[j,i-1]
			l4[j,i] = (inn.po[j,i-1]-inn.po[j-1,i-1])*fdy[j,i-1]
			m1[j,i] = (inn.pzo[j+1,i]-inn.pzo[j,i])*fdy[j,i]
			m1[j,i] = (inn.pzo[j,i]-inn.pzo[j-1,i])*fdy[j,i]
			m3[j,i] = (inn.pzo[j+1,i-1]-inn.pzo[j,i-1])*fdy[j,i-1]
			m4[j,i] = (inn.pzo[j,i-1]-inn.pzo[j-1,i-1])*fdy[j,i-1]

v.contri_po 	= 0.25*(w1+w2+w3+w4)
v.contri_pzo 	= 0.25*(v1+v2+v3+v4)
v.vg 			= v.contri_po + v.contri_pzo 
v.contri_po 	= np.ma.array(v.contri_po, mask = inn.vke.mask)
v.contri_pzo 	= np.ma.array(v.contri_pzo, mask = inn.vke.mask)
v.vg 			= np.ma.array(v.vg, mask = inn.vke.mask)

u.contri_po 	= 0.25*(l1+l2+l3+l4)
u.contri_pzo 	= 0.25*(m1+m2+m3+m4)
u.ug 			= u.contri_po + u.contri_pzo 
u.contri_po 	= np.ma.array(u.contri_po, mask = inn.uko.mask)
u.contri_pzo 	= np.ma.array(u.contri_pzo, mask = inn.uko.mask)
u.ug 			= np.ma.array(u.ug, mask = inn.uko.mask)



# Now write netCDF file for v
f1 = Dataset('box_vg_py.nc','w',clobber=True, format='NETCDF4') #'w' stands for write
xs 	= f1.createDimension('x', lat_3.shape[1])
ys 	= f1.createDimension('y', lat_3.shape[0])
lons 	= f1.createVariable('lon_3', np.float64,('y','x'))
lats 	= f1.createVariable('lat_3', np.float64,('y','x'))
data_v1		= f1.createVariable('vg',np.float64,('y','x'))
data_v2		= f1.createVariable('contri_po',np.float64,('y','x'))
data_v3		= f1.createVariable('contri_pzo',np.float64,('y','x'))

lats[:,:]	= lat_3.copy()
lons[:,]	= lon_3.copy()
data_v1[:,:]	= v.vg.copy()
data_v2[:,:]	= v.contri_po.copy()
data_v3[:,:]	= v.contri_pzo.copy()
#lons[:] = np.arange(48, 51.1, 0.1)

f1.Conventions = 'CF-1.0'

f1.close()

# Now for u
f2 = Dataset('box_ug_py.nc','w',clobber=True, format='NETCDF4') #'w' stands for write
xs 	= f2.createDimension('x', lat_2.shape[1])
ys 	= f2.createDimension('y', lat_2.shape[0])
lons 	= f2.createVariable('lon_2', np.float64,('y','x'))
lats 	= f2.createVariable('lat_2', np.float64,('y','x'))
data_u1		= f2.createVariable('ug',np.float64,('y','x'))
data_u2		= f2.createVariable('contri_po',np.float64,('y','x'))
data_u3		= f2.createVariable('contri_pzo',np.float64,('y','x'))

lats[:,:]	= lat_2.copy()
lons[:,]	= lon_2.copy()
data_u1[:,:]	= u.ug.copy()
data_u2[:,:]	= u.contri_po.copy()
data_u3[:,:]	= u.contri_pzo.copy()
#lons[:] = np.arange(48, 51.1, 0.1)

f2.Conventions = 'CF-1.0'

f2.close()
