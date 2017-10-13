from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
class Args: pass
read = Args()
pas = Args()




filename='../atlantic_moc.nc'
fh_amoc 		= Dataset(filename,mode='r')
amoc_ 	= fh_amoc.variables["atlantic_moc"]

depth_	= fh_amoc.variables["depth_2"]
depth 	= depth_[:].copy()
tmp = np.linspace(-89,90,num=180)
y1,z1 = np.meshgrid(tmp,depth)

read.amoc	= amoc_[0,:,:,0].copy()
f3 = plt.figure()
plt.set_cmap('coolwarm')
v = np.linspace(-3e10, 3e10, 100, endpoint=True)
plt.contourf(y1,-z1,read.amoc[:,:],v,extend="both")
plt.colorbar()
