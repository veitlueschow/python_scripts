from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
import numpy.ma as ma
import xz_tools as xz_tools
class Args: pass

  
# for "box" use : cdo selindexbox,1950,2600,510,580 ( if you start from "_box" files)
# for use from whole: cdo selindexbox,1950,2600,1010,1080:

# cdo selindexbox,1950,2600,1010,1080 ../div_U+rho+.nc div_U+rho+_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../v+rho+.nc v+rho+_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../w+rho+.nc w+rho+_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../u+rho+.nc u+rho+_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../div_U+rho+.nc div_U+rho+_stripe.nc
#cdo selindexbox,1950,2600,1010,1080 ../U_grad_rho.nc U_grad_rho_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 /work/mh0256/m300522/data_storm/tape/1960s/vke_tm.nc vke_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 /work/mh0256/m300522/data_storm/tape/1960s/rhopoto_tm.nc rhopoto_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../dz_rhopoto.nc dz_rhopoto_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../dx_rhopoto.nc dx_rhopoto_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../dy_rhopoto.nc dy_rhopoto_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../Pe2Pm_hor.nc Pe2Pm_hor_stripe.nc
# cdo selindexbox,1950,2600,1010,1080 ../kappa/kappa_nakamura.nc kappan_stripe.nc

#--------------------------Read Stuff---------------------------------------------------

vrho = tools.netread_data('v+rho+_stripe.nc','vrho_eddy') # vrho_eddy
lat3,lon3,depth3 = tools.netread_grid('v+rho+_stripe.nc','lat_3','lon_3','depth_2')
urho = tools.netread_data('u+rho+_stripe.nc','urho_eddy') # urho_eddy
lat2,lon2,depth2 = tools.netread_grid('u+rho+_stripe.nc','lat_2','lon_2','depth_2')
wrho = tools.netread_data('w+rho+_stripe.nc','wrho_eddy') # wrho_eddy

vke = tools.netread_data('vke_stripe.nc','vke') # vke
uko = tools.netread_data('uko_stripe.nc','uko') # uko
wo_p = tools.netread_data('wo_p_stripe.nc','wo') # wo at p-point

rho = tools.netread_data('rhopoto_stripe.nc','rhopoto') # density
lat,lon,depth = tools.netread_grid('rhopoto_stripe.nc','lat','lon','depth_2')

divUrho = tools.netread_data("div_U+rho+_stripe.nc","div_Urho_eddy") # eddy flux divergence
U_grad_rho = tools.netread_data("U_grad_rho_stripe.nc","mean_advection_rho") # mean density advection

dx_rho = tools.netread_data('dx_rhopoto_stripe.nc','dx_rhopoto') # zonal density gradient
dy_rho = tools.netread_data('dy_rhopoto_stripe.nc','dy_rhopoto') # meridional density gradient
dz_rho = tools.netread_data('dz_rhopoto_stripe.nc','dz_rhopoto') # vertical density gradient

Pe2Pm_hor = tools.netread_data('Pe2Pm_hor_stripe.nc','urho_eddy') # Potential energy conversion from eddy to mean

kappan = tools.netread_data('kappan_stripe.nc','kappa_nakamura') # kappa from Nakamura 2000

###--Choose the depth levels that you are interested in and get indices--
ztop 	= 1200.
zbot	= 3000.

tmp1 = min(depth[:], key=lambda x:abs(x-ztop))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
ktop = tmp[0]
tmp1 = min(depth[:], key=lambda x:abs(x-zbot))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
kbot = tmp[0]
tmp1 = min(depth[:], key=lambda x:abs(x-2000.))
tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
k2k = tmp[0]

###----find maximum of DWBC in 2000 m (at k2k) for each meridional slice--

maxpos = np.zeros((vke.shape[1]))
for j in range(vke.shape[1]):
		tmp1 = vke[k2k,j,:].min()
		tmp2 = np.where(vke[k2k,j,:]==tmp1)[0]
		maxpos[j] = tmp2[0]



###-------------compute meridional mean----------------------------
lx=400 # for average
rx=600 # for average
llon=470 # for plot
rlon=510 # for plot

mean_rho, data_rho = xz_tools.meri_average(vke,maxpos,rho,ktop,kbot,lx,rx,llon,rlon)
mean_vke, data_vke = xz_tools.meri_average(vke,maxpos,vke,ktop,kbot,lx,rx,llon,rlon)
#mean_urho, data_urho = xz_tools.meri_average(vke,maxpos,urho,ktop,kbot,lx,rx,llon,rlon)
mean_divUrho, data_divUrho = xz_tools.meri_average(vke,maxpos,divUrho,ktop,kbot,lx,rx,llon,rlon)
mean_U_grad_rho, data_U_grad_rho = xz_tools.meri_average(vke,maxpos,U_grad_rho,ktop,kbot,lx,rx,llon,rlon)
mean_wrho, data_wrho = xz_tools.meri_average(vke,maxpos,wrho,ktop,kbot,lx,rx,llon,rlon)
mean_uko, data_uko = xz_tools.meri_average(vke,maxpos,uko,ktop,kbot,lx,rx,llon,rlon)
mean_wo_p, data_wo_p = xz_tools.meri_average(vke,maxpos,wo_p,ktop,kbot,lx,rx,llon,rlon)
mean_dx_rho, data_dx_rho = xz_tools.meri_average(vke,maxpos,dx_rho,ktop,kbot,lx,rx,llon,rlon)
mean_dz_rho, data_dz_rho = xz_tools.meri_average(vke,maxpos,dz_rho,ktop,kbot,lx,rx,llon,rlon)
mean_dy_rho, data_dy_rho = xz_tools.meri_average(vke,maxpos,dy_rho,ktop,kbot,lx,rx,llon,rlon)
mean_Pe2Pm_hor, data_Pe2Pm_hor = xz_tools.meri_average(vke,maxpos,Pe2Pm_hor,ktop,kbot,lx,rx,llon,rlon)
mean_kappan, data_kappan = xz_tools.meri_average(vke,maxpos,kappan,ktop,kbot,lx,rx,llon,rlon)

###-------------construct cartesian grid----------------------------
x,z = np.meshgrid(lon[0,llon:rlon],depth[ktop:kbot])
x3,z3 = np.meshgrid(lon3[0,llon:rlon],depth3[ktop:kbot])
x2,z2 = np.meshgrid(lon2[0,llon:rlon],depth3[ktop:kbot])



###-----------get isopycnal normal vector ----------------------------

nrho_x = data_dx_rho[:,:].copy()
nrho_z = data_dz_rho[:,:].copy()
abs_n = np.double(np.sqrt(nrho_x*nrho_x+nrho_z*nrho_z))

nrho_x = np.double(nrho_x / abs_n)
nrho_z = np.double(nrho_z / abs_n)

###-------------get diapycnal velocties from that---------------------
w_dia = data_divUrho * nrho_z
u_dia = data_divUrho * nrho_x


### ---------------plot...---------------------------------------------
high=2e-8
b = aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_divUrho,high,'Eddy Density Fluxes between 10S and 18S')
b.canvas.set_window_title('Eddy Fluxes of Density between 10S and 18S')
b.savefig("10S-18S_eddy_flux_divergence.png")


#a = aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_U_grad_rho,high,'Mean density advection between 10S and 18S')
#a.canvas.set_window_title('Mean advection of density between 10S and 16S')
#a.savefig("10S-18S_density_advection.png")

c = aplot.myplot_quiver(x3,z3,data_rho,x,z,220*u_dia,w_dia,'Diapycnal Eddy fluxes')
c.savefig("10S-18S_diapycnal_fluxes.png")

high = 1e-11
d = aplot.myplot_multi2(x,z,data_rho,x,z,data_vke,x,z,data_Pe2Pm_hor,high,'Conversion eddy to mean potential energy 10S and 18S')
d.canvas.set_window_title('Conversion eddy to mean potential energy 10S and 16S')
d.savefig("10S-18S_Conv_eddy2mean.png")
