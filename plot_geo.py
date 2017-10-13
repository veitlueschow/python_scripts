from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
class Args: pass
inn = Args()

fh1 = Dataset('outx',mode='r')
fh2 = Dataset('diffx',mode='r')
fh3 = Dataset('box_uko',mode='r')

ug_ 	= fh1.variables["po"]
diffu_ 	= fh2.variables["uko"]
uko_	= fh3.variables["uko"]


inn.ug		= ug_[0,:,:].copy()
inn.diffu	= diffu_[0,:,:].copy()
inn.uko		= uko_[0,:,:].copy()
