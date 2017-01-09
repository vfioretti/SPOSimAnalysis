"""
 histoSPOProtonExit.py  -  description
 ---------------------------------------------------------------------------------
 analysis of the soft proton distribution at the pore exit
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python histoSPOProtonExit N_in theta_0 angle_dis energy_0 model
 ---------------------------------------------------------------------------------
 Parameters:
 - N_in = number of simulated particles
 - theta_0 = incoming angle
 - energy_0 = incoming energy
 - angle_bin = bin angle
 - model = scattering model
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2016/06/01: creation date
"""



import pyfits
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Import the input parameters
arg_list = sys.argv
N_in = int(arg_list[1])
theta_0 = float(arg_list[2])
angle_dis = int(arg_list[3])
if (angle_dis == 0): disname = "planar"
if (angle_dis == 1): disname = "cos"
energy_0 = float(arg_list[4])

model = arg_list[5]

# Reading simulation output
name_file_in = "./simulation_output/OUTPUT_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_in)+"_"+model+".dat"
f_in = open(name_file_in, 'r')
	
event_out = []
energy_out = []
mdx_out = []
mdy_out = []
mdz_out = []
theta_out = []
phi_out = []

for line in f_in:
	line = line.strip()
	if not line.startswith("#"):
		columns = line.split()
		# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle
		columns[0] = float(columns[0])
		columns[1] = float(columns[1])
		columns[2] = float(columns[2])
		columns[3] = float(columns[3])
		columns[4] = float(columns[4])
		columns[5] = float(columns[5])
		columns[6] = float(columns[6])
	
		event_out.append(columns[0])
		energy_out.append(columns[1])
		mdx_out.append(columns[2])
		mdy_out.append(columns[3])
		mdz_out.append(columns[4])
		theta_out.append(columns[5])
		phi_out.append(columns[6])



# Plot the results
gs = gridspec.GridSpec(1, 1, height_ratios=[4,1]) 
gs.update(hspace=0.0)
ax = plt.subplot(gs[0])


# Theta histogram
vecThetaOut_max = np.max(theta_out)
vecThetaOut_min = np.min(theta_out)
angle_max = 100.

theta_bin = 1
#n_bins_theta = int(vecThetaOut_max - vecThetaOut_min)/theta_bin
n_bins_theta = int((angle_max - vecThetaOut_min)/theta_bin)
if n_bins_theta < 1:
	n_bins_theta = 50

N_array_theta, bin_array_theta = np.histogram(theta_out, bins = n_bins_theta)
N_array_theta_norm = np.zeros(len(N_array_theta))
err_N_array_theta = np.zeros(len(N_array_theta))
err_angle = np.zeros(len(N_array_theta))
angle_array = np.zeros(len(N_array_theta))
left_angle_array = np.zeros(len(N_array_theta))

for jn in xrange(len(N_array_theta)):
	err_angle[jn] = (bin_array_theta[jn+1] - bin_array_theta[jn])/2.
	N_array_theta_norm[jn] = (float(N_array_theta[jn])/np.float(N_in))
	angle_array[jn] = bin_array_theta[jn] + err_angle[jn]
	left_angle_array[jn] = bin_array_theta[jn]
	err_N_array_theta[jn] = (np.sqrt(float(N_array_theta[jn]))/np.float(N_in))


# out angle distribution
ax.bar(left_angle_array, N_array_theta_norm, width=2.*err_angle, edgecolor="lightblue", facecolor='lightblue', lw = 0.01, label='BoGEMMS', alpha=0.7)
ax.errorbar(angle_array, N_array_theta_norm, xerr=err_angle, yerr=err_N_array_theta, capsize=0, fmt='none', lw = 1, ecolor='black')


# Write output to file
name_fileout = "./results/HISTO_"+str(int(energy_0))+"keV_"+str(theta_0)+"deg_"+str(N_in)+"_"+model+".dat"
print "Writing in "+name_fileout
f_out = open(name_fileout, 'wb')
for iel in xrange(len(angle_array)):
	# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle 
	f_out.write(str(angle_array[iel])+" "+str(err_angle[iel])+" "+str(N_array_theta_norm[iel])+" "+str(err_N_array_theta[iel])+" "+str(N_array_theta[iel])+"\n")


ax.set_yscale("log")
ax.set_xlabel("Theta [deg.]")
ax.set_ylabel("# detected particles/ #generated particles")
ax.legend(numpoints=1, loc=1)
title = str(N_in)+" events, "+str(energy_0)+" keV, "+str(theta_0)+" deg., "+model+" model"
ax.set_title(title)

f_out.close()
plt.show()
hdulist.close()
