"""
 histoSPOProtonExit.py  -  description
 ---------------------------------------------------------------------------------
 analysis of the soft proton distribution at the pore exit
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python histoSPOProtonExit N_each N_tot theta_0 angle_dis energy_0 model
 ---------------------------------------------------------------------------------
 Parameters:
 - N_each = number of simulated particles for each run
 - N_tot = total number of simulated particles
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
N_each = int(arg_list[1])
N_tot = int(arg_list[2])
theta_0 = float(arg_list[3])
angle_dis = int(arg_list[4])
if (angle_dis == 0): disname = "planar"
if (angle_dis == 1): disname = "cos"
energy_0 = float(arg_list[5])

model = arg_list[6]

N_runs = N_tot/N_each

event_out = []
energy_out = []
mdx_out = []
mdy_out = []
mdz_out = []
theta_out = []
phi_out = []
N_out = []
N_in_real = []
    
for jrun in xrange(N_runs):
    N_run = N_each*(jrun + 1)
	# Reading simulation output
    name_file_in = "./simulation_output/OUTPUT_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_run)+"_"+model+".dat"
    f_in = open(name_file_in, 'r')

    for line in f_in:
	    line = line.strip()
	    if not line.startswith("#"):
		    columns = line.split()
		    # angle_x err_angle_x Eff err_Eff N_out N_in solid_angle
		    columns[0] = int(columns[0])
		    columns[1] = float(columns[1])
		    columns[2] = float(columns[2])
		    columns[3] = float(columns[3])
		    columns[4] = float(columns[4])
		    columns[5] = float(columns[5])
		    columns[6] = float(columns[6])
		    columns[7] = int(columns[7])
		    columns[8] = int(columns[8])
	
		    event_out.append(columns[0])
		    energy_out.append(columns[1])
		    mdx_out.append(columns[2])
		    mdy_out.append(columns[3])
		    mdz_out.append(columns[4])
		    theta_out.append(columns[5])
		    phi_out.append(columns[6])
		    N_out.append(columns[7])
		    N_in_real.append(columns[8])

N_out = N_out[0]
N_in_real = N_in_real[0]

# Plot the results
plt.figure(1,figsize=[12,7])
gs1 = gridspec.GridSpec(1, 2, height_ratios=[4,1]) 
gs1.update(wspace=0.2)
ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])


# Theta histogram
vecThetaOut_max = np.max(theta_out)
vecThetaOut_min = np.min(theta_out)
angle_max = 30.

theta_bin = 1
n_bins_theta = int((vecThetaOut_max - vecThetaOut_min)/theta_bin)
#n_bins_theta = int((angle_max - vecThetaOut_min)/theta_bin)
if n_bins_theta < 1:
	n_bins_theta = 50

N_array_theta, bin_array_theta = np.histogram(theta_out, bins = n_bins_theta)
N_array_theta_norm = np.zeros(len(N_array_theta))
err_N_array_theta = np.zeros(len(N_array_theta))
err_angle_theta = np.zeros(len(N_array_theta))
angle_array_theta = np.zeros(len(N_array_theta))
left_angle_array_theta = np.zeros(len(N_array_theta))

for jn in xrange(len(N_array_theta)):
	err_angle_theta[jn] = (bin_array_theta[jn+1] - bin_array_theta[jn])/2.
	N_array_theta_norm[jn] = (float(N_array_theta[jn])/np.float(N_in_real))
	angle_array_theta[jn] = bin_array_theta[jn] + err_angle_theta[jn]
	left_angle_array_theta[jn] = bin_array_theta[jn]
	err_N_array_theta[jn] = (np.sqrt(float(N_array_theta[jn]))/np.float(N_in_real))


# out angle distribution
ax1.bar(left_angle_array_theta, N_array_theta_norm, width=2.*err_angle_theta, edgecolor="lightblue", facecolor='lightblue', lw = 0.01, alpha=0.7)
ax1.errorbar(angle_array_theta, N_array_theta_norm, xerr=err_angle_theta, yerr=err_N_array_theta, capsize=0, fmt='o', lw = 1, ecolor='black')

# Write output to file - theta
name_fileout = "./results/HISTOTHETA_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_tot)+"_"+model+".dat"
print "Writing in "+name_fileout
f_out = open(name_fileout, 'wb')
for iel in xrange(len(angle_array_theta)):
	# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle 
	f_out.write(str(angle_array_theta[iel])+" "+str(err_angle_theta[iel])+" "+str(N_array_theta_norm[iel])+" "+str(err_N_array_theta[iel])+" "+str(N_array_theta[iel])+"\n")

f_out.close()

ax1.set_yscale("log")
ax1.set_xlabel("Theta [deg.]")
ax1.set_ylabel("#detected particles/ #generated particles")
ax1.legend(numpoints=1, loc=1)


# Phi histogram
vecPhiOut_max = np.max(phi_out)
vecPhiOut_min = np.min(phi_out)
angle_max = 30.

phi_bin = 10
#n_bins_phi = int(vecphiOut_max - vecphiOut_min)/phi_bin
n_bins_phi = int((vecPhiOut_max - vecPhiOut_min)/phi_bin)
if n_bins_phi < 1:
	n_bins_phi = 50

N_array_phi, bin_array_phi = np.histogram(phi_out, bins = n_bins_phi)
N_array_phi_norm = np.zeros(len(N_array_phi))
err_N_array_phi = np.zeros(len(N_array_phi))
err_angle_phi = np.zeros(len(N_array_phi))
angle_array_phi = np.zeros(len(N_array_phi))
left_angle_array_phi = np.zeros(len(N_array_phi))

for jn in xrange(len(N_array_phi)):
	err_angle_phi[jn] = (bin_array_phi[jn+1] - bin_array_phi[jn])/2.
	N_array_phi_norm[jn] = (float(N_array_phi[jn])/np.float(N_in_real))
	angle_array_phi[jn] = bin_array_phi[jn] + err_angle_phi[jn]
	left_angle_array_phi[jn] = bin_array_phi[jn]
	err_N_array_phi[jn] = (np.sqrt(float(N_array_phi[jn]))/np.float(N_in_real))


# out angle distribution
ax2.bar(left_angle_array_phi, N_array_phi_norm, width=2.*err_angle_phi, edgecolor="yellow", facecolor='yellow', lw = 0.01, alpha=0.7)
ax2.errorbar(angle_array_phi, N_array_phi_norm, xerr=err_angle_phi, yerr=err_N_array_phi, capsize=0, fmt='o', lw = 1, ecolor='black')

# Write output to file - phi
name_fileout = "./results/HISTOPHI_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_tot)+"_"+model+".dat"
print "Writing in "+name_fileout
f_out = open(name_fileout, 'wb')
for iel in xrange(len(angle_array_phi)):
	# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle 
	f_out.write(str(angle_array_phi[iel])+" "+str(err_angle_phi[iel])+" "+str(N_array_phi_norm[iel])+" "+str(err_N_array_phi[iel])+" "+str(N_array_phi[iel])+"\n")

f_out.close()

ax2.set_yscale("log")
ax2.set_xlabel("Phi [deg.]")
#ax2.set_ylabel("# detected particles/ #generated particles")
ax2.legend(numpoints=1, loc=1)

title = str(N_in_real)+" events, "+str(energy_0)+" keV, "+str(theta_0)+" deg., "+model+" model"
plt.suptitle(title)

f_in.close()
plt.show()

