"""
 computeSPOProtonExit.py  -  description
 ---------------------------------------------------------------------------------
 computing the soft proton distribution at the pore exit
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python computeSPOProtonExit filedir N_file N_in theta_0 angle_dis energy_0 model
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - N_in = number of simulated particles
 - theta_0 = incoming angle
 - angle_dis = 0 [planar], 1 [cos]
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
filedir = arg_list[1]
N_fits = int(arg_list[2])
N_in = int(arg_list[3])
theta_0 = float(arg_list[4])
angle_dis = float(arg_list[5])
if (angle_dis == 0): disname = "planar"
if (angle_dis == 1): disname = "cos"
energy_0 = float(arg_list[6])

model = arg_list[7]

HPlate_vol_id_min = 0 
HPlate_vol_id_max = 3 
HRib_vol_id_min = 4 
HRib_vol_id_max = 7 

sphere_vol_id = 10000000


N_in_real = 0

vecEventIDOut = []
vecThetaOut = []
vecPhiOut = []
vecEnergyOut = []
vecMDXOut = []
vecMDYOut = []
vecMDZOut = []


for jfits in xrange(N_fits):
    
    print '%%%%%%%%%%%%%% READING BoGEMMS FILE: ', jfits+1
    hdulist = pyfits.open(filedir+'/xyz.'+str(jfits)+'.fits.gz')
    
    tbdata = hdulist[1].data
    cols = hdulist[1].columns
    
    evt_id = tbdata.field('EVT_ID')
    trk_id = tbdata.field('TRK_ID')
    vol_id = tbdata.field('VOLUME_ID')
    ene_ent = tbdata.field('E_KIN_ENT')
    ene_exit = tbdata.field('E_KIN_EXIT')
    MDZ_ent = tbdata.field('MDZ_ENT')
    MDY_ent = tbdata.field('MDY_ENT')
    MDX_ent = tbdata.field('MDX_ENT')
    Z_ent = tbdata.field('Z_ENT')
    Y_ent = tbdata.field('Y_ENT')
    X_ent = tbdata.field('X_ENT')
    MDZ_exit = tbdata.field('MDZ_EXIT')
    MDY_exit = tbdata.field('MDY_EXIT')
    MDX_exit = tbdata.field('MDX_EXIT')
    Z_exit = tbdata.field('Z_EXIT')
    Y_exit = tbdata.field('Y_EXIT')
    X_exit = tbdata.field('X_EXIT')
    part_id = tbdata.field('PARTICLE_ID')
    
    # ---------------------------------

    where_pore = np.where((((vol_id >= HPlate_vol_id_min) & (vol_id <= HPlate_vol_id_max)) | ((vol_id >= HRib_vol_id_min) & (vol_id <= HRib_vol_id_max))) & (part_id == 2212))
    evt_id_pore = evt_id[where_pore]
    vol_id_pore = vol_id[where_pore]
    ene_pore = ene_exit[where_pore]
    mdz_pore = MDZ_exit[where_pore]
    mdy_pore = MDY_exit[where_pore]
    mdx_pore = MDX_exit[where_pore]
        
    where_sphere = np.where((vol_id == sphere_vol_id) & (part_id == 2212) & (Z_ent < 0))
    evt_id_sphere = evt_id[where_sphere]
    ene_sphere = ene_ent[where_sphere]
    mdz_sphere = MDZ_ent[where_sphere]
    mdy_sphere = MDY_ent[where_sphere]
    mdx_sphere = MDX_ent[where_sphere]

    jev_pore = 0
    while (1):
		same_ev = np.where(evt_id_sphere == evt_id_pore[jev_pore])
		same_ev_pore = np.where(evt_id_pore == evt_id_pore[jev_pore])
		N_in_real+=1
		same_ev = same_ev[0]
		same_ev_pore = same_ev_pore[0]
		
		if (same_ev):
			ene_pore_extract = ene_pore[same_ev_pore]
			mdz_pore_extract = mdz_pore[same_ev_pore]
			mdy_pore_extract = mdy_pore[same_ev_pore]
			mdx_pore_extract = mdx_pore[same_ev_pore]
		
			ene_sphere_extract = ene_sphere[same_ev]
			mdz_sphere_extract = mdz_sphere[same_ev]
			mdy_sphere_extract = mdy_sphere[same_ev]
			mdx_sphere_extract = mdx_sphere[same_ev] 
			
			for jsame in xrange(len(ene_pore_extract)):
				if (ene_pore_extract[jsame] == ene_sphere_extract):
				    if (mdz_pore_extract[jsame] == mdz_sphere_extract):
				        if (mdy_pore_extract[jsame] == mdy_sphere_extract):
				            if (mdx_pore_extract[jsame] == mdx_sphere_extract):
				                 vecEventIDOut.append(evt_id_pore[jev_pore])
				                 vecEnergyOut.append(ene_sphere_extract)
				                 vecMDXOut.append(mdx_sphere_extract)
				                 vecMDYOut.append(mdy_sphere_extract)
				                 vecMDZOut.append(mdz_sphere_extract)
				                 vecThetaOut.append((180./np.pi)*np.arccos(-mdz_sphere_extract))
				                 vecPhiOut.append((180./np.pi)*np.arctan(mdy_sphere_extract/mdx_sphere_extract))
		
		len_same_ev_pore = len(same_ev_pore)
		last_evt_id = same_ev_pore[len_same_ev_pore - 1]
		if (last_evt_id < (len(evt_id_pore)-1)):
			jev_pore = same_ev_pore[len_same_ev_pore - 1] + 1
		else:
			break

print "vecEventIDOut", vecEventIDOut
    		    		
# SPHERE
N_out = len(vecEnergyOut)

print "Exiting protons/entering protons = ", float(N_out)/float(N_in)
print "Real number of protons entering the pore = ", N_in_real


# Write output to file
name_fileout = "./simulation_output/OUTPUT_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_in)+"_"+model+".dat"
print "Writing in "+name_fileout
f_out = open(name_fileout, 'wb')
f_out.write("EventID Energy[keV] MDX MDY MDZ Theta[deg] Phi[deg] \n)
for iel in xrange(len(vecEventIDOut)):
	# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle 
	f_out.write(str(vecEventIDOut[iel])+" "+str(vecEnergyOut[iel])+" "+str(vecMDXOut[iel])+" "+str(vecMDYOut[iel])+" "+str(vecMDZOut[iel])+" "+str(vecThetaOut[iel])+" "+str(vecPhiOut[iel])+"\n")



f_out.close()
hdulist.close()
