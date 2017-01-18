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



from astropy.io import fits
import numpy as np
import math
import sys, os
import matplotlib.pyplot as plt
from matplotlib import gridspec

# Import the input parameters
arg_list = sys.argv
filedir = arg_list[1]
N_fits = int(arg_list[2])
N_in = int(arg_list[3])
theta_0 = float(arg_list[4])
angle_dis = int(arg_list[5])
if (angle_dis == 0): disname = "planar"
if (angle_dis == 1): disname = "cos"
energy_0 = float(arg_list[6])

model = arg_list[7]

HPlate_vol_id_min = 0 
HPlate_vol_id_max = 3 
HRib_vol_id_min = 4 
HRib_vol_id_max = 7 

sphere_vol_id = 100000
ent_vol_id = 200000
exit_vol_id = 300000


N_in_real = 0

vecEventIDOut = []
vecThetaOut = []
vecPhiOut = []
vecEnergyOut = []
vecMDXOut = []
vecMDYOut = []
vecMDZOut = []
vecEventXYZ = []
vecX = []
vecY = []
vecZ = []
nScatt = []

for jfits in xrange(N_fits):
    
    print '%%%%%%%%%%%%%% READING BoGEMMS FILE: ', jfits+1
    hdulist = fits.open(filedir+'/xyz.'+str(jfits)+'.fits.gz')
    
    tbdata = hdulist[1].data
    cols = hdulist[1].columns
    
    evt_id = tbdata.field('EVT_ID')
    trk_id = tbdata.field('TRK_ID')
    vol_id = tbdata.field('VOLUME_ID')
    moth_id = tbdata.field('MOTHER_ID')
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
    z_pore = Z_exit[where_pore]
    y_pore = Y_exit[where_pore]
    x_pore = X_exit[where_pore]
        
    where_sphere = np.where((vol_id == sphere_vol_id) & (part_id == 2212) & (Z_ent < 0) & (moth_id == 0))
    evt_id_sphere = evt_id[where_sphere]
    ene_sphere = ene_ent[where_sphere]
    mdz_sphere = MDZ_ent[where_sphere]
    mdy_sphere = MDY_ent[where_sphere]
    mdx_sphere = MDX_ent[where_sphere]
    z_sphere = Z_ent[where_sphere]
    y_sphere = Y_ent[where_sphere]
    x_sphere = X_ent[where_sphere]

    where_ent = np.where((vol_id == ent_vol_id) & (part_id == 2212) & (moth_id == 0))
    evt_id_ent = evt_id[where_ent]
    ene_ent = ene_ent[where_ent]
    mdz_ent = MDZ_ent[where_ent]
    mdy_ent = MDY_ent[where_ent]
    mdx_ent = MDX_ent[where_ent]
    z_ent = Z_ent[where_ent]
    y_ent = Y_ent[where_ent]
    x_ent = X_ent[where_ent]

    N_in_real += len(evt_id_ent)


    where_exit = np.where((vol_id == exit_vol_id) & (part_id == 2212) & (moth_id == 0))
    evt_id_exit = evt_id[where_exit]
    ene_exit = ene_exit[where_exit]
    mdz_exit = MDZ_exit[where_exit]
    mdy_exit = MDY_exit[where_exit]
    mdx_exit = MDX_exit[where_exit]
    z_exit = Z_exit[where_exit]
    y_exit = Y_exit[where_exit]
    x_exit = X_exit[where_exit]

    if (where_sphere[0].size):
		jev_sphere = 0
		while (1):
			same_ev = np.where(evt_id_sphere == evt_id_sphere[jev_sphere])
			same_ev_ent = np.where(evt_id_ent == evt_id_sphere[jev_sphere])
			same_ev_exit = np.where(evt_id_exit == evt_id_sphere[jev_sphere])
			same_ev_pore = np.where(evt_id_pore == evt_id_sphere[jev_sphere])
		
			same_ev = same_ev[0]
			same_ev_ent = same_ev_ent[0]
			same_ev_exit = same_ev_exit[0]
			same_ev_pore = same_ev_pore[0]
		
			if (same_ev_ent) and (same_ev_exit):
				ene_ent_extract = ene_ent[same_ev_ent]
				mdz_ent_extract = mdz_ent[same_ev_ent]
				mdy_ent_extract = mdy_ent[same_ev_ent]
				mdx_ent_extract = mdx_ent[same_ev_ent]
				z_ent_extract = z_ent[same_ev_ent]
				y_ent_extract = y_ent[same_ev_ent]
				x_ent_extract = x_ent[same_ev_ent]

				ene_exit_extract = ene_exit[same_ev_exit]
				mdz_ent_extract = mdz_exit[same_ev_exit]
				mdy_ent_extract = mdy_exit[same_ev_exit]
				mdx_ent_extract = mdx_exit[same_ev_exit]
				mdz_exit_extract = mdz_exit[same_ev_exit]
				mdy_exit_extract = mdy_exit[same_ev_exit]
				mdx_exit_extract = mdx_exit[same_ev_exit]
				for jmd in xrange(len(mdz_ent_extract)):
				    if (mdz_ent_extract[jmd] != mdz_exit_extract[jmd]): 
				       print "Remizovich in the exit dummy volume (MDZ)! - EventID = "+str(evt_id_sphere[jev_sphere])
				    else:
				       if (mdy_ent_extract[jmd] != mdy_exit_extract[jmd]):
				          print "Remizovich in the exit dummy volume (MDY)! - EventID = "+str(evt_id_sphere[jev_sphere])
				       else:
				          if (mdx_ent_extract[jmd] != mdx_exit_extract[jmd]):
				              print "Remizovich in the exit dummy volume (MDX)! - EventID = "+str(evt_id_sphere[jev_sphere])
				
				z_exit_extract = z_exit[same_ev_exit]
				y_exit_extract = y_exit[same_ev_exit]
				x_exit_extract = x_exit[same_ev_exit]
					
				ene_sphere_extract = ene_sphere[jev_sphere]
				mdz_sphere_extract = mdz_sphere[jev_sphere]
				mdy_sphere_extract = mdy_sphere[jev_sphere]
				mdx_sphere_extract = mdx_sphere[jev_sphere] 
				z_sphere_extract = z_sphere[jev_sphere]
				y_sphere_extract = y_sphere[jev_sphere]
				x_sphere_extract = x_sphere[jev_sphere] 
			
				if (len(same_ev_pore)):
					ene_pore_extract = ene_pore[same_ev_pore]
					mdz_pore_extract = mdz_pore[same_ev_pore]
					mdy_pore_extract = mdy_pore[same_ev_pore]
					mdx_pore_extract = mdx_pore[same_ev_pore]
					z_pore_extract = z_pore[same_ev_pore]
					y_pore_extract = y_pore[same_ev_pore]
					x_pore_extract = x_pore[same_ev_pore]
				else:
					print "Proton leakage at EventID "+str(evt_id_sphere[jev_sphere])+"!"
			
				for jsame in xrange(len(ene_exit_extract)):
					if (ene_exit_extract[jsame] == ene_sphere_extract):
						if (mdz_exit_extract[jsame] == mdz_sphere_extract):
							if (mdy_exit_extract[jsame] == mdy_sphere_extract):
								if (mdx_exit_extract[jsame] == mdx_sphere_extract):
									 vecEventIDOut.append(evt_id_sphere[jev_sphere])
									 vecEnergyOut.append(ene_sphere_extract)
									 vecMDXOut.append(mdx_sphere_extract)
									 vecMDYOut.append(mdy_sphere_extract)
									 vecMDZOut.append(mdz_sphere_extract)
									 vecThetaOut.append((180./np.pi)*np.arccos(-mdz_sphere_extract))
									 if ((mdx_sphere_extract >= 0.) & (mdy_sphere_extract >= 0.)): 
									     vecPhiOut.append((180./np.pi)*np.arctan(mdy_sphere_extract/mdx_sphere_extract)) - 180.)
									 if ((mdx_sphere_extract < 0.) & (mdy_sphere_extract >= 0.)): 
									     vecPhiOut.append((180./np.pi)*np.arctan(mdy_sphere_extract/mdx_sphere_extract))
									 if ((mdx_sphere_extract < 0.) & (mdy_sphere_extract < 0.)): 
									     vecPhiOut.append((180./np.pi)*np.arctan(mdy_sphere_extract/mdx_sphere_extract))
									 if ((mdx_sphere_extract >= 0.) & (mdy_sphere_extract < 0.)): 
									     vecPhiOut.append((180./np.pi)*np.arctan(mdy_sphere_extract/mdx_sphere_extract) + 180.))

									 if (len(same_ev_pore)):
										 xlist = [x_ent_extract.tolist(), x_pore_extract.tolist(), x_exit_extract.tolist(), [x_sphere_extract]]
										 ylist = [y_ent_extract.tolist(), y_pore_extract.tolist(), y_exit_extract.tolist(), [y_sphere_extract]]
										 zlist = [z_ent_extract.tolist(), z_pore_extract.tolist(), z_exit_extract.tolist(), [z_sphere_extract]]
										 nScatt.append(len(x_pore_extract))
										 for jxyz in xrange(len(xlist)):
											elx = xlist[jxyz]
											ely = ylist[jxyz]
											elz = zlist[jxyz]
											for jel in xrange(len(elx)):
											   vecEventXYZ.append(evt_id_sphere[jev_sphere])
											   vecX.append(elx[jel])
											   vecY.append(ely[jel])
											   vecZ.append(elz[jel])
									 
										
									 else:
										 xlist = [x_ent_extract.tolist(), x_exit_extract.tolist(), [x_sphere_extract]]
										 ylist = [y_ent_extract.tolist(), y_exit_extract.tolist(), [y_sphere_extract]]
										 zlist = [z_ent_extract.tolist(), z_exit_extract.tolist(), [z_sphere_extract]]
										 nScatt.append(0)
										 for jxyz in xrange(len(xlist)):
											elx = xlist[jxyz]
											for jel in xrange(len(elx)):
											   vecEventXYZ.append(evt_id_sphere[jev_sphere])
											   vecX.append(elx[jel])
											   vecY.append(elx[jel])
											   vecZ.append(elx[jel])
		
			len_same_ev = len(same_ev)
			last_evt_id = same_ev[len_same_ev - 1]
			if (last_evt_id < (len(evt_id_sphere)-1)):
				jev_sphere = same_ev[len_same_ev - 1] + 1
			else:
				break
    

    		    		
N_out = len(vecEnergyOut)


print "Real number of protons entering the pore = ", N_in_real
print "Exiting protons/entering protons = ", float(N_out)/float(N_in_real)

path_simoutput = './simulation_output'
if not os.path.exists(path_simoutput):
	os.makedirs(path_simoutput)

# Write output to file
name_fileout = path_simoutput+"/OUTPUT_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_in)+"_"+model+".dat"
print "Writing in "+name_fileout
f_out = open(name_fileout, 'wb')
f_out.write("# N_in: "+str(N_in_real)+" \n")
f_out.write("# Exiting protons/entering protons: "+str(float(N_out)/float(N_in_real))+" \n")
f_out.write("# EventID Energy[keV] MDX MDY MDZ Theta[deg] Phi[deg] nScatterings N_out N_in_real \n")
for iel in xrange(len(vecEventIDOut)):
	# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle 
	f_out.write(str(vecEventIDOut[iel])+" "+str(vecEnergyOut[iel])+" "+str(vecMDXOut[iel])+" "+str(vecMDYOut[iel])+" "+str(vecMDZOut[iel])+" "+str(vecThetaOut[iel])+" "+str(vecPhiOut[iel])+" "+str(nScatt[iel])+" "+str(N_out)+" "+str(N_in_real)+"\n")

# Write positions to FITS file
name_fileout = path_simoutput+"/XYZ_"+str(int(energy_0))+"keV_"+disname+str(theta_0)+"deg_"+str(N_in)+"_"+model+".fits.gz"
print "Writing in "+name_fileout
xyz_col1 = fits.Column(name='EVENT_ID', format='1J', array=vecEventXYZ)
xyz_col2 = fits.Column(name='X', format='1D', array=vecX)
xyz_col3 = fits.Column(name='Y', format='1D', array=vecY)
xyz_col4 = fits.Column(name='Z', format='1D', array=vecZ)

xyz_cols = fits.ColDefs([xyz_col1, xyz_col2, xyz_col3, xyz_col4])
xyz_tbhdu = fits.BinTableHDU.from_columns(xyz_cols)
xyz_tbhdu.writeto(name_fileout, clobber=1)



f_out.close()
hdulist.close()
