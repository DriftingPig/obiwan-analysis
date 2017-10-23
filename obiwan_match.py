'''
This code matches obiwan-tractor file with the obiwan-simcat file and the DR_5 file, to get and validate the randoms within the tractor files
'''

import os 
import glob
import numpy as n
from astropy.io import fits
from math import *
# obiwan outputs two patches currently :

# ra 175 outputs
# /global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175
top_dir_175 = "/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175/elg/"
brick_pathes_175 = n.array(glob.glob(top_dir_175 + "???/*"))
brick_pathes_175.sort()
brick_list_175 = n.array([ os.path.basename(br) for br in brick_pathes_175 ])

# ra 125 outputs
# /global/cscratch1/sd/kaylanb/obiwan_out/obiwan_elg_9deg/elg/
top_dir_125 = "/global/cscratch1/sd/kaylanb/obiwan_out/obiwan_elg_9deg/elg/" 
brick_pathes_125 = n.array(glob.glob(top_dir_125 + "???/*"))
brick_pathes_125.sort()
brick_list_125 = n.array([ os.path.basename(br) for br in brick_pathes_125 ])

# choose one of the two dirs above 
top_dir = top_dir_175
brick_pathes = brick_pathes_175
brick_list = brick_list_175

# choose a brick to process in the list above
index = 0
anglelim=0.00028

            
def select_all(index, top_dir = top_dir, brick_pathes = brick_pathes, brick_list = brick_list):
    brick_name = brick_list[index]
    brick_path = brick_pathes[index]
    #simcat files
    path_2_tractor_sim_files = n.array(glob.glob(os.path.join(brick_path, "*rs*", "obiwan", "simcat-elg-"+brick_name+".fits")))
    path_2_tractor_sim_files.sort()
    #tractor files
    path_2_tractor_random_files = n.array(glob.glob(os.path.join(brick_path, "*rs*", "tractor", "tractor-"+brick_name+".fits")))
    path_2_tractor_random_files.sort()
    #DR_5 data
    top_dir_dr5 = "/global/project/projectdirs/cosmo/data/legacysurvey/dr5/tractor"
    path_2_tractor_data_file = os.path.join( top_dir_dr5, brick_name[:3], "tractor-"+brick_name+".fits")
    if len(path_2_tractor_random_files)==0:
        return False, 0., False, 0., 0., 0., 0.
    if len(path_2_tractor_sim_files)==0:
        raise ValueError("lack of simcat files! brickname = ",brickname)
    if os.path.isfile(path_2_tractor_data_file) is False:
        raise ValueError("lack of data files in DR5! brickname = ",brickname)
        #return False, 0., False, 0., 0., 0., 0. #use this command if the data files are not complete
    #select data files
    flag_D, datN, datS = select_ELG(path_2_tractor_data_file)
    if flag_D:
        print("number of galaxies NGC selection, SGC selection:", len(datN), len(datS))

    flags = [flag_D]
    ii=0
    for path in path_2_tractor_random_files[:] :
        flag_R, rdN_i, rdS_i = select_ELG(path)
        flags.append(flag_R)
        if ii==0:
            rdN = rdN_i
            rdS = rdS_i
        else:
            rdN = n.hstack((rdN,rdN_i))
            rdS = n.hstack((rdS,rdS_i))
        ii+=1
    ii=0
    print("number of randoms NGC selection, SGC selection:", len(rdN), len(rdS))
    for path in path_2_tractor_sim_files[:] :
        dat=fits.open(path)[1].data
        if ii==0:
            sim = dat
        else:
            sim = n.hstack((sim,dat))
        ii+=1
    if len(rdN)==0 and len(rdS)==0:
        flag_R=False
    else:
        flag_R=True
    if flag_D is False and flag_R is False:
        print("no ELGs (DR5 and Random) in this brick,",brick_name)
        return False, sim , False, datN, datS, rdN, rdS #FF:no ELGs in this brick
    if flag_D is False:
        print ("no DR5 ELGs in this brick,",brick_name)
        return False, sim , True, datN, datS, rdN, rdS #FT: no DR5 ELGs in this brick
    if flag_R is False:
        print("no random ELGs in this brick,",brick_name)
        return True, sim , False, datN, datS, rdN, rdS #TF: no random ELGs in this brick
    return True, sim, True, datN, datS, rdN, rdS #TT: all good
        
def select_ELG( path_2_tractor_file ):
	"""
	Given the path to a tractor catalog, it returns two sub catalogs with the eBOSS ELG selections applied (NGC and SGC).
	"""
	# opens the tractor file
	hdu=fits.open(path_2_tractor_file) 
	dat=hdu[1].data
	hdu.close()
	# the color color selection
	g     = 22.5 - 2.5 * n.log10(dat['flux_g'] / dat['mw_transmission_g'])
	r_mag = 22.5 - 2.5 * n.log10(dat['flux_r'] / dat['mw_transmission_r'])
	z_mag = 22.5 - 2.5 * n.log10(dat['flux_z'] / dat['mw_transmission_z'])
	gr = g - r_mag
	rz = r_mag - z_mag
	color_sgc = (g>21.825)&(g<22.825)&(-0.068*rz+0.457<gr)&(gr< 0.112*rz+0.773) &(0.218*gr+0.571<rz)&(rz<-0.555*gr+1.901)
	color_ngc = (g>21.825)&(g<22.9)  &(-0.068*rz+0.457<gr)&(gr< 0.112*rz+0.773) &(0.637*gr+0.399<rz)&(rz<-0.555*gr+1.901) 
	# the junk rejection criterion
	noJunk = (dat['brick_primary']) & (dat['anymask_g']==0) & (dat['anymask_r']==0) & (dat['anymask_z']==0) #& (dat['TYCHO2INBLOB']==False)
	# the low depth region rejection
	value_g=dat['psfdepth_g']
	value_r=dat['psfdepth_r']
	value_z=dat['psfdepth_z']
	gL = 62.79716079 
	rL = 30.05661087
	zL_ngc = 11.0
	zL_sgc = 12.75  
	depth_selection_ngc = (value_g > gL) & (value_r > rL) & (value_z > zL_ngc)
	depth_selection_sgc = (value_g > gL) & (value_r > rL) & (value_z > zL_sgc)
	# final selection boolean array :
	selection_sgc =(noJunk)&(color_sgc)&(depth_selection_sgc)
	selection_ngc =(noJunk)&(color_ngc)&(depth_selection_ngc)
	# returns the catalogs of ELGs
	if len(selection_sgc.nonzero()[0])>0 or  len(selection_ngc.nonzero()[0])>0 :
		flag = True
		return flag, dat[selection_ngc], dat[selection_sgc]
	else :
		flag = False
		return flag, dat[selection_ngc], dat[selection_sgc]

def sort_brick(index):
    from astropy.table import Table
    N_bricks = 0
    flag_D, elg_sim, flag_R, datN, datS, rdN, rdS = select_all(index)
    #if sim_flag:
        #print len(datN),len(datS)
    if elg_sim == 0 :
        return False, 0., 0., 0., 0., 0., 0., 0.
    
    elg_sim=n.array(elg_sim)
    coldefs_elg_sim = fits.ColDefs(elg_sim)
    sim_hdu = fits.BinTableHDU.from_columns(coldefs_elg_sim)
    tb_sim = Table.read(sim_hdu)
    tb_sim.sort('ra')
    
    datN=n.array(datN)
    coldefs_datN = fits.ColDefs(datN)
    datN_hdu = fits.BinTableHDU.from_columns(coldefs_datN)
    tb_datN = Table.read(datN_hdu)
    tb_datN.sort('ra')

    datS=n.array(datS)
    coldefs_datS = fits.ColDefs(datS)
    datS_hdu = fits.BinTableHDU.from_columns(coldefs_datS)
    tb_datS = Table.read(datS_hdu)
    tb_datS.sort('ra')

    rdN=n.array(rdN)
    coldefs_rdN = fits.ColDefs(rdN)
    rdN_hdu = fits.BinTableHDU.from_columns(coldefs_rdN)
    tb_rdN = Table.read(rdN_hdu)
    tb_rdN.sort('ra')
            
    rdS=n.array(rdS)
    coldefs_rdS = fits.ColDefs(rdS)
    rdS_hdu = fits.BinTableHDU.from_columns(coldefs_rdS)
    tb_rdS = Table.read(rdS_hdu)
    tb_rdS.sort('ra')
    return True, tb_sim,tb_datN,tb_datS,tb_rdN,tb_rdS, datN, datS


def AngD(tb1,i,tb2,j):
    phi1 = tb1['ra'][i]*pi/180.0
    phi2 = tb2['ra'][j]*pi/180.0
    thi1 = tb1['dec'][i]*pi/180.0
    thi2 = tb2['dec'][j]*pi/180.0
    if cos(thi1)*cos(thi2)*cos(phi1-phi2)+sin(thi1)*sin(thi2)>=1.:
        return 0.
    return (acos(cos(thi1)*cos(thi2)*cos(phi1-phi2)+sin(thi1)*sin(thi2)))*180./pi

test = 0
def obiwan_random_match(index):#or 175
    flag, tb_sim,tb_datN,tb_datS,tb_rdN,tb_rdS, datN, datS = sort_brick(index)


    if flag is False:
        print flag
        return False, 0, 0, 0, 0, 0, 0, 0, 0,0.,0.
    '''
    #get the z column for simcat catalogues
    from astropy.table import Column
    import sys
    sys.path.append('/global/cscratch1/sd/huikong/obiwan_code/obiwan')
    from obiwan.db_tools import redshifts_for_ids
    z_list=[]
    for i in range(0,len(tb_sim)):
        ids = [tb_sim['id'][i]]
        rdsft = redshifts_for_ids(ids,db_table='obiwan_elg_ra175')
        z_list.append(rdsft[1][0])
    col_z = Column(n.array(z_list), name='z')
    tb_sim.add_column(col_z, index=0)
    import matplotlib.pyplot as plt
    plt.hist(tb_sim['z'])
    plt.show()
    '''
    
    from astropy.table import Column
    up_flag = 0
    update_flag = 1
    rdN_id_list = -n.ones(len(tb_rdN))
    matched = False
    rdN_total_match = 0
    import math
    for i in range(0,len(tb_rdN)):
        j=up_flag
        while j<len(tb_sim):
            if update_flag:
                if (tb_rdN['ra'][i]-tb_sim['ra'][j])*math.fabs(math.cos(((tb_rdN['dec'][i]+tb_sim['dec'][j])/2)*pi/180.))<anglelim*2 :
                    update_flag = 0
                    up_flag = j
                    j-=1
            else:
                if AngD(tb_rdN,i,tb_sim,j)<anglelim:
                        if matched :
                            #print 'duplicated result encountered! '
                            if AngD(tb_rdN,i,tb_sim,j)<rdN_id_list[i]:
                                rdN_id_list[i] = AngD(tb_rdN,i,tb_sim,j)
                        else:
                            rdN_id_list[i] = AngD(tb_rdN,i,tb_sim,j)
                            matched = True
                            rdN_total_match+=1
                            if len(n.where(n.array(rdN_id_list)>=0)[0]) != rdN_total_match :
                                raise ValueError('do no match!')
                            
                if (tb_sim['ra'][j]-tb_rdN['ra'][i])*math.fabs(math.cos(((tb_sim['dec'][j]+tb_rdN['dec'][i])/2)*pi/180.))>anglelim*2 :
                    update_flag = 1
                    matched = False
                    break
                if j == len(tb_sim)-1 :
                    update_flag = 1
                    matched = False
            j+=1
                    
    up_flag = 0
    update_flag = 1
    rdN_id_list2 = -n.ones(len(tb_rdN))
    matched = False
    rdN_total_match2 = 0    
    for i in range(0,len(tb_rdN)):
        j=up_flag
        while j<len(tb_datN):
            if update_flag:
                if (tb_rdN['ra'][i]-tb_datN['ra'][j])*math.fabs(math.cos(((tb_rdN['dec'][i]+tb_datN['dec'][j])/2)*pi/180.))<anglelim :
                    update_flag = 0
                    up_flag = j
                    j-=1
            else:
                    if AngD(tb_rdN,i,tb_datN,j)<anglelim:
                            if matched :
                                #print 'duplicated result encountered#! '+str(i)
                                if AngD(tb_rdN,i,tb_datN,j)<rdN_id_list2[i]:
                                    rdN_id_list2[i] = AngD(tb_rdN,i,tb_datN,j)
                            else:
                                rdN_id_list2[i] = AngD(tb_rdN,i,tb_datN,j)
                                matched = True
                                rdN_total_match2+=1
                    if (tb_datN['ra'][j]-tb_rdN['ra'][i])*math.fabs(math.cos(((tb_rdN['dec'][i]+tb_datN['dec'][j])/2)*pi/180.))>anglelim*2 :
                        update_flag = 1
                        matched = False
                        break
                    if j == len(tb_datN)-1 :
                        update_flag = 1
                        matched = False
            j+=1
                        
 
    up_flag = 0
    update_flag = 1
    rdS_id_list = -n.ones(len(tb_rdS))
    matched = False
    rdS_total_match = 0
    import math
    for i in range(0,len(tb_rdS)):
        j=up_flag
        while j<len(tb_sim):
            if update_flag:
                if (tb_rdS['ra'][i]-tb_sim['ra'][j])*math.fabs(math.cos(((tb_rdS['dec'][i]+tb_sim['dec'][j])/2)*pi/180.))<anglelim :
                    update_flag = 0
                    up_flag = j
                    j-=1
            else:
                if AngD(tb_rdS,i,tb_sim,j)<anglelim:
                        if matched :
                            #print 'duplicated result encountered! '+str(i)
                            if AngD(tb_rdS,i,tb_sim,j)<rdS_id_list[i] :
                                    rdS_id_list[i] = AngD(tb_rdS,i,tb_sim,j)
                        else:
                            rdS_id_list[i] = AngD(tb_rdS,i,tb_sim,j)
                            matched = True
                            rdS_total_match+=1
                if (tb_sim['ra'][j]-tb_rdS['ra'][i])*math.fabs(math.cos(((tb_rdS['dec'][i]+tb_sim['dec'][j])/2)*pi/180.))>anglelim*2 :
                    update_flag = 1
                    matched = False
                    break
                if j == len(tb_sim)-1 :
                    update_flag = 1
                    matched = False
            j+=1
                    
    up_flag = 0
    update_flag = 1
    rdS_id_list2 = -n.ones(len(tb_rdS))
    matched = False
    rdS_total_match2 = 0    
    for i in range(0,len(tb_rdS)):
        j=up_flag
        while j<len(tb_datS):
            if update_flag:
                if (tb_rdS['ra'][i]-tb_datS['ra'][j])*math.fabs(math.cos(((tb_rdS['dec'][i]+tb_datS['dec'][j])/2)*pi/180.))<anglelim :
                    update_flag = 0
                    up_flag = j
                    j-=1
            else:
                    if AngD(tb_rdS,i,tb_datS,j)<anglelim:
                            if matched :
                                #print 'duplicated result encountered##! '+str(i)
                                if AngD(tb_rdS,i,tb_datS,j)<rdS_id_list2[i] :
                                    rdS_id_list2[i] = AngD(tb_rdS,i,tb_datS,j)
                            else:
                                rdS_id_list2[i] = AngD(tb_rdS,i,tb_datS,j)
                                matched = True
                                rdS_total_match2+=1
                    if (tb_datS['ra'][j]-tb_rdS['ra'][i])*math.fabs(math.cos(((tb_rdS['dec'][i]+tb_datS['dec'][j])/2)*pi/180.))>anglelim*2 :
                        update_flag = 1
                        matched = False
                        break
                    if j == len(tb_datS)-1 :
                        update_flag = 1
                        matched = False
            j+=1
                        
    
    col_sim_flag = Column(rdN_id_list, name='sim_flag')
    tb_rdN.add_column(col_sim_flag, index=0)
    col_dat_flag = Column(rdN_id_list2, name='dat_flag')
    tb_rdN.add_column(col_dat_flag, index=0)
    
    col_sim_flag = Column(rdS_id_list, name='sim_flag')
    tb_rdS.add_column(col_sim_flag, index=0)
    col_dat_flag = Column(rdS_id_list2, name='dat_flag')
    tb_rdS.add_column(col_dat_flag, index=0)
    
    print rdN_total_match2#dat match
    print rdN_total_match#sim match
    print len(tb_rdN)

    print rdS_total_match2#dat match
    print rdS_total_match#sim match
    print len(tb_rdS)
    
    return True, tb_rdN, tb_rdS, rdN_total_match2, rdN_total_match, len(tb_rdN), rdS_total_match2, rdS_total_match, len(tb_rdS),datN, datS
                    

    
flag=False
N_bricks=0
count = 0
while flag is False:
    print '#'+str(count)
    flag, tb_rdN, tb_rdS, rdN_total_match2, rdN_total_match, len_tb_rdN, rdS_total_match2, rdS_total_match, len_tb_rdS, datN, datS = obiwan_random_match(count)
    count += 1

print count
for ii in range(count,len(brick_list)):
    print '#'+str(ii)
    flag, tb_rdN_i,tb_rdS_i, rdN_total_match2_i, rdN_total_match_i, len_tb_rdN_i, rdS_total_match2_i, rdS_total_match_i, len_tb_rdS_i, datN_i, datS_i = obiwan_random_match(ii)
    if flag:
        N_bricks += 1
        tb_rdN = n.hstack((tb_rdN,tb_rdN_i))
        tb_rdS = n.hstack((tb_rdS,tb_rdS_i))
        rdN_total_match2+=rdN_total_match2_i
        rdN_total_match+= rdN_total_match_i
        len_tb_rdN+=len_tb_rdN_i
        rdS_total_match2+=rdS_total_match2_i
        rdS_total_match+=rdS_total_match_i
        len_tb_rdS+=len_tb_rdS_i
        datN = n.hstack((datN,datN_i))
        datS = n.hstack((datS,datS_i))
        selection = tb_rdN['sim_flag']>=0
        if len(tb_rdN[selection]) != rdN_total_match :
            print 'sim_flag='+str(len(tb_rdN[selection]))+' total_match='+str(rdN_total_match)
            print len(tb_rdN)
            raise ValueError('do not match!2')
        
print str(rdN_total_match2)+' '+str(rdN_total_match)+' '+str(len_tb_rdN)+'\n'
print str(rdS_total_match2)+' '+str(rdS_total_match)+' '+str(len_tb_rdS)+'\n'  
coldefs_elg_rdN_flagged = fits.ColDefs(tb_rdN)
elg_rdN_flagged_hdu = fits.BinTableHDU.from_columns(coldefs_elg_rdN_flagged)
elg_rdN_flagged_hdu.writeto('./output/elg_rdN_flagged175.fits',overwrite = True)

coldefs_elg_rdS_flagged = fits.ColDefs(tb_rdS)
elg_rdS_flagged_hdu = fits.BinTableHDU.from_columns(coldefs_elg_rdS_flagged)
elg_rdS_flagged_hdu.writeto('./output/elg_rdS_flagged175.fits',overwrite = True)

coldefs_datN = fits.ColDefs(datN)
elg_datN_hdu = fits.BinTableHDU.from_columns(coldefs_datN)
elg_datN_hdu.writeto('./output/elg_N175.fits',overwrite = True)

coldefs_datS = fits.ColDefs(datS)
elg_datS_hdu = fits.BinTableHDU.from_columns(coldefs_datS)
elg_datS_hdu.writeto('./output/elg_S175.fits',overwrite = True)

