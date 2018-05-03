#################################################
#SMA CASA Continuum Reduction Script
#################################################
'''CASA script to be used for the flagging, calibration and imaging of SMA Continuum Data
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) Calibrated Split MS for each side-band (and full band) -- [target]_[obsDate]_[band]_[sideband].ms
        (2) Continuum images in each side-band (and full band) -- [target]_[obsDate]_[band]_[sideband]_clean1.image(.tt0).pbcor
        (3) File of flux densities from image/UV plane fitting -- imfit_results.txt/uvfit_results.txt
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
       - All output images are also converted to fits format (just append .fits to end of images above)
       - All input logged in user_input.log.
       - This script is intended to be used with raw data that has been converted to CASA MS format;
	   	 There are 2 methods to get from raw data to CASA MS
         (1) -see instructions in how_to_sma_scripts.txt (recommended!)
         (2) Alternatively, the old way to convert is through miriad. 
         If you have can create raw miriad files (from idl2miriad task in MIR),
         run miriad bash script (miriad2fits.sh) and run CASA script (fits2casa.py).
WARNING: If you have SMA data calibrated in MIR or MIRIAD, use instructions here,
https://www.cfa.harvard.edu/rtdc/SMAdata/process/casa/convertcasa/
to convert to CASA MS for imaging.

Written by: Alex J. Tetarenko
Last Updated: April 18 2018
Works in CASA-5.1.1 now!'''

print '##################################################'
print 'Welcome to Alexs SMA Continuum Reduction Script'
print '##################################################\n'

#packages you may need
from astropy.io import fits
import numpy as np
import re
import sys
import imp
import os
import linecache
import find
import warnings
import webbrowser
import datetime
warnings.filterwarnings('ignore')
from utils_CASA import imfit_point, phselfcal, writeDict
from ekoch_casa_tools import set_imagermode,has_field,set_cellsize,set_imagesize,find_expected_beams,getBaselinePercentile,get_mosaic_info
import uvmultifit as uvm
from astropy.io import ascii

#define output directory
my_dir='/mnt/bigdata/tetarenk/SMA_maxi1820/raw_data/final_MS/rec230/'
if not os.path.isdir(my_dir):
	os.system('sudo mkdir '+my_dir)
	os.system('sudo chown ubuntu '+my_dir)
	os.system('sudo chmod -R u+r '+my_dir) 
	os.system('sudo chmod -R u+w '+my_dir)
	os.system('sudo chmod -R u+x '+my_dir)
print 'You have set your output directory to ', my_dir
print 'All output images & intermediate data products are put in this directory.\n'

#param file location
param_dir_file='/mnt/bigdata/tetarenk/SMA_maxi1820/raw_data/final_MS/rec230/params_sma.txt'
print 'You have set your param file to ', param_dir_file
print 'Please make sure all parameters are correct, they will change for each data set!\n'

#make dictionary to write all user input variables to for log
dict_log=[]

#################################################
#Defining Parameters Section
#################################################
print 'Reading in parameters...'
#read in param file
def getVar(filename):
	'''Easy way to read in variables from parameter file'''
	f=open(filename)
	global data_params
	data_params=imp.load_source('data_params','',f)
	f.close()
getVar(param_dir_file)
#data set params
ms_name_lsb=data_params.ms_name_lsb
ms_name_usb=data_params.ms_name_usb
obsDate=data_params.obsDate
target=data_params.target
fields_lsb=data_params.fields_lsb
spw_lsb=data_params.spw_lsb
scans_lsb=data_params.scans_lsb
fields_usb=data_params.fields_usb
spw_usb=data_params.spw_usb
scans_usb=data_params.scans_usb
remakems=data_params.remakems
doImage=data_params.doImage
bandsIM=data_params.bandsIM
do_ant_correct=data_params.do_ant_correct
ant_corr_file=data_params.ant_corr_file
#general image params
use_auto=data_params.use_auto
mythreshold=data_params.mythreshold
myimsize=data_params.myimsize
mycell=data_params.mycell
myniter=data_params.myniter
mynterms=data_params.mynterms
multiscale=data_params.multiscale
robust=data_params.robust
weighting=data_params.weighting
#mask options
mymask=data_params.mymask
#uv fitting
uv_fit=data_params.uv_fit

#write all variables to log dictionary
dict_log.extend([('ms_name_lsb',ms_name_lsb),('ms_name_usb',ms_name_usb),('obsDate',obsDate),('target',target),\
	('fields_lsb',fields_lsb),('spw_lsb',spw_lsb),('scans_lsb',scans_lsb),('fields_usb',fields_usb),('spw_usb',spw_usb),\
	('scans_usb',scans_usb),('remakeMS',remakems),('mythreshold',mythreshold),('myimsize',myimsize),('mycell',mycell),('ant file',ant_corr_file),('ant_correct',do_ant_correct),\
	('myniter',myniter),('mynterms',mynterms),('mymask',mymask),('uv_fit',uv_fit),('doImage',doImage),('use_auto',use_auto),('bandsIM',bandsIM)])
#################################################


#################################################
#Examining and Splitting Data Section
#################################################
#listobs first
if not os.path.isfile(my_dir+obsDate+'_lsb_listfile.txt'):
	print 'Making lsb listobs file...'
	listobs(ms_name_lsb,listfile=my_dir+obsDate+'_lsb_listfile.txt')
else:
	print 'Listobs lsb file already exists'
if not os.path.isfile(my_dir+obsDate+'_usb_listfile.txt'):
	print 'Making usb listobs file...'
	listobs(ms_name_usb,listfile=my_dir+obsDate+'_usb_listfile.txt')
else:
	print 'Listobs usb file already exists'
seelo=raw_input('Do you want to see the listobs? y or n-->')
if seelo=='y':
	os.system('pluma '+my_dir+obsDate+'_lsb_listfile.txt &')
	os.system('pluma '+my_dir+obsDate+'_usb_listfile.txt &')
	raw_input('Please press enter to continue when you are done.')
else:
	print 'Okay. Moving on.'

#make split data sets of lsb and usb
if not os.path.isdir(my_dir+'data_'+obsDate+'_'+'lsb.ms'):
	print 'Making split MS for lsb band...'
	split(vis=ms_name_lsb,outputvis=my_dir+'data_'+obsDate+'_'+'lsb.ms',\
		spw=spw_lsb,datacolumn='data',scan=scans_lsb,field=fields_lsb)
elif remakems=='T':
	print 'Remaking split MS for lsb band...'
	os.system('rm -rf '+my_dir+'data_'+obsDate+'_'+'lsb.ms')
	os.system('rm -rf '+my_dir+'data_'+obsDate+'_'+'lsb.ms.flagversions')
	split(vis=ms_name_lsb,outputvis=my_dir+'data_'+obsDate+'_'+'lsb.ms',\
		spw=spw_lsb,datacolumn='data',scan=scans_lsb,field=fields_lsb)
else:
	print 'Split MS for lsb already exists.'
if not os.path.isdir(my_dir+'data_'+obsDate+'_'+'usb.ms'):
	print 'Making split MS for usb band...'
	split(vis=ms_name_usb,outputvis=my_dir+'data_'+obsDate+'_'+'usb.ms',\
		spw=spw_usb,datacolumn='data',scan=scans_usb,field=fields_usb)
elif remakems=='T':
	print 'Remaking split MS for usb band...'
	os.system('rm -rf '+my_dir+'data_'+obsDate+'_'+'usb.ms')
	os.system('rm -rf '+my_dir+'data_'+obsDate+'_'+'usb.ms.flagversions')
	split(vis=ms_name_usb,outputvis=my_dir+'data_'+obsDate+'_'+'usb.ms',\
		spw=spw_usb,datacolumn='data',scan=scans_usb,field=fields_usb)
else:
	print 'Split MS for usb already exists.'
#################################################

print 'Beginning reduction process...'
print 'Running listobs.'
print 'Please examine files for fields ids & spw ranges.'
print 'Also check channel resolution same for all spws. If not may need mstransform.'
ms_namel=my_dir+'data_'+obsDate+'_'+'lsb.ms'
ms_nameu=my_dir+'data_'+obsDate+'_'+'usb.ms'
ms_namel_prefix=ms_namel.strip('.ms')
ms_nameu_prefix=ms_nameu.strip('.ms')
os.system('rm -rf '+ms_namel_prefix+'_listfile.txt')
listobs(ms_namel,listfile=ms_namel_prefix+'_listfile.txt')
os.system('pluma '+ms_namel_prefix+'_listfile.txt &')
os.system('rm -rf '+ms_nameu_prefix+'_listfile.txt')
listobs(ms_nameu,listfile=ms_nameu_prefix+'_listfile.txt')
os.system('pluma '+ms_nameu_prefix+'_listfile.txt &')
raw_input('Please press enter to continue when you are done.')
#define variables specific to each band
print 'Please enter the following data set specifics:'
bpf_cal=raw_input('Please enter field id for bandpass cal (if >1 seperate by comma)-->')
bpf_lst=bpf_cal.split(',')
second_cal=raw_input('Please enter field id for secondary cal (if >1 seperate by comma)-->')
second_lst=second_cal.split(',')
flux_cal=raw_input('Please enter field id for flux cal (if >1 seperate by comma)-->')
flux_lst=flux_cal.split(',')
target_id=raw_input('Please enter field id for target-->')
#timerbp=raw_input('Please enter the time range of the bandpass cal scan(s). e.g., 11:43:30.0~11:47:27.0-->')
field_lst=[]
[field_lst.append(i) for i in bpf_lst]
[field_lst.append(i) for i in second_lst]
[field_lst.append(i) for i in flux_lst]
spw_low=raw_input('Please enter lower side-band spw range. e.g., 0~15; Do not include dummy spw!-->')
first_lsb_spw=int(spw_low.split('~')[0])
spw_high=raw_input('Please enter upper side-band spw range. e.g., 0~15; Do not include dummy spw!-->')
first_usb_spw=int(spw_high.split('~')[0])
band='1mm'
band_low='219GHz'
band_high='231GHz'
date=obsDate+'_'+band
ant_plot=my_dir+'antennas_'+date+'.png'
cal_table_prefixl=my_dir+target+'_'+date+'_lsb'
cal_table_prefixu=my_dir+target+'_'+date+'_usb'
split_low=cal_table_prefixl+'_'+band_low+'.ms'
split_high=cal_table_prefixu+'_'+band_high+'.ms'
split_full=my_dir+target+'_'+date+'_comb.ms'
#add to log dictionary
dict_log.append(('bpf_cal',bpf_cal))
dict_log.append(('second_cal',second_cal))
dict_log.append(('flux_cal',flux_cal))
dict_log.append(('target_id',target_id))
dict_log.append(('spw_low',spw_low))
dict_log.append(('spw_high',spw_high))
dict_log.append(('band',band))
dict_log.append(('band_low',band_low))
dict_log.append(('band_high',band_high))
#########################################
#Reference antenna and flagging
#########################################
#plot antenna positions to determne reference antenna near centre of array that 
#is present throughout whole observation and has no problems
print 'Plotting antenna positions. Please choose a reference antenna near centre of array.'
os.system('rm -rf '+ant_plot)
plotants(vis=ms_namel,figfile=ant_plot)
ref_ant=raw_input('Please enter reference antenna. e.g., 2-->')
dict_log.append(('ref_ant',ref_ant))
skipflag=raw_input('Is the data set already flagged?y or n-->')
dict_log.append(('skip_flag',skipflag))
if skipflag=='n':
	#flagging
	print 'Starting Flagging...'
	#dummy scans
	flag_dummy=raw_input('Do you want to flag the dummy scans?y or n-->')
	dict_log.append(('flag_dummy',flag_dummy))
	if flag_dummy=='y':
		dum_scanl,dumscanu=raw_input('Please enter dummy scan ids for lsb and usb. e.g., 1 2-->').split(' ')
		dict_log.append(('dum_scanl',dum_scanl))
		dict_log.append(('dum_scanu',dum_scanu))
		print 'Flagging dummy scans.'
		flagdata(vis=ms_namel, flagbackup=True, mode='manual', scan=dum_scanl)
		flagdata(vis=ms_nameu, flagbackup=True, mode='manual', scan=dum_scanu)
	else:
		print 'Dummy scans not flagged.'
	#first integration
	print 'Flagging first integration...'
	inte=raw_input('Enter the integration time. e.g., 30.0-->')
	dict_log.append(('integration time',inte))
	flagdata(vis=ms_namel, mode='quack', quackinterval=float(inte), quackmode='beg',quackincrement=True)
	flagdata(vis=ms_nameu, mode='quack', quackinterval=float(inte), quackmode='beg',quackincrement=True)
	#dummy cont spw
	dum0=raw_input('Is the dummy cont spw present?y or n-->')
	dict_log.append(('dum_spw_flag',dum0))
	if dum0=='y':
		print 'Flagging dummy spw in both side-bands...'
		flagdata(vis=ms_namel,spw='0',field='',antenna='')
		flagdata(vis=ms_nameu,spw='0',field='',antenna='')
	#end channels
	flag_end=raw_input('Do you want to flag the end channels?y or n-->')
	dict_log.append(('flag_end_channels',flag_end))
	if flag_end=='y':
		beg=raw_input('Please enter beginning channels to flag. e.g., 0~3-->')
		endd=raw_input('Please enter end channels to flag. e.g., 60~63-->')
		dict_log.append(('beginning_channels',beg))
		dict_log.append(('end_channels',endd))
		print 'Flagging end channels...'
		flagdata(vis=ms_namel,spw=spw_low+':'+beg,field='',antenna='')
		flagdata(vis=ms_nameu,spw=spw_high+':'+endd,field='',antenna='')
	else:
		print 'End channels not flagged.'
	print ' (1) Plotting Amp vs time. Look for obvious bad data and ipointing data.'
	print 'LSB first...'
	plotms(vis=ms_namel,xaxis="time",yaxis="amp",coloraxis="field",iteraxis="antenna")
	raw_input('Please press enter when ready to continue.')
	print 'USB next...'
	plotms(vis=ms_nameu,xaxis="time",yaxis="amp",coloraxis="field",iteraxis="antenna")
	raw_input('Please press enter when ready to continue.')
	#ipointing data
	flagipoint=raw_input('Do you want to flag ipointing data?y or n-->')
	if flagipoint=='y':
		dict_log.append(('flag_ipoint',flagipoint))
		print 'Flagging ipointing data...'
		point_lsbf,point_lsbt=raw_input('Please enter field and timerange of ipointing in lsb. e.g. 2 2015/06/17/17:08:44.0~2015/06/17/17:14:30-->').split(' ')
		point_usbf,point_usbt=raw_input('Please enter field and timerange of ipointing in usb. e.g. 2 2015/06/17/17:08:44.0~2015/06/17/17:14:30-->').split(' ')
		dict_log.append(('ipointl_field',point_lsbf))
		dict_log.append(('ipointl_time',point_lsbt))
		dict_log.append(('ipointu_field',point_usbf))
		dict_log.append(('ipointu_time',point_usbt))
		flagdata(ms_namel,field=point_lsbf,timerange=point_lsbt)
		flagdata(ms_nameu,field=point_usbf,timerange=point_usbt)
	else:
		print 'No ipointing data to flag.'
		dict_log.append(('ipointl_field','none'))
		dict_log.append(('ipointl_time','none'))
		dict_log.append(('ipointu_field','none'))
		dict_log.append(('ipointu_time','none'))
	print '(2) Baselines versus antenna to look for bad antennas,bad channels/spws...'
	lastf,fant=raw_input('Please enter last field id and first antenna. e.g., 2 0-->').split(' ')
	dict_log.append(('last_field',lastf))
	dict_log.append(('first_antenna',fant))
	print 'LSB...'
	plotms(vis=ms_namel,field=lastf,spw='', antenna=fant,xaxis='antenna2',yaxis='amp')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotms(vis=ms_nameu,field=lastf,spw='', antenna=fant,xaxis='antenna2',yaxis='amp')
	raw_input('Please press enter when ready to continue.')
	print '(3) Bandpass for reference antenna (second_cal)...'
	print 'Look for bad spws and RFI spikes'
	print 'LSB...'
	plotms(vis=ms_namel,field=second_cal,spw='',antenna=ref_ant+'&*',xaxis='frequency', \
		yaxis='amp',iteraxis='baseline',coloraxis='spw')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotms(vis=ms_nameu,field=second_cal,spw='',antenna=ref_ant+'&*',xaxis='frequency', \
		yaxis='amp',iteraxis='baseline',coloraxis='spw')
	raw_input('Please press enter when ready to continue.')
	print '(4) Data stream plot...'
	print 'LSB...'
	plotms(vis=ms_namel,field='',timerange='',antenna=fant,spw='',xaxis='time',yaxis='antenna2',\
		plotrange=[-1,-1,0,8],coloraxis='field')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotms(vis=ms_nameu,field='',timerange='',antenna=fant,spw='',xaxis='time',yaxis='antenna2',\
		plotrange=[-1,-1,0,8],coloraxis='field')
	raw_input('Please press enter when ready to continue.')
	print 'Flagging...'
	#from log
	##bad ants
	flag_badant=raw_input('Do you want to flag any known bad antennas?y or n-->')
	dict_log.append(('flag_badant',flag_badant))
	if flag_badant=='y':
		badal=raw_input('Please enter bad ants to flag in lsb. e.g. 2,3-->')
		badau=raw_input('Please enter bad ants to flag in usb. e.g. 2,3-->')
		dict_log.append(('badant_lsb',badal))
		dict_log.append(('badant_usb',badau))
		flagdata(vis=ms_namel, flagbackup=True, mode='manual', antenna=badal)
		flagdata(vis=ms_nameu, flagbackup=True, mode='manual', antenna=badau)
	raw_input('Please press enter when ready to continue.')
	##bad scans
	flag_badscan=raw_input('Do you want to flag any known bad scans?y or n-->')
	dict_log.append(('flag_badscan',flag_badscan))
	if flag_badscan=='y':
		badsl=raw_input('Please enter bad scans to flag in lsb. e.g. 2,3-->')
		badsu=raw_input('Please enter bad scans to flag in usb. e.g. 2,3-->')
		dict_log.append(('badscan_lsb',badsl))
		dict_log.append(('badscan_usb',badsu))
		flagdata(vis=ms_namel, flagbackup=True, mode='manual', scan=badsl)
		flagdata(vis=ms_nameu, flagbackup=True, mode='manual', scan=badsu)
	raw_input('Please press enter when ready to continue.')
	#phase jumps; ant,field,timerange
	flag_pj=raw_input('Do you want to flag any known phase jumps?y or n-->')
	dict_log.append(('flag_phasejump',flag_pj))
	if flag_pj=='y':
		badpjl=raw_input('Please enter bad ants,fields,timeranges to flag in lsb. e.g. 2;3;2015/07/02/12:15:00.0~2015/07/02/12:30:00.0-->').split(' ')
		badpju=raw_input('Please enter bad ants,fields,timeranges to flag in usb. e.g. 2;3;2015/07/02/12:15:00.0~2015/07/02/12:30:00.0-->').split(' ')
		dict_log.append(('bad_phasejumpl',badpjl))
		dict_log.append(('bad_phasejumpl',badpju))
		for i in range(0,len(badpjl)):
			strg_pjl=badpjl[i].split(';')
			flagdata(vis=ms_namel, mode='manual', antenna=strg_pjl[0],field=strg_pjl[1],timerange=strg_pjl[2])
		for i in range(0,len(badpju)):
			strg_pju=badpju[i].split(';')
			flagdata(vis=ms_nameu, mode='manual', antenna=strg_pju[0],field=strg_pju[1],timerange=strg_pju[2])
		raw_input('Please press enter when ready to continue.')
	#other bad data ant,spw,field
	badasflsb=raw_input('Please enter bad ant,spw,and field to flag in lsb(enter if none). e.g., 2,3;5:4~9;3 ;5;3-->').split(' ')
	badasfusb=raw_input('Please enter bad ant,spw,and field to flag in usb (enter if none). e.g., 2,3;5:4~9;3 ;5;3-->').split(' ')
	dict_log.append(('bad_antspefield_lsb',badasflsb))
	dict_log.append(('bad_antspefield_usb',badasfusb))
	if '' in badasflsb:
		print 'Nothing to flag.'
	else:
		print 'Flagging selected lsb data.'
		for i in range(0,len(badasflsb)):
			strglsb=badasflsb[i].split(';')
			flagdata(vis=ms_namel,flagbackup=True, mode='manual', antenna=strglsb[0],spw=strglsb[1],field=strglsb[2])
	if '' in badasfusb:
		print 'Nothing to flag.'
	else:
		print 'Flagging selected usb data.'
		for i in range(0,len(badasfusb)):
			strgusb=badasfusb[i].split(';')
			flagdata(vis=ms_nameu,flagbackup=True, mode='manual', antenna=strgusb[0],spw=strgusb[1],field=strgusb[2])
	print 'Final check of flagged data...'
	print 'LSB...'
	plotms(vis=ms_namel,field=second_cal,spw='', antenna=ref_ant,xaxis='frequency',yaxis='amp')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotms(vis=ms_nameu,field=second_cal,spw='', antenna=ref_ant,xaxis='frequency',yaxis='amp')
	raw_input('Please press enter when ready to continue.')
	flag_again=raw_input('Do you need to do more flagging? y or n-->')
	while flag_again=='y':
		count_f=1
		badasflsb=raw_input('Please enter bad ant,spw,field,scan/timerange to flag in lsb (enter if none). e.g., 2,3;5:4~9;3;10:48:00~10:56:00 ;5;3;4,5-->').split(' ')
		badasfusb=raw_input('Please enter bad ant,spw,field,scan/timerange to flag in usb (enter if none). e.g., 2,3;5:4~9;3;10:48:00~10:56:00 ;5;3;4,5-->').split(' ')
		dict_log.append((ms_namel_prefix+'_flag_antspwfield_lsb'+str(count_f),badasflsb))
		dict_log.append((ms_nameu_prefix+'_flag_antspwfield_usb'+str(count_f),badasfusb))
		if '' in badasflsb:
			print 'Nothing to flag.'
		else:
			print 'Flagging selected lsb data.'
			for i in range(0,len(badasflsb)):
				strglsb=badasflsb[i].split(';')
				if ':' in strglsb[3]:
					flagdata(vis=ms_namel,flagbackup=True, mode='manual', antenna=strglsb[0],spw=strglsb[1],field=strglsb[2],timerange=strglsb[3])
				else:
					flagdata(vis=ms_namel,flagbackup=True, mode='manual', antenna=strglsb[0],spw=strglsb[1],field=strglsb[2],scan=strglsb[3])
		if '' in badasfusb:
			print 'Nothing to flag.'
		else:
			print 'Flagging selected usb data.'
			for i in range(0,len(badasfusb)):
				strgusb=badasfusb[i].split(';')
				if ':' in strgusb[3]:
					flagdata(vis=ms_nameu,flagbackup=True, mode='manual', antenna=strgusb[0],spw=strgusb[1],field=strgusb[2],timerange=strgusb[3])
				else:
					flagdata(vis=ms_nameu,flagbackup=True, mode='manual', antenna=strgusb[0],spw=strgusb[1],field=strgusb[2],scan=strgusb[3])
		print 'Plotting...'
		print 'LSB...'
		plotms(vis=ms_namel,field=second_cal,spw='', antenna=ref_ant,xaxis='frequency',yaxis='amp')
		raw_input('Please press enter when ready to continue.')
		print 'USB...'
		plotms(vis=ms_nameu,field=second_cal,spw='', antenna=ref_ant,xaxis='frequency',yaxis='amp')
		raw_input('Please press enter when ready to continue.')
		count_f=count_f+1
		flag_again=raw_input('Do you need to do more flagging? y or n-->')
	intera=raw_input('Flagging is finished. Do you want to do interactive calibration?y or n-->')
	dict_log.append(('interactive',intera))
	if intera=='n':
		print 'You have chosen to not do interactive calibration.'
		print 'No plots will be made and no additional flagging is required.'
		print 'Please go do something else for a while.'
		print 'You will not be prompted until the check for latent baseline issues after bandpass.'
		print '*********************************************************************'
	else:
		print 'The rest of the script is interactive. Please stay by the computer.'
		print '*********************************************************************'
else:
	intera=raw_input('Flagging is finished. Do you want to do interactive calibration?y or n-->')
	dict_log.append(('interactive',intera))
	if intera=='n':
		print 'You have chosen to not do interactive calibration.'
		print 'No plots will be made and no additional flagging is required.'
		print 'Please go do something else for a while.'
		print 'You will not be prompted until the check for latent baseline issues after bandpass.'
		print '*********************************************************************'
	else:
		print 'The rest of the script is interactive. Please stay by the computer.'
		print '*********************************************************************'
########################################

#####################################
#Antenna position corrections
#####################################
if do_ant_correct=='T':
	print 'Performing antenna position corrections...'
	ant_pos_file_read=ascii.read(ant_corr_file,delimiter=' ',data_start=0,names=['ant_lst','X','Y','Z'],guess=False)
	ant_lst=",".join(ant_pos_file_read['ant_lst'])
	ant_lst_offset=[]
	for i in range(0,len(ant_pos_file_read['ant_lst'])):
		ant_lst_offset.append(1e-3*tt['X'][i])
		ant_lst_offset.append(1e-3*tt['Y'][i])
		ant_lst_offset.append(1e-3*tt['Z'][i])
	print 'Selected antennas for corrections: ',ant_lst
	print 'Offsets in meters (X,Y,Z looped over all antennas): ', ant_lst_offset
	gencal(caltable=cal_table_prefixl+'antpos',vis=ms_namel, caltype='antpos',parameter=ant_lst_offset,antenna=ant_lst)
	gencal(caltable=cal_table_prefixu+'antpos',vis=ms_nameu, caltype='antpos',parameter=ant_lst_offset,antenna=ant_lst)
	applycal(vis=ms_namel,gaintable=cal_table_prefixl+'antpos',gainfield='')
	applycal(vis=ms_nameu,gaintable=cal_table_prefixu+'antpos',gainfield='')
#####################################

#####################################
#Bandpass and Phase solutions
#####################################
#initial phase solution to take away major phase variations with time
os.system('rm -rf '+cal_table_prefixl+'phase_int.cal')
os.system('rm -rf '+cal_table_prefixu+'phase_int.cal')
print 'Initial phase solution before bandpass...'
gaincal(vis=ms_namel,caltable=cal_table_prefixl+'phase_int.cal',field=bpf_cal,solint="int",\
	calmode='p',refant=ref_ant,gaintype="G",minsnr=2.0,spw=spw_low,combine="spw")
gaincal(vis=ms_nameu,caltable=cal_table_prefixu+'phase_int.cal',field=bpf_cal,solint="int",\
	calmode="p",refant=ref_ant,gaintype="G",minsnr=2.0,spw=spw_high,combine="spw")
if intera=='y':
	print 'Plotting solutions...'
	print 'LSB...'
	plotcal(caltable=cal_table_prefixl+'phase_int.cal',xaxis="time",yaxis="phase",subplot=331,\
		iteration="antenna",plotrange=[0,0,-180,180],overplot=False,clearpanel='Auto',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'phase_intl.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB ...'
	plotcal(caltable=cal_table_prefixu+'phase_int.cal',xaxis="time",yaxis="phase",subplot=331,\
		iteration="antenna",plotrange=[0,0,-180,180],overplot=False,clearpanel='Auto',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'phase_intu.png')
	raw_input('Please press enter when ready to continue.')
#now bandpass solution
os.system('rm -rf '+cal_table_prefixl+'bandpassB.cal')
os.system('rm -rf '+cal_table_prefixu+'bandpassB.cal')
print 'Bandpass solution...'
bandpass(vis=ms_namel, caltable=cal_table_prefixl+'bandpassB.cal', field=bpf_cal,\
	spw=spw_low, solint='inf', combine='scan',bandtype='B', refant=ref_ant,solnorm=True,\
	 gaintable=cal_table_prefixl+'phase_int.cal',spwmap=[first_lsb_spw])
bandpass(vis=ms_nameu, caltable=cal_table_prefixu+'bandpassB.cal', field=bpf_cal,\
	spw=spw_high, solint='inf', combine='scan',bandtype='B', refant=ref_ant,solnorm=True,\
	 gaintable=cal_table_prefixu+'phase_int.cal',spwmap=[first_usb_spw])
if intera=='y':
	print 'Plotting solutions...'
	print 'LSB...'
	plotcal(caltable=cal_table_prefixl+'bandpassB.cal',xaxis='freq',yaxis='amp',\
		field=bpf_cal,antenna='',spw='',timerange='',subplot=311,\
		overplot=False,clearpanel='Auto',iteration='antenna',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'bandpassBl.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotcal(caltable=cal_table_prefixu+'bandpassB.cal',xaxis='freq',yaxis='amp',\
		field=bpf_cal,antenna='',spw='',timerange='',subplot=311,\
		overplot=False,clearpanel='Auto',iteration='antenna',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'bandpassBu.png')
	raw_input('Please press enter when ready to continue.')
#redo phase after applying bandpass, esspecially good if had variations in intial one
os.system('rm -rf '+cal_table_prefixl+"phase_all.cal")
os.system('rm -rf '+cal_table_prefixu+"phase_all.cal")
print 'Redo phase solution after bandpass...'
gaincal(vis=ms_namel,caltable=cal_table_prefixl+"phase_all.cal",field=bpf_cal,solint="int",calmode="p",\
	refant=ref_ant,gaintype="G",minsnr=2.0,spw=spw_low,combine="spw",gaintable=cal_table_prefixl+"bandpassB.cal")
gaincal(vis=ms_nameu,caltable=cal_table_prefixu+"phase_all.cal",field=bpf_cal,solint="int",calmode="p",\
	refant=ref_ant,gaintype="G",minsnr=2.0,spw=spw_high,combine="spw",gaintable=cal_table_prefixu+"bandpassB.cal")
if intera=='y':
	print 'Plotting solutions...'
	print 'LSB...'
	plotcal(caltable=cal_table_prefixl+'phase_all.cal',xaxis='time',yaxis='phase',field=bpf_cal,\
		antenna='',spw='',timerange='',subplot=311,overplot=False,clearpanel='Auto',iteration='antenna',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'phase_alll.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotcal(caltable=cal_table_prefixu+'phase_all.cal',xaxis='time',yaxis='phase',field=bpf_cal,\
		antenna='',spw='',timerange='',subplot=311,overplot=False,clearpanel='Auto',iteration='antenna',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'phase_all.png')
	raw_input('Please press enter when ready to continue.')
#####################################

#####################################
##Amp solutions
#####################################
#solve for amp while applying bandpass and phase
os.system('rm -rf '+cal_table_prefixl+'ap.cal')
os.system('rm -rf '+cal_table_prefixu+'ap.cal')
print 'Amp solutions...'
gaincal(vis=ms_namel, caltable=cal_table_prefixl+'ap.cal',field=bpf_cal, spw=spw_low,\
 gaintype='G', minsnr=2.0,refant=ref_ant, calmode='a',solint='300s', combine='spw',\
 gaintable=[cal_table_prefixl+"phase_all.cal",cal_table_prefixl+"bandpassB.cal"],spwmap=[[first_lsb_spw],[]])
gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'ap.cal',field=bpf_cal, spw=spw_high,\
 gaintype='G', minsnr=2.0,refant=ref_ant, calmode='a',solint='300s', combine='spw',\
 gaintable=[cal_table_prefixu+"phase_all.cal",cal_table_prefixu+"bandpassB.cal"],spwmap=[[first_usb_spw],[]])
if intera=='y':
	print 'Plotting solutions...'
	print 'LSB...'
	plotcal(caltable=cal_table_prefixl+'ap.cal',xaxis='time',yaxis='amp',field=bpf_cal,antenna='',\
		spw='',timerange='',subplot=331,overplot=False,clearpanel='Auto',iteration='antennas',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'apl.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotcal(caltable=cal_table_prefixu+'ap.cal',xaxis='time',yaxis='amp',field=bpf_cal,antenna='',\
		spw='',timerange='',subplot=331,overplot=False,clearpanel='Auto',iteration='antennas',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'apu.png')
	raw_input('Please press enter when ready to continue.')
#####################################

#####################################
#apply bandpass,phase,amp to bp cal
#####################################
print 'Applying Bandpass, phase and amp solutions to bp cal...'
applycal(vis=ms_namel,spw=spw_low, field=bpf_cal,
	gaintable=[cal_table_prefixl+"phase_all.cal",cal_table_prefixl+'ap.cal',
	cal_table_prefixl+'bandpassB.cal'],
	spwmap=[[first_lsb_spw],[first_lsb_spw],[]],gainfield=[bpf_cal,bpf_cal,bpf_cal])
applycal(vis=ms_nameu,spw=spw_high, field=bpf_cal,
	gaintable=[cal_table_prefixu+"phase_all.cal",cal_table_prefixu+'ap.cal',
	cal_table_prefixu+'bandpassB.cal'],
	spwmap=[[first_usb_spw],[first_usb_spw],[]],gainfield=[bpf_cal,bpf_cal,bpf_cal])
#####################################

#####################################
#Check for latent baseline issues
#####################################
#check no latent baseline issues, that need to be calibrated out
#have to do blcal if you see them
if intera=='n':
	print '**********************************************'
	print 'Entering interactive part of script again...'
	print '**********************************************'
print 'Checking for latent baseline issues...'
print 'Plotting LSB...'
plotms(vis=ms_namel,xaxis="freq",yaxis="amp",
       ydatacolumn="corrected",iteraxis="baseline",
       spw=spw_low,field=bpf_cal,averagedata=True,avgtime='99999')
raw_input('Please press enter when ready to continue.')
print 'Plotting USB...'
plotms(vis=ms_nameu,xaxis="freq",yaxis="amp",
       ydatacolumn="corrected",iteraxis="baseline",
       spw=spw_high,field=bpf_cal,averagedata=True,avgtime='99999')
raw_input('Please press enter when ready to continue.')
doblcal=raw_input('Are there latent baseline issues?y or n-->')
dict_log.append(('latent_baseline',doblcal))
if doblcal=='y':
	os.system('rm -rf '+cal_table_prefixl+'blcal1.cal')
	os.system('rm -rf '+cal_table_prefixl+'blcal2.cal')
	os.system('rm -rf '+cal_table_prefixu+'blcal1.cal')
	os.system('rm -rf '+cal_table_prefixu+'blcal2.cal')
	print 'If you flagged any full spws they need to be selected out here.'
	print 'LSB...'
	spwflaglsb=raw_input('Please enter spw range not including those fully flagged in lsb. e.g., 1~22,24~50-->')
	dict_log.append(('unflagged_spwlsb',spwflaglsb))
	blcal(vis=ms_namel, caltable=cal_table_prefixl+'blcal1.cal',field=bpf_cal,\
	 spw=spwflaglsb,solint='inf', combine='spw,scan',\
	 gaintable=[cal_table_prefixl+"phase_all.cal",cal_table_prefixl+'ap.cal',\
	 cal_table_prefixl+'bandpassB.cal'], calmode='a',spwmap=[[first_lsb_spw],[first_lsb_spw],[]])
	blcal(vis=ms_namel, caltable=cal_table_prefixl+'blcal2.cal',field=bpf_cal,\
	 spw=spwflaglsb,solint='inf', combine='scan',\
	 gaintable=[cal_table_prefixl+"phase_all.cal",cal_table_prefixl+'ap.cal',\
	 cal_table_prefixl+'blcal1.cal',cal_table_prefixl+'bandpassB.cal'], calmode='a',spwmap=[[first_lsb_spw],[first_lsb_spw],[first_lsb_spw],[]])
	print 'USB...'
	spwflagusb=raw_input('Please enter spw range not including those fully flagged in usb. e.g., 1~22,24~50-->')
	dict_log.append(('unflagged_spwusb',spwflagusb))
	blcal(vis=ms_nameu, caltable=cal_table_prefixu+'blcal1.cal',field=bpf_cal,\
		spw=spwflagusb,solint='inf', combine='spw,scan',\
		gaintable=[cal_table_prefixu+"phase_all.cal",cal_table_prefixu+'ap.cal',\
		cal_table_prefixu+'bandpassB.cal'], calmode='a',spwmap=[[first_usb_spw],[first_usb_spw],[]])
	blcal(vis=ms_nameu, caltable=cal_table_prefixu+'blcal2.cal',field=bpf_cal,\
		spw=spwflagusb,solint='inf', combine='scan',\
		gaintable=[cal_table_prefixu+"phase_all.cal",cal_table_prefixu+'ap.cal',\
		cal_table_prefixu+'blcal1.cal',cal_table_prefixu+'bandpassB.cal'], calmode='a',spwmap=[[first_usb_spw],[first_usb_spw],[first_usb_spw],[]])
	print 'Plotting solutions...'
	print 'Should be close to 1.'
	print 'LSB...'
	plotcal(caltable=cal_table_prefixl+'blcal2.cal',xaxis='freq',yaxis='amp',field=bpf_cal,\
		antenna='',spw='',timerange='',subplot=511,overplot=False,clearpanel='Auto',\
		iteration='antenna',plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'blcal2l.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotcal(caltable=cal_table_prefixu+'blcal2.cal',xaxis='freq',yaxis='amp',field=bpf_cal,\
		antenna='',spw='',timerange='',subplot=511,overplot=False,clearpanel='Auto',\
		iteration='antenna',plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'blcal2u.png')
	raw_input('Please press enter when ready to continue.')
	#reapply solutions to bp cal
	print 'Reapplying solutions to BP cal...'
	applycal(vis=ms_namel,spw=spw_low, field=bpf_cal,gaintable=[cal_table_prefixl+"phase_all.cal",\
		cal_table_prefixl+'ap.cal',cal_table_prefixl+'blcal2.cal',cal_table_prefixl+'bandpassB.cal'],\
		spwmap=[[first_lsb_spw],[first_lsb_spw],[],[]],gainfield=[bpf_cal,bpf_cal,bpf_cal,bpf_cal])
	applycal(vis=ms_nameu,spw=spw_high, field=bpf_cal,gaintable=[cal_table_prefixu+"phase_all.cal",\
		cal_table_prefixu+'ap.cal',cal_table_prefixu+'blcal2.cal',cal_table_prefixu+'bandpassB.cal'],\
		spwmap=[[first_usb_spw],[first_usb_spw],[],[]],gainfield=[bpf_cal,bpf_cal,bpf_cal,bpf_cal])
	print 'Plotting solutions...'
	print 'LSB...'
	plotms(vis=ms_namel,xaxis="freq",yaxis="amp",ydatacolumn="corrected",iteraxis="baseline",spw=spw_low,\
		field=bpf_cal,averagedata=True,avgtime='99999')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotms(vis=ms_nameu,xaxis="freq",yaxis="amp",ydatacolumn="corrected",iteraxis="baseline",spw=spw_high,\
		field=bpf_cal,averagedata=True,avgtime='99999')
	raw_input('Please press enter when ready to continue.')
#####################################

#####################################
#Apply bandpass sols to other cals
#####################################
print 'Applying bandpass solutions to other cals...'
if doblcal=='y':
	applycal(vis=ms_namel,spw=spw_low, field=second_cal+','+flux_cal,\
		gaintable=[cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],\
		spwmap=[[],[]],gainfield=[bpf_cal,bpf_cal])
	applycal(vis=ms_nameu,spw=spw_high, field=second_cal+','+flux_cal,\
		gaintable=[cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
		spwmap=[[],[]],gainfield=[bpf_cal,bpf_cal])
else:
	applycal(vis=ms_namel,spw=spw_low, field=second_cal+','+flux_cal,\
		gaintable=[cal_table_prefixl+'bandpassB.cal'],\
		spwmap=[[]],gainfield=[bpf_cal])
	applycal(vis=ms_nameu,spw=spw_high, field=second_cal+','+flux_cal,\
		gaintable=[cal_table_prefixu+'bandpassB.cal'],\
		spwmap=[[]],gainfield=[bpf_cal])
#####################################

#####################################
#Set fluxscale
#####################################
print 'Setting flux scale...'
for i in range(0,len(flux_lst)):
	setjy(vis=ms_namel, field=flux_lst[i], spw=spw_low,standard='Butler-JPL-Horizons 2012')
	setjy(vis=ms_nameu, field=flux_lst[i], spw=spw_high,standard='Butler-JPL-Horizons 2012')
#####################################

#####################################
#final gain and flux calibration
#####################################
#final gain and flux calibration for all cals
#solve for 2 soln, short and long, short fed to amp step, long for target
os.system('rm -rf '+cal_table_prefixl+'allp.cal')
os.system('rm -rf '+cal_table_prefixl+'allpscan.cal')
os.system('rm -rf '+cal_table_prefixu+'allp.cal')
os.system('rm -rf '+cal_table_prefixu+'allpscan.cal')
print 'Final gain and flux calibration...'
print 'Solving for 2 solutions. Short for amp cal, long for target.'
print 'Starting with the short solution...'
if doblcal=='y':
	gaincal(vis=ms_namel, caltable=cal_table_prefixl+'allp.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spwflaglsb, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='int', combine='spw',\
		gaintable=[cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],spwmap=[[],[]])
	gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allp.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spwflagusb, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='int', combine='spw',\
		gaintable=[cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],spwmap=[[],[]])
else:
	gaincal(vis=ms_namel, caltable=cal_table_prefixl+'allp.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spw_low, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='int', combine='spw',\
		gaintable=[cal_table_prefixl+'bandpassB.cal'],spwmap=[[]])
	gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allp.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spw_high, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='int', combine='spw',\
		gaintable=[cal_table_prefixu+'bandpassB.cal'],spwmap=[[]])
print 'Now the longer solutions...'
if doblcal=='y':
	gaincal(vis=ms_namel, caltable=cal_table_prefixl+'allpscan.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
	 spw=spwflaglsb, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='inf', combine='spw',\
	 gaintable=[cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],spwmap=[[],[]])
	gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allpscan.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
	 spw=spwflagusb, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='inf', combine='spw',\
	 gaintable=[cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],spwmap=[[],[]])
else:
	gaincal(vis=ms_namel, caltable=cal_table_prefixl+'allpscan.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spw_low, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='inf', combine='spw',\
		gaintable=[cal_table_prefixl+'bandpassB.cal'],spwmap=[[]])
	gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allpscan.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spw_high, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='p',solint='inf', combine='spw',\
		gaintable=[cal_table_prefixu+'bandpassB.cal'],spwmap=[[]])
if intera=='y':
	print 'Plotting solutions...'
	print 'LSB short solution...'
	plotcal(caltable=cal_table_prefixl+'allp.cal',xaxis='time',yaxis='phase',\
		field=flux_cal+','+bpf_cal+','+second_cal,antenna='',spw='',timerange='',\
		subplot=511,overplot=False,clearpanel='Auto',iteration='antenna,field',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'allpl.png')
	raw_input('Please press enter when ready to continue.')
	print 'LSB long solution...'
	plotcal(caltable=cal_table_prefixl+'allpscan.cal',xaxis='time',yaxis='phase',\
		field=flux_cal+','+bpf_cal+','+second_cal,antenna='',spw='',timerange='',\
		subplot=511,overplot=False,clearpanel='Auto',iteration='antenna,field',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'allpscanl.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB short solution...'
	plotcal(caltable=cal_table_prefixu+'allp.cal',xaxis='time',yaxis='phase',\
		field=flux_cal+','+bpf_cal+','+second_cal,antenna='',spw='',timerange='',\
		subplot=511,overplot=False,clearpanel='Auto',iteration='antenna,field',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'allpu.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB long solution...'
	plotcal(caltable=cal_table_prefixu+'allpscan.cal',xaxis='time',yaxis='phase',\
		field=flux_cal+','+bpf_cal+','+second_cal,antenna='',spw='',timerange='',\
		subplot=511,overplot=False,clearpanel='Auto',iteration='antenna',\
		plotsymbol='o',plotcolor='blue',showgui=True,figfile=my_dir+'allpscanu.png')
	raw_input('Please press enter when ready to continue.')
#amp with applying short phase
os.system('rm -rf '+cal_table_prefixl+'allap.cal')
os.system('rm -rf '+cal_table_prefixu+'allap.cal')
print 'Amp cal while applying short phase...'
doif2=raw_input('Do you want to only use IF2 for flux cal in the USB?y or n-->')
dict_log.append(('only_if2_flux_usb',doif2))

if doblcal=='y':
	gaincal(vis=ms_namel, caltable=cal_table_prefixl+'allap.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
		spw=spwflaglsb, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='ap',solint='300s', combine='spw',\
		gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],\
		spwmap=[[first_lsb_spw],[],[]])
	if doif2=='y':
		doif2sp=raw_input('Enter spw range of IF2. e.g., 25~50-->')
		dict_log.append(('if2_spw_usb',doif2sp))
		gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allap.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
			spw=doif2sp, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='ap',solint='300s', combine='spw',\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
			spwmap=[[first_usb_spw],[],[]])
	else:
		gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allap.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
			spw=spwflagusb, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='ap',solint='300s', combine='spw',\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
			spwmap=[[first_usb_spw],[],[]])

else:
	gaincal(vis=ms_namel, caltable=cal_table_prefixl+'allap.cal',field=flux_cal+','+bpf_cal+','+second_cal,
		spw=spw_low, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='ap',solint='300s', combine='spw',\
		gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'bandpassB.cal'],\
		spwmap=[[first_lsb_spw],[]])
	if doif2=='y':
		doif2sp=raw_input('Enter spw range of IF2. e.g., 25~50-->')
		dict_log.append(('if2_spw_usb',doif2sp))
		gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allap.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
			spw=doif2sp, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='ap',solint='300s', combine='spw',\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'bandpassB.cal'],\
			spwmap=[[first_usb_spw],[]])
	else:
		gaincal(vis=ms_nameu, caltable=cal_table_prefixu+'allap.cal',field=flux_cal+','+bpf_cal+','+second_cal,\
			spw=spw_high, gaintype='G', minsnr=2.0,refant=ref_ant, calmode='ap',solint='300s', combine='spw',\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'bandpassB.cal'],\
			spwmap=[[first_usb_spw],[]])
if intera=='y':
	print 'Plotting solutions...'
	print 'LSB...'
	plotcal(caltable=cal_table_prefixl+'allap.cal',xaxis='time',yaxis='amp',field=flux_cal+','+bpf_cal+','+second_cal,\
		antenna='',spw='',timerange='',subplot=331,overplot=False,clearpanel='Auto',iteration='antennas',plotsymbol='o',\
		plotcolor='blue',showgui=True,figfile=my_dir+'allapl.png')
	raw_input('Please press enter when ready to continue.')
	print 'USB...'
	plotcal(caltable=cal_table_prefixu+'allap.cal',xaxis='time',yaxis='amp',field=flux_cal+','+bpf_cal+','+second_cal,\
		antenna='',spw='',timerange='',subplot=331,overplot=False,clearpanel='Auto',iteration='antennas',plotsymbol='o',\
		plotcolor='blue',showgui=True,figfile=my_dir+'allapu.png')
	raw_input('Please press enter when ready to continue.')
#####################################

#####################################
#Derive absolute flux scale
#####################################
#derive absolute flux scale with setjy above
os.system('rm -rf '+cal_table_prefixl+'flux.cal')
os.system('rm -rf '+cal_table_prefixu+'flux.cal')
if2startspw=int(doif2sp.split('~')[0])
dict_log.append(('if2_spw_start',if2startspw))
print 'Deriving absolute flux scale...'
fluxscale(vis=ms_namel,caltable=cal_table_prefixl+'allap.cal',\
	fluxtable=cal_table_prefixl+'flux.cal',reference=flux_cal,\
	transfer=bpf_cal+','+second_cal)
if doif2=='y':
	fluxscale(vis=ms_nameu,caltable=cal_table_prefixu+'allap.cal',\
		fluxtable=cal_table_prefixu+'flux.cal',reference=flux_cal,\
		transfer=bpf_cal+','+second_cal,refspwmap=[if2startspw])
else:
	fluxscale(vis=ms_nameu,caltable=cal_table_prefixu+'allap.cal',\
		fluxtable=cal_table_prefixu+'flux.cal',reference=flux_cal,\
		transfer=bpf_cal+','+second_cal)
#####################################

#####################################
#Apply Calibration
#####################################
print 'Applying calbration...'
print 'BP cal(s)...'
for i in range(0,len(bpf_lst)):
	if doblcal=='y':
		applycal(vis=ms_namel,spw=spw_low, field=bpf_lst[i],\
			gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
			cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],\
			spwmap=[[first_lsb_spw],[first_lsb_spw],[],[]],gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i],bpf_lst[i]])
		if doif2=='y':
			applycal(vis=ms_nameu,spw=spw_high, field=bpf_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
				spwmap=[[first_usb_spw],[if2startspw],[],[]],gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i],bpf_lst[i]])
		else:
			applycal(vis=ms_nameu,spw=spw_high, field=bpf_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
				spwmap=[[first_usb_spw],[first_usb_spw],[],[]],gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i],bpf_lst[i]])
	else:
		applycal(vis=ms_namel,spw=spw_low, field=bpf_lst[i],\
			gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
			cal_table_prefixl+'bandpassB.cal'],\
			spwmap=[[first_lsb_spw],[first_lsb_spw],[]],gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i]])
		if doif2=='y':
			applycal(vis=ms_nameu,spw=spw_high, field=bpf_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal'],\
				spwmap=[[first_usb_spw],[if2startspw],[]],gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i]])
		else:
			applycal(vis=ms_nameu,spw=spw_high, field=bpf_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal'],\
				spwmap=[[first_usb_spw],[first_usb_spw],[]],gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i]])
print 'Flux cal(s)...'
for i in range(0,len(flux_lst)):
	if doblcal=='y':
		applycal(vis=ms_namel,spw=spw_low, field=flux_lst[i],\
			gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
			cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],\
			spwmap=[[first_lsb_spw],[first_lsb_spw],[],[]],gainfield=[flux_lst[i],flux_lst[i],bpf_cal,bpf_cal])
		if doif2=='y':
			applycal(vis=ms_nameu,spw=spw_high, field=flux_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
				spwmap=[[first_usb_spw],[if2startspw],[],[]],gainfield=[flux_lst[i],flux_lst[i],bpf_cal,bpf_cal])
		else:
			applycal(vis=ms_nameu,spw=spw_high, field=flux_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
				spwmap=[[first_usb_spw],[first_usb_spw],[],[]],gainfield=[flux_lst[i],flux_lst[i],bpf_cal,bpf_cal])
	else:
		applycal(vis=ms_namel,spw=spw_low, field=flux_lst[i],\
			gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
			cal_table_prefixl+'bandpassB.cal'],\
			spwmap=[[first_lsb_spw],[first_lsb_spw],[]],gainfield=[flux_lst[i],flux_lst[i],bpf_cal])
		if doif2=='y':
			applycal(vis=ms_nameu,spw=spw_high, field=flux_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal'],\
				spwmap=[[first_usb_spw],[if2startspw],[]],gainfield=[flux_lst[i],flux_lst[i],bpf_cal])
		else:
			applycal(vis=ms_nameu,spw=spw_high, field=flux_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal'],\
				spwmap=[[first_usb_spw],[first_usb_spw],[]],gainfield=[flux_lst[i],flux_lst[i],bpf_cal])
print 'Second cal(s)...'
for i in range(0,len(second_lst)):
	if doblcal=='y':
		applycal(vis=ms_namel,spw=spw_low, field=second_lst[i],\
			gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
			cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],\
			spwmap=[[first_lsb_spw],[first_lsb_spw],[],[]], gainfield=[second_lst[i],second_lst[i],bpf_cal,bpf_cal])
		if doif2=='y':
			applycal(vis=ms_nameu,spw=spw_high, field=second_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
				spwmap=[[first_usb_spw],[if2startspw],[],[]], gainfield=[second_lst[i],second_lst[i],bpf_cal,bpf_cal])
		else:
			applycal(vis=ms_nameu,spw=spw_high, field=second_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
				spwmap=[[first_usb_spw],[first_usb_spw],[],[]], gainfield=[second_lst[i],second_lst[i],bpf_cal,bpf_cal])
	else:
		applycal(vis=ms_namel,spw=spw_low, field=second_lst[i],\
			gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
			cal_table_prefixl+'bandpassB.cal'],\
			spwmap=[[first_lsb_spw],[first_lsb_spw],[]], gainfield=[second_lst[i],second_lst[i],bpf_cal])
		if doif2=='y':
			applycal(vis=ms_nameu,spw=spw_high, field=second_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal'],\
				spwmap=[[first_usb_spw],[if2startspw],[]], gainfield=[second_lst[i],second_lst[i],bpf_cal])
		else:
			applycal(vis=ms_nameu,spw=spw_high, field=second_lst[i],\
				gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
				cal_table_prefixu+'bandpassB.cal'],\
				spwmap=[[first_usb_spw],[first_usb_spw],[]], gainfield=[second_lst[i],second_lst[i],bpf_cal])
print 'Target...'
if doblcal=='y':
	applycal(vis=ms_namel,spw=spw_low, field=target_id,\
		gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
		cal_table_prefixl+'bandpassB.cal',cal_table_prefixl+'blcal2.cal'],\
		spwmap=[[first_lsb_spw],[first_lsb_spw],[],[]], gainfield=[second_cal,second_cal,bpf_cal,bpf_cal])
	if doif2=='y':
		applycal(vis=ms_nameu,spw=spw_high, field=target_id,\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
			cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
			spwmap=[[first_usb_spw],[if2startspw],[],[]], gainfield=[second_cal,second_cal,bpf_cal,bpf_cal])
	else:
		applycal(vis=ms_nameu,spw=spw_high, field=target_id,\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
			cal_table_prefixu+'bandpassB.cal',cal_table_prefixu+'blcal2.cal'],\
			spwmap=[[first_usb_spw],[first_usb_spw],[],[]], gainfield=[second_cal,second_cal,bpf_cal,bpf_cal])
else:
	applycal(vis=ms_namel,spw=spw_low, field=target_id,\
		gaintable=[cal_table_prefixl+'allp.cal',cal_table_prefixl+'flux.cal',\
		cal_table_prefixl+'bandpassB.cal'],\
		spwmap=[[first_lsb_spw],[first_lsb_spw],[]], gainfield=[second_cal,second_cal,bpf_cal])
	if doif2=='y':
		applycal(vis=ms_nameu,spw=spw_high, field=target_id,\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
			cal_table_prefixu+'bandpassB.cal'],\
			spwmap=[[first_usb_spw],[if2startspw],[]], gainfield=[second_cal,second_cal,bpf_cal])
	else:
		applycal(vis=ms_nameu,spw=spw_high, field=target_id,\
			gaintable=[cal_table_prefixu+'allp.cal',cal_table_prefixu+'flux.cal',\
			cal_table_prefixu+'bandpassB.cal'],\
			spwmap=[[first_usb_spw],[first_usb_spw],[]], gainfield=[second_cal,second_cal,bpf_cal])
#####################################

#####################################
##split out target
#####################################
os.system('rm -rf '+split_low)
os.system('rm -rf '+split_high)
os.system('rm -rf '+split_full)
print 'Splitting out target data for both side-bands...'
split(vis=ms_namel,outputvis=split_low,\
	datacolumn='corrected',field=target_id,antenna='',spw=spw_low)
split(vis=ms_nameu,outputvis=split_high,\
	datacolumn='corrected',field=target_id,antenna='',spw=spw_high)
print 'Making concatenated full LSB+USB MS...'
concat(vis=[split_low,split_high],concatvis=split_full)
##########################################
if doImage=='T':
	###########################################
	##Imaging
	###########################################
	dopscl='n'
	dopscu='n'
	dopscb='n'
	print 'Imaging...'
	if use_auto=='T':
		myimsize=set_imagesize(split_low,0,'0')
		mycell=set_cellsize(split_low,0)
	if 'L' in bandsIM:
		print 'Lower side-band...'
		if mymask=='':
			os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_low+'_clean1*')
			print 'Using interactive mode so you can make a mask...'
			print 'Cleaning...'
			clean(vis=split_low, imagename=my_dir+target+'_'+date+'_'+band_low+'_clean1',field='',spw='',interactive=True,\
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
				mode='mfs',niter=0,nterms=mynterms,multiscale=multiscale,robust=robust)
		else:
			os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_low+'_clean1*')
			print 'Cleaning...'
			clean(vis=split_low, imagename=my_dir+target+'_'+date+'_'+band_low+'_clean1',field='',mask=mymask,spw='',interactive=False,\
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
				mode='mfs',niter=myniter,nterms=mynterms,multiscale=multiscale,robust=robust)
		if mynterms>1:
			imagenl=my_dir+target+'_'+date+'_'+band_low+'_clean1.image.tt0'
		else:
			imagenl=my_dir+target+'_'+date+'_'+band_low+'_clean1.image'
		print 'Correcting for PB...'
		os.system('rm -rf '+imagenl+'.pbcor')
		os.system('rm -rf '+imagenl+'.pbcor.fits')
		immath(imagename=[imagenl,my_dir+target+'_'+date+'_'+band_low+'_clean1.flux'],\
			expr='IM0/IM1',outfile=imagenl+'.pbcor')
		print 'Making fits image...'
		exportfits(imagename=imagenl+'.pbcor',fitsimage=imagenl+'.pbcor.fits')
		imagenl=imagenl+'.pbcor'
		fluxl,errl,unitl,freql,errl_real=imfit_point(imagenl,my_dir,'I')
		print 'Lower side-band flux density of ',fluxl,' +/- ',errl, unitl
		print 'Local RMS in Lower sideband image is: ',errl_real,' Jy'
		dopscl=raw_input('Do you want to do phase selfcal?y or n-->')
		dict_log.append(('phself_lsb',dopscl))
		if dopscl=='y':
			selfcal_low,scim_low=phselfcal(split_low,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
		date,band_low,'y','',multiscale,robust,weighting)
			fluxl_sc,errl_sc,unitl_sc,freql_sc,errl_real_sc=imfit_point(scim_low,my_dir)
	if 'U' in bandsIM:
		print 'Upper side-band...'
		if use_auto=='T':
			myimsize=set_imagesize(split_high,0,'0')
			mycell=set_cellsize(split_high,0)
		if mymask=='':
			os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_high+'_clean1*')
			print 'Using interactive mode so you can make a mask...'
			print 'Cleaning...'
			clean(vis=split_high, imagename=my_dir+target+'_'+date+'_'+band_high+'_clean1',field='',spw='',interactive=True,\
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
				mode='mfs',niter=0,nterms=mynterms,multiscale=multiscale,robust=robust)
		else:
			os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_high+'_clean1*')
			print 'Cleaning...'
			clean(vis=split_high, imagename=my_dir+target+'_'+date+'_'+band_high+'_clean1',field='',mask=mymask,spw='',interactive=False,\
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
				mode='mfs',niter=myniter,nterms=mynterms,multiscale=multiscale,robust=robust)
		if mynterms>1:
			imagenu=my_dir+target+'_'+date+'_'+band_high+'_clean1.image.tt0'
		else:
			imagenu=my_dir+target+'_'+date+'_'+band_high+'_clean1.image'
		print 'Correcting for PB...'
		os.system('rm -rf '+imagenu+'.pbcor')
		os.system('rm -rf '+imagenu+'.pbcor.fits')
		immath(imagename=[imagenu,my_dir+target+'_'+date+'_'+band_high+'_clean1.flux'],\
			expr='IM0/IM1',outfile=imagenu+'.pbcor')
		print 'Making fits image...'
		exportfits(imagename=imagenu+'.pbcor',fitsimage=imagenu+'.pbcor.fits')
		imagenu=imagenu+'.pbcor'
		fluxu,erru,unitu,frequ,erru_real=imfit_point(imagenu,my_dir,'I')
		print 'Upper side-band flux density of ',fluxu,' +/- ',erru, unitu
		print 'Local RMS in upper sideband image is: ',erru_real,' Jy'
		dopscu=raw_input('Do you want to do phase selfcal?y or n-->')
		dict_log.append(('phself_usb',dopscu))
		if dopscu=='y':
			selfcal_high,scim_high=phselfcal(split_high,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
		date,band_low,'y','',multiscale,robust,weighting)
			fluxu_sc,erru_sc,unitu_sc,frequ_sc,erru_real_sc=imfit_point(scim_high,my_dir)
	if 'B' in bandsIM:
		print 'Combined side-band...'
		if use_auto=='T':
			myimsize=set_imagesize(split_full,0,'0')
			mycell=set_cellsize(split_full,0)
		if mymask=='':
			os.system('rm -rf '+my_dir+target+'_'+date+'_both_clean1*')
			print 'Using interactive mode so you can make a mask...'
			print 'Cleaning...'
			clean(vis=[split_low,split_high], imagename=my_dir+target+'_'+date+'_both_clean1',field='',spw='',interactive=True,\
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
				mode='mfs',niter=0,nterms=mynterms,multiscale=multiscale,robust=robust)
		else:
			os.system('rm -rf '+my_dir+target+'_'+date+'_both_clean1*')
			print 'Cleaning...'
			clean(vis=[split_low,split_high], imagename=my_dir+target+'_'+date+'_both_clean1',field='',mask=mymask,spw='',interactive=False,\
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
				mode='mfs',niter=myniter,nterms=mynterms,multiscale=multiscale,robust=robust)
		if mynterms>1:
			imagenb=my_dir+target+date+'_both_clean1.image.tt0'
		else:
			imagenb=my_dir+target+date+'_both_clean1.image'
		print 'Correcting for PB...'
		os.system('rm -rf '+imagenb+'.pbcor')
		os.system('rm -rf '+imagenb+'.pbcor.fits')
		immath(imagename=[imagenb,my_dir+target+date+'_both_clean1.flux'],\
			expr='IM0/IM1',outfile=imagenb+'.pbcor')
		print 'Making fits image...'
		exportfits(imagename=imagenb+'.pbcor',fitsimage=imagenb+'.pbcor.fits')
		imagenb=imagenb+'.pbcor'
		fluxb,errb,unitb,freqb,errb_real=imfit_point(imagenb,my_dir,'I')
		print 'Combined side-band flux density of ',fluxb,' +/- ',errb, unitb
		print 'Local RMS in combined sideband image is: ',errb_real,' Jy'
		dopscb=raw_input('Do you want to do phase selfcal?y or n-->')
		dict_log.append(('phself_both',dopscb))
		if dopscb=='y':
			selfcal_both,scim_both=phselfcal(split_full,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
		date,band_low,'y','',multiscale,robust,weighting)
			fluxb_sc,errb_sc,unitb_sc,freqb_sc,errb_real_sc=imfit_point(scim_both,my_dir)

	#writing imfit result to file
	print 'Writing imfit results to file...'
	resul_file=open(my_dir+'imfit_results.txt','w')
	resul_file.write('#band freq flux err unit err_real\n')
	if dopscl=='n' and dopscu=='n' and dopscb=='n':
		if 'L' in bandsIM:
			resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freql,fluxl,errl,unitl,errl_real))
		if 'U' in bandsIM:
			resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,frequ,fluxu,erru,unitu,erru_real))
		if 'B' in bandsIM:
			resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freqb,fluxb,errb,unitb,errb_real))
	else:
		resul_file.write('Pre-selfcal:\n')
		if 'L' in bandsIM:
			resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freql,fluxl,errl,unitl,errl_real))
		if 'U' in bandsIM:
			resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,frequ,fluxu,erru,unitu,erru_real))
		if 'B' in bandsIM:
			resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freqb,fluxb,errb,unitb,errb_real))
		resul_file.write('Post selfcal:\n')
		if 'L' in bandsIM:
			if dopscl=='y':
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freql_sc,fluxl_sc,errl,unitl_sc,errl_real_sc))
		if 'U' in bandsIM:
			if dopscu=='y':
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,frequ_sc,fluxu_sc,erru_sc,unitu_sc,erru_real_sc))
		if 'B' in bandsIM:
			if dopscb=='y':
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freqb_sc,fluxb_sc,errb_sc,unitb_sc,errb_real_sc))
	resul_file.close()
	###########################################

	###########################################
	#UVfitting
	###########################################
	if uv_fit=='T':
		print 'Performing UV fitting...'
		comp_uv='delta'
		stokes_param='LL'
		print 'Lower side-band...'
		comblsb=split_low.strip('.ms')
		mstransform(vis=split_low, outputvis=comblsb+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
		fitfulluv_low=uvm.uvmultifit(vis=comblsb+'_mstrans.ms', spw='', column = "data", \
			uniform=False, model=[comp_uv],stokes = stokes_param,outfile=my_dir+'lsbmodelfit.dat',\
			var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
		src_uv_init_low=fitfulluv_low.result['Parameters'][2]
		src_uv_err_low=fitfulluv_low.result['Uncertainties'][2]
		print 'Upper side-band...'
		combusb=split_high.strip('.ms')
		mstransform(vis=split_high, outputvis=combusb+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
		fitfulluv_high=uvm.uvmultifit(vis=combusb+'_mstrans.ms', spw='', column = "data", \
			uniform=False, model=[comp_uv],stokes = stokes_param,outfile=my_dir+'usbmodelfit.dat',\
			var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
		src_uv_init_high=fitfulluv_high.result['Parameters'][2]
		src_uv_err_high=fitfulluv_high.result['Uncertainties'][2]
		print 'Combined side-band...'
		combboth=split_full.strip('.ms')
		mstransform(vis=split_full, outputvis=combboth+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
		fitfulluv=uvm.uvmultifit(vis=combboth+'_mstrans.ms', spw='', column = "data", \
			uniform=False, model=[comp_uv],stokes = stokes_param,outfile=my_dir+'bothmodelfit.dat',\
			var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
		src_uv_init=fitfulluv.result['Parameters'][2]
		src_uv_err=fitfulluv.result['Uncertainties'][2]

		print 'Writing uvfit results to file...'
		resuluv_file=open(my_dir+'uvfit_results.txt','w')
		resuluv_file.write('{0} {1} {2} {3}\n'.format(band,'lower',src_uv_init_low,src_uv_err_low))
		resuluv_file.write('{0} {1} {2} {3}\n'.format(band,'upper',src_uv_init_high,src_uv_err_high))
		resuluv_file.write('{0} {1} {2} {3}\n'.format(band,'combined',src_uv_init,src_uv_err))
		resuluv_file.close()

if doImage=='T':
	print 'All side-band data sets have been reduced and imaged.'
else:
	print 'All side-band data sets have been reduced. No imaging was performed.'

###########################################

print 'Cleaning up...'
os.system('rm -rf casa*.log')
os.system('rm -rf *.last')
print 'Writing user_input log file...'
writeDict(dict_log, my_dir+'user_input_'+obsDate+'.logg',str(datetime.datetime.now()))
print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'

