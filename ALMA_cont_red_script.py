#################################################
#ALMA CASA Continuum Reduction Script
#################################################
'''CASA script to be used for the flagging, calibration and imaging of ALMA Continuum Data (Cycle 4 and above)
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) Calibrated Split MS for full band (and each spw) -- data_[target]_[obsDate]_[band]_[TID]_([spw]_)calibrated.ms
        (2) Continuum images for full band (and each epw) -- [target]_[obsDate]_[band]_[TID]_([spw]_)im.image.pbcor
        (3) File of flux densities from image(stokes)/UV plane fitting -- imfit_results_[band].txt/uvfit_results_[band].txt
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
       - All output images are also converted to fits format (just append .fits to end of images above)
       - This script is intended to be used with raw data; can be asdm or ms format.
       - All input logged in user_input.logg.
Written by: Alex J. Tetarenko
Last Updated: Dec 22 2017
Works in CASA-5.1.1 now!'''

print '##################################################'
print 'Welcome to Alexs ALMA Continuum Reduction Script'
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
from utils_CASA import imfit_point, phselfcal, writeDict,source_dict_create,setjy_parse,last_field_parse
from ekoch_casa_tools import set_imagermode,has_field,set_cellsize,set_imagesize,find_expected_beams,getBaselinePercentile,get_mosaic_info
import uvmultifit as uvm
from astropy.io import ascii
import analysisUtils as au

#define output directory
my_dir='/mnt/bigdata/tetarenk/ALMA_maxi1535/my_red/ep1/b6/'
if not os.path.isdir(my_dir):
	os.system('sudo mkdir '+my_dir)
	os.system('sudo chown ubuntu '+my_dir)
	os.system('sudo chmod -R u+r '+my_dir) 
	os.system('sudo chmod -R u+w '+my_dir)
	os.system('sudo chmod -R u+x '+my_dir)
print 'You have set your output directory to ', my_dir
print 'All output images & intermediate data products are put in this directory.\n'

#param file location
param_dir_file='/mnt/bigdata/tetarenk/ALMA_maxi1535/params_alma.txt'
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
raw_type=data_params.raw_type
d_name=data_params.d_name
obsDate=data_params.obsDate
target=data_params.target
band=data_params.band
doImage=data_params.doImage
bandsIM=data_params.bandsIM
#general image params
mythreshold=data_params.mythreshold
myimsize=data_params.myimsize
mycell=data_params.mycell
myniter=data_params.myniter
mynterms=data_params.mynterms
mystokes=data_params.mystokes
multiscale=data_params.multiscale
robust=data_params.robust
weighting=data_params.weighting
outlierfile=data_params.outlierfile
#mask options
mymask=data_params.mymask
#uv fitting
uv_fit=data_params.uv_fit

#write all variables to log dictionary
dict_log.extend([('raw_data_type',raw_type),('d_name_input',d_name),('obsDate', obsDate),('target',target),\
	('bands_input',band),('do imaging',doImage),('bands to image',bandsIM),('weighting',weighting),('outfile',outlierfile),\
	('mythreshold',mythreshold),('myimsize',myimsize),('mycell',mycell),('myniter',myniter),('mynterms',mynterms),\
	('mystokes',mystokes),('mymask',mymask),('uv_fit',uv_fit),('multiscale',multiscale),('robust',robust)])
#################################################

#################################################
#Convert & Split Data Section
#################################################
#convert raw data to MS if needed
if raw_type=='asdm':
	print 'Converting raw data from asdm to MS format...'
	os.system('rm -rf '+my_dir+'data_MS.ms')
	os.system('rm -rf '+my_dir+'data_MS.ms.flagversions')
	os.system('rm -rf '+my_dir+'data_'+target+'_'+obsDate+'_'+band+'pretw.ms')
	os.system('rm -rf '+my_dir+'data_'+target+'_'+obsDate+'_'+band+'pretw.ms.flagversions')
	importasdm(asdm = d_name, vis=my_dir+'data_MS', asis='*')
	print 'Making split MS for ', band,' band...'
	split(vis=my_dir+'data_MS',outputvis=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_'+'pretw.ms',datacolumn='data')
	os.system('rm -rf '+my_dir+'data_MS.ms')
	os.system('rm -rf '+my_dir+'data_MS.ms.flagversions')
else:
	print 'Raw data is already in MS format'
	print 'Splitting out raw data into new MS...'
	os.system('rm -rf '+my_dir+'data_'+target+'_'+obsDate+'_'+band+'pretw.ms')
	os.system('rm -rf '+my_dir+'data_'+target+'_'+obsDate+'_'+band+'pretw.ms.flagversions')
	split(vis=d_name,outputvis=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_'+'pretw.ms',datacolumn='data')

ms_name=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_'+'pretw.ms'
flagmanager(vis=ms_name,mode='save',versionname=target+'_'+obsDate+'_'+band+'_initial',comment='after initial split')

#listobs first
if not os.path.isfile(my_dir+obsDate+'_listfile.txt'):
	print 'Making listobs file...'
	listobs(ms_name,listfile=my_dir+obsDate+'_listfile.txt')
else:
	print 'Listobs file already exists'
print 'Opening listobs file...'
os.system('pluma '+my_dir+obsDate+'_listfile.txt &')
raw_input('Please press enter when ready to continue.')

#enter some data specifics
num_spw=raw_input('Please enter the total number of spws in the data set. e.g. 25--> ')
science_spw=au.getScienceSpws(ms_name,intent='OBSERVE_TARGET#ON_SOURCE',returnString=True).split(',')
bpf_cal=raw_input('Please enter bandpass cals (if >1 seperate by comma)--> ')
bpf_lst=bpf_cal.split(',')
bp_scan=raw_input('Please enter bandpass cal scan(s), e.g., 2,3--> ')
flux_cal=raw_input('Please enter flux cals (if >1 seperate by comma)--> ')
flux_lst=flux_cal.split(',')
second_cal=raw_input('Please enter second cals (if >1 seperate by comma)--> ')
second_lst=second_cal.split(',')
prim_cal,check_src=raw_input('Of these second cals, which is primary and which is check src, e.g., prim,check--> ').split(',')
target_id=raw_input('Please enter target id (if >1 seperate by comma)--> ')
target_lst=target_id.split(',')

cal_lst=list(set(bpf_lst+second_lst+flux_lst))
full_lst=list(set(bpf_lst+second_lst+flux_lst+target_lst))
cal_lst.sort()
full_lst.sort()

dict_log.append(('number of spws',num_spw))
dict_log.append(('bandpass cals',bpf_cal))
dict_log.append(('flux cals',flux_cal))
dict_log.append(('second cals',second_cal))
dict_log.append(('targets',target_id))
dict_log.append(('prim cal',prim_cal))
dict_log.append(('check cal',check_src))

#ant postions
print 'Plotting antenna positions. Please choose a reference antenna near centre of array.'
os.system('rm -rf '+my_dir+'antennas_'+obsDate+'.png')
plotants(vis=ms_name,figfile=my_dir+'antennas_'+obsDate+'.png')
ref_ant=raw_input('Please enter reference antenna. e.g., DV23-->')
dict_log.append(('ref_ant',ref_ant))
#################################################

#################################################
#Tsys and WVR Section
#################################################
#Each ms include tables that contain the Tsys and WVR measurements associated with 
#the visibility data. We must create calibration tables from these measurements.

cal_table_prefix=my_dir+obsDate+'_'+band
wvr_cal=my_dir+obsDate+'_'+band+'.wvr'
tsys_cal=my_dir+obsDate+'_'+band+'.tsys'

#wvr
print 'Computing WVR corrections...'
os.system('rm -rf '+wvr_cal)
wvrgcal(vis=ms_name, caltable=wvr_cal, segsource=True, toffset=-1)

#tsys
print 'Computing Tsys corrections...'
os.system('rm -rf '+tsys_cal)
gencal(vis=ms_name, caltable=tsys_cal, caltype='tsys')

#figure out which sources have tsys measurments and which dont
print 'Showing listobs for Tsys measurments...'
os.system('rm -rf '+my_dir+obsDate+'_listfile_tsys.txt')
listobs(vis=ms_name, verbose=False, selectdata=True, intent='CALIBRATE_ATMOS*',\
	listfile=my_dir+obsDate+'_listfile_tsys.txt')
os.system('pluma '+my_dir+obsDate+'_listfile_tsys.txt &')
raw_input('Please press enter when ready to continue.')

tsys_fields=raw_input('Please enter fields with Tsys measurments, e.g., 1,2,3--> ').split(',')
tsys_fields_wo=raw_input('Please enter fields without Tsys measurments, e.g., 1,2,3--> ').split(',')
tsys_spws=raw_input('Please enter multi-channel spws with Tsys measurments, e.g., 17,18,19--> ').split(',')
dict_log.append(('tsys fields',tsys_fields))
dict_log.append(('no tsys fields',tsys_fields_wo))
dict_log.append(('tsys spws',tsys_spws))

#tsys vs time
print 'Plotting Tsys vs time...cycle through ants'
plotcal(caltable=tsys_cal, xaxis='time',yaxis='tsys',plotrange=[0,0,100,500],\
	iteration='antenna',subplot=421,poln='',spw=tsys_spws[0],showgui=True)
raw_input('Please press enter when ready to continue.')
#tsys vs freq
print 'Plotting Tsys vs freq...cycle through ants'
plotcal(caltable=tsys_cal, xaxis='freq', yaxis='tsys', spw='',subplot=111, iteration='antenna')
raw_input('Please press enter when ready to continue.')
tsys_detail=raw_input('Do you want a detailed look at Tsys?y or n?--> ')
if tsys_detail=='y':
	#tsys vs freq
	print 'Check Tsys solutions. Look out for bad ants/times for flagging.'
	print 'Plotting Tsys vs freq...'
	for field_TS in tsys_fields:
		for spw_TS in tsys_spws:
			print ' Field:', field_TS, ' SPW: ',spw_TS,
			plotbandpass(caltable=tsys_cal, xaxis='freq',yaxis='amp',showtsky=True,\
				subplot=22,field=field_TS,spw=spw_TS, interactive=True,plotrange=[0,0,0,0])
			raw_input('Please press enter when ready to continue.')
	print 'Plotting Tsys vs freq...all scans'
	for field_TS in tsys_fields:
		for spw_TS in tsys_spws:
			print ' Field:', field_TS, ' SPW: ',spw_TS,
			plotbandpass(caltable=tsys_cal, xaxis='freq',yaxis='amp',showtsky=True,\
				subplot=22,field=field_TS,spw=spw_TS, interactive=True,plotrange=[0,0,0,0],overlay='time')
			raw_input('Please press enter when ready to continue.')
	print 'Plotting Tsys vs freq...all ants'
	for field_TS in tsys_fields:
		for spw_TS in tsys_spws:
			print ' Field:', field_TS, ' SPW: ',spw_TS,
			plotbandpass(caltable=tsys_cal, xaxis='freq',yaxis='amp',showtsky=True,\
				subplot=22,field=field_TS,spw=spw_TS, interactive=True,plotrange=[0,0,0,0],overlay='antenna')
			raw_input('Please press enter when ready to continue.')

#apply tsys and wvr corrections
print 'Applying corrections for sources with Tsys...'
#creat spw map if needed for those spws without tsys measurments
tsysmap = range(int(num_spw))
print 'The science spws are: ', ",".join(science_spw)
map_spw=raw_input('Please enter the tsys spws to map these to if needed. e.g., 1,2,5,7 (enter if none)--> ').split(',')
dict_log.append(('map spws',map_spw))
if len(map_spw)==len(science_spw):
	for s in range(0,len(map_spw)):
		if map_spw[s] != '':
			tsysmap[int(science_spw[s])] = int(map_spw[s])
#apply
for field_TS in tsys_fields:
	applycal(vis=ms_name, spw=",".join(science_spw), field=field_TS, gainfield=field_TS,\
		gaintable=[tsys_cal, wvr_cal], spwmap=[tsysmap,[]],interp=['linear,spline','nearest'],\
		flagbackup=False, calwt=True)

print 'Applying corrections for sources without Tsys...'
#enter replacments fields for those sources without tsys measurments
print 'The fields without tsys are: ', ",".join(tsys_fields_wo)
nocal_replace=raw_input('Please enter the replacement fields closest in time/space, e.g., 1,2,3--> ').split(',')
dict_log.append(('tsys replacment fields',nocal_replace))
#apply
for ts in range(0,len(tsys_fields_wo)):
	applycal(vis=ms_name, spw=",".join(science_spw), field=tsys_fields_wo[ts],\
		gaintable=[tsys_cal, wvr_cal], spwmap=[tsysmap,[]],gainfield=[nocal_replace[ts],tsys_fields_wo[ts]],\
		interp=['linear,spline','nearest'],flagbackup=False,calwt=True)

print 'Check Tsys/WVR corrected data...'
print 'Cycle through spws, looks for bad channels.'
strongF=raw_input('Please enter the fields to check (usually the brightest cals), e.g., 1,2,3--> ')
dict_log.append(('brightest cals for tys check',strongF))
print 'Before...'
plotms(vis=ms_name,spw=",".join(science_spw),xaxis='frequency',yaxis='amp',field=strongF,\
	antenna='*&*',avgtime='1e8',avgscan=True,coloraxis='field',iteraxis='spw',\
	xselfscale=True,ydatacolumn='data',showgui=True)
raw_input('Please press enter when ready to continue.')
print 'After...'
plotms(vis=ms_name,spw=",".join(science_spw),xaxis='frequency',yaxis='amp',field=strongF,\
	antenna='*&*',avgtime='1e8',avgscan=True,coloraxis='field',iteraxis='spw',\
	xselfscale=True,ydatacolumn='corrected',showgui=True)
raw_input('Please press enter when ready to continue.')

#split out science data
print 'Splitting out science data...'
os.system('rm -rf ' + my_dir+'data_'+target+'_'+obsDate+'_'+band+'_'+'posttw.ms')
split(vis=ms_name,outputvis=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_'+'posttw.ms',\
	datacolumn='corrected',spw=",".join(science_spw))
ms_name=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_'+'posttw.ms'
#listobs
print 'Opening listobs for split data set...'
os.system('rm -rf '+my_dir+obsDate+'_listfile_aftertw.txt')
listobs(ms_name,listfile=my_dir+obsDate+'_listfile_aftertw.txt')
os.system('pluma '+my_dir+obsDate+'_listfile_aftertw.txt &')
raw_input('Please press enter when ready to continue.')
flagmanager(vis=ms_name,mode='save',versionname=target+'_'+obsDate+'_'+band+'_tsyswvr',comment='after tsys/wvr')
#################################################


#################################################
#Flagging and Inspection Section
#################################################
print 'Starting apriori flagging...'
#autocorrelation data
print 'Flagging autocorrelation data...'
flagdata(vis=ms_name, mode='manual', autocorr=True, action = 'apply', flagbackup=True)
#shadowed data
print 'Flagging shadowed data...'
flagdata(vis=ms_name, mode='shadow',tolerance=12.0, action='apply', flagbackup=False)
#pointing/atmosphere/dummy scan data
print 'Flagging pointing and atmosphere data...'
flagdata(vis=ms_name,mode='manual',intent='*POINTING*', flagbackup=False)
flagdata(vis=ms_name,mode='manual',intent='*ATMOSPHERE*', flagbackup=False)
dum_scan=raw_input('Please enter any dummy scans to flag, e.g., 1,2,3 (enter if none)--> ')
dict_log.append(('flagged dummy scans',dum_scan))
flagdata(vis=ms_name, flagbackup=True, mode='manual', scan=dum_scan)
#bad ants from log or wvr/tsys
badants=raw_input('Please enter bad ants from log or wvr/tsys correction, e.g., DA43,DA52 (enter if none)--> ')
dict_log.append(('flagged bad ants',badants))
if badants=='':
	print 'No bad ants to flag.'
else:
	print 'Flagging bad ants: ',badants
	flagdata(vis=ms_name, mode='manual', antenna=badants, action='apply', flagbackup=False)
#end channels
flag_end=raw_input('Do you want to flag end channels? y or n--> ')
dict_log.append(('flag end channels',flag_end))
if flag_end=='y':
	beg=raw_input('Please enter beginning channels to flag. e.g., 0~3-->')
	endd=raw_input('Please enter end channels to flag. e.g., 124~127-->')
	dict_log.append(('beg channels flagged',beg))
	dict_log.append(('end channels flagged',endd))
	print 'Flagging beginning channels ',beg,' and end channels ',endd,' ...'
	flagdata(vis=ms_name,spw='0~'+str(len(science_spw)-1)+':'+beg,field='',antenna='')
	flagdata(vis=ms_name,spw='0~'+str(len(science_spw)-1)+':'+endd,field='',antenna='')
else:
	print 'End channels not flagged.'

print 'Starting detailed flagging...'
print 'Look for obvious bad data.'
print 'Keep track of ants/spws/channels to flag. You will be prompted for their values after plotting.'
print '(1) Amp vs time...iterating over spw'
print 'Look for bad or low amp data.'
plotms(vis=ms_name,spw='',xaxis='time',yaxis='amp',field='',avgchannel='128',\
	coloraxis='field',iteraxis='spw',showgui=True)
raw_input('Please press enter when ready to continue.')
print '(2) Phase vs freq (ref_ant only,bpf cal only)...iterating over spw'
print 'Significant delay errors or baseline errors will show up as phase slopes here.'
plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='phase',field=bpf_cal,antenna=ref_ant,\
	avgtime='1e8',avgscan=True,coloraxis='baseline',iteraxis='spw',xselfscale=True,showgui=True)
raw_input('Please press enter when ready to continue.')
print '(2) Amp vs freq (second cals and target only)...iterating over spw'
print 'Look for birdies.'
plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=second_cal+','+target_id,\
	avgtime='1e8',avgscan=True,coloraxis='field',iteraxis='spw',xselfscale=True,yselfscale=True,showgui=True)
raw_input('Please press enter when ready to continue.')
print '(3) BP and flux cal spectral features...'
print 'Look for absorbtion and emission features.'
plotms(vis=ms_name,spw='',xaxis='channel',yaxis='amp',field=bpf_cal,\
	avgtime='1e8',coloraxis='field',iteraxis='spw',showgui=True)
raw_input('Please press enter when ready to continue.')
if flux_cal != bpf_cal:
	plotms(vis=ms_name,spw='',xaxis='channel',yaxis='amp',field=flux_cal,\
		avgtime='1e8',coloraxis='field',iteraxis='spw',showgui=True)
raw_input('Please press enter when ready to continue.')
print 'Flagging...'
badasf=raw_input('Please enter bad ant,spw,field,corr,and timerange to flag (enter if none). e.g., DV10,DA12;5:4~9;3;YY;9:52:10.0~9:53:10.0 ;5;3;;-->').split(' ')
dict_log.append(('flags',badasf))
if '' in badasf:
	print 'Nothing to flag.'
else:
	print 'Flagging selected data.'
	for i in range(0,len(badasf)):
		strg=badasf[i].split(';')
		flagdata(vis=ms_name,flagbackup=True, mode='manual', antenna=strg[0],spw=strg[1],field=strg[2],correlation=strg[3],timerange=strg[4])
print 'Final check of flagged data...'
plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=second_cal+','+target_id,\
	avgtime='1e8',avgscan=True,coloraxis='field',iteraxis='spw',xselfscale=True,yselfscale=True,showgui=True)
raw_input('Please press enter when ready to continue.')
flag_again=raw_input('Do you need to do more flagging? y or n-->')
dict_log.append(('flag again',flag_again))
while flag_again=='y':
	badasf2=raw_input('Please enter bad ant,spw,field,corr,and timerange to flag (enter if none). e.g., DV10,DA12;5:4~9;3;YY;9:52:10.0~9:53:10.0 ;5;3;;-->').split(' ')
	raw_input('press enter to flag selected data.')
	dict_log.append(('flags2',badasf2))
	if '' in badasf2:
		print 'Nothing to flag.'
	else:
		print 'Flagging selected data.'
		for i in range(0,len(badasf2)):
			strg2=badasf2[i].split(';')
			flagdata(vis=ms_name,flagbackup=True, mode='manual', antenna=strg2[0],spw=strg2[1],field=strg2[2],correlation=strg2[3],timerange=strg2[4])
	print 'Plotting...'
	plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=second_cal+','+target_id,\
		avgtime='1e8',avgscan=True,coloraxis='field',iteraxis='spw',xselfscale=True,yselfscale=True,showgui=True)
	raw_input('Please press enter when ready to continue.')
	raw_input('press enter to continue.')
	flag_again=raw_input('Do you need to do more flagging? y or n-->')
#save flags
flagmanager(vis=ms_name,mode='save',versionname=target+'_'+obsDate+'_'+band+'_flagging',comment='after flagging')
writeDict(dict_log, my_dir+'user_input_'+obsDate+'_flag.logg',str(datetime.datetime.now()))
#################################################

#################################################
#Antenna positions correction
#################################################
#make list of pre gain tables to be applied in tandem below
gt_lst=[]

do_ant=raw_input('Do you wish to do antenna positions corrections? y or n?--> ')
dict_log.append(('ant corrections',do_ant))
if do_ant=='y':
	#antpos
	cor_ants=raw_input('Please enter ants to be corrected, e.g., DA48,DA65--> ')
	raw_input('Please press enter when ready to continue.')
	cor_pos=raw_input('Please enter set of three number corrections for each antenna in tandem, e.g., -1,-2,3,4,5,6--> ').split(',')
	cor_lst=[float(i) for i in cor_pos]
	gencal(caltable=cal_table_prefix+'.antpos',
		   vis=ms_name, caltype='antpos',
		   parameter=cor_lst, antenna=cor_ants)
	gt_lst.append(cal_table_prefix+'.antpos')
#################################################

#################################################
#Set up flux model
#################################################
print 'Setting the flux model...'
flux_man=raw_input('Do you need to manually set the flux cal values? y or n-->')
dict_log.append(('manually set flux model',flux_man))
if flux_man=='n':
	setjy(vis=ms_name,field=flux_cal,usescratch=False,standard='Butler-JPL-Horizons 2012',scalebychan=False)
	print 'Plotting flux model...'
	plotms(vis=ms_name,field=flux_cal,xaxis='uvdist',yaxis='amp',coloraxis='spw',ydatacolumn='model',showgui=True,avgtime='1e8')
	raw_input('Please press enter when ready to continue.')
else:
	f_man=raw_input('Please enter cal name--> ')
	#nu_man=raw_input('Please enter central freq in GHz for each science spw, e.g., 90.496,92.434--> ').split(',')
	nu_man=au.getScienceFrequencies(ms_name)
	#ds=raw_input('Please enter date of obs, e.g. 20170911--> ')
	ds=au.getObservationStartDate(ms_name).split(' ')[0]
	dict_log.append(('manual field',f_man))
	dict_log.append(('manual nus',nu_man))
	#dict_log.append(('manual date',ds))
	for i in range(0,len(science_spw)):
		sj_val=au.getALMAFlux(f_man,str(nu_man[i]/1e9)+' GHz',date=ds.replace('-',''))
		dict_log.append(('manual si',sj_val['spectralIndex']))
		setjy(fluxdensity=[float(sj_val['fluxDensity']), 0.0, 0.0, 0.0], scalebychan=True,\
			vis=ms_name, spw=str(i),spix=float(sj_val['spectralIndex']),\
			field=flux_cal, reffreq=str(nu_man[i]/1e9)+' GHz',\
			selectdata=True, standard='manual', usescratch=True)
#################################################

#################################################
#Bandpass cal
#################################################

#first look at bandpass cal in time
print 'Examining bandpass cal...'
print 'First look for middle channels for initial phase cal before bandpass...'
plotms(vis=ms_name, spw='0~'+str(len(science_spw)-1), antenna=ref_ant, xaxis='freq', yaxis='amp',iteraxis='spw',\
	coloraxis='spw', symbolshape = 'circle',avgtime='1e8',avgscan=True,field=bpf_cal,scan=bp_scan)
raw_input('Please press enter when ready to continue.')
bp_chan=raw_input('Please enter middle channels, e.g., 0~3:60~80--> ')
dict_log.append(('BP channels',bp_chan))
print 'Initial phase solution before bandpass...'
gaincal(vis=ms_name,caltable=cal_table_prefix+'.bpphase.gcal',field=bpf_cal,\
	spw=bp_chan,refant=ref_ant,\
	calmode='p',solint='int',minsnr=2.0,minblperant=4,intent='CALIBRATE_BANDPASS#ON_SOURCE')
gt_lst.append(cal_table_prefix+'.bpphase.gcal')
print 'Plotting solutions...'
print 'XX'
plotcal(caltable=cal_table_prefix+'.bpphase.gcal',xaxis='time',yaxis='phase',spw='',antenna='',\
	iteration='antenna',subplot=421,plotrange=[0,0,-180,180],poln='X',showgui=True)
raw_input('Please press enter when ready to continue.')
print 'YY'
plotcal(caltable=cal_table_prefix+'.bpphase.gcal',xaxis='time',yaxis='phase',spw='',antenna='',\
	iteration='antenna',subplot=421,plotrange=[0,0,-180,180],poln='X',showgui=True)
raw_input('Please press enter when ready to continue.')

print 'BP solution...'
bandpass(vis=ms_name,caltable=cal_table_prefix+'.bandpass.bcal',field=bpf_cal,spw='',combine='',\
	refant=ref_ant,solint='inf',solnorm=True,minblperant=4, bandtype='B',\
	fillgaps=3,gaintable=gt_lst,intent='CALIBRATE_BANDPASS#ON_SOURCE')
print 'Plotting solutions...'
for ii in range(0,len(science_spw)):
	print 'Amp for Spw ', ii
	plotbandpass(cal_table_prefix+'.bandpass.bcal',xaxis='freq',yaxis='amp', spw=str(i),antenna='',\
		subplot=42, showatm=True,interactive=True)
	raw_input('Please press enter when ready to continue.')
	print 'Phase for Spw ', ii
	plotbandpass(cal_table_prefix+'.bandpass.bcal',xaxis='freq',yaxis='phase', spw=str(i),antenna='',\
		subplot=42, showatm=True,interactive=True)
	raw_input('Please press enter when ready to continue.')
gt_lst.remove(cal_table_prefix+'.bpphase.gcal')
gt_lst.append(cal_table_prefix+'.bandpass.bcal')
#################################################

#################################################
#antenna-based phase and amplitude gain cal
#################################################
#(1) we will do amp and phase separetly as the phase changes on a much shorter timescale than the amplitude
#(2) also do seperate olution for application to target later. We do this as phase-scatter within a scan can
#dominate the interpolation between calibrator scans, so solve for the phase on the scan time rather then int.
print 'Gain calibration...'

#phase for cals
print 'Short phase solution for cals...'
gaincal(vis=ms_name,caltable=cal_table_prefix+'.intphase.gcal',field=",".join(cal_lst),\
	spw='',refant=ref_ant,calmode='p',solint='int',minsnr=2.0,minblperant=4,gaintable=gt_lst)
#phase for target
print 'Longer phase solution for target...'
gaincal(vis=ms_name,caltable=cal_table_prefix+'.scanphase.gcal',field=",".join(cal_lst),\
	spw='',refant=ref_ant,calmode='p',solint='inf',minsnr=2.0,minblperant=4,gaintable=gt_lst)
gt_lst.append(cal_table_prefix+'.intphase.gcal')
#amp
print 'Amp solution...'
gaincal(vis=ms_name,caltable=cal_table_prefix+'.amp.gcal',field=",".join(cal_lst),\
	spw='',refant=ref_ant,calmode='ap',solint='inf',minsnr=2.0,minblperant=4,gaintable=gt_lst)

print 'Plotting phase solutions...'
print 'First short solution...'
plotcal(caltable=cal_table_prefix+'.intphase.gcal',xaxis='time',yaxis='phase',\
	antenna='',spw='',field=",".join(cal_lst),iteration='antenna',subplot=421,\
	plotrange=[0,0,-180,180],poln='X',showgui=True)
raw_input('Please press enter when ready to continue.')
plotcal(caltable=cal_table_prefix+'.intphase.gcal',xaxis='time',yaxis='phase',\
	antenna='',spw='',field=",".join(cal_lst),iteration='antenna',subplot=421,\
	plotrange=[0,0,-180,180],poln='Y',showgui=True)
raw_input('Please press enter when ready to continue.')
print 'Second long solution...'
plotcal(caltable=cal_table_prefix+'.scanphase.gcal',xaxis='time',yaxis='phase',\
	antenna='',spw='',field=",".join(cal_lst),iteration='antenna',subplot=421,\
	plotrange=[0,0,-180,180],poln='X',showgui=True)
raw_input('Please press enter when ready to continue.')
plotcal(caltable=cal_table_prefix+'.scanphase.gcal',xaxis='time',yaxis='phase',\
	antenna='',spw='',field=",".join(cal_lst),iteration='antenna',subplot=421,\
	plotrange=[0,0,-180,180],poln='Y',showgui=True)
raw_input('Please press enter when ready to continue.')
print 'Look for residual phase error...'
print 'Should not be scatter of more than a few degrees'
plotcal(caltable=cal_table_prefix+'.amp.gcal',xaxis='time',\
	yaxis='phase',antenna='',spw='',field=",".join(cal_lst),\
	plotrange=[0,0,-1,1],iteration='antenna',subplot=421,showgui=True)
raw_input('Please press enter when ready to continue.')
print 'Plotting amp solutions...'
plotcal(caltable=cal_table_prefix+'.amp.gcal',xaxis='time',yaxis='amp',\
	antenna='',iteration='antenna',subplot=421,spw='',poln='X',\
	plotrange=[0,0,0.0,0],showgui=True)
raw_input('Please press enter when ready to continue.')
plotcal(caltable=cal_table_prefix+'.amp.gcal',xaxis='time',yaxis='amp',\
	antenna='',iteration='antenna',subplot=421,spw='',poln='Y',\
	plotrange=[0,0,0.0,0],showgui=True)
raw_input('Please press enter when ready to continue.')
gt_lst.append(cal_table_prefix+'.amp.gcal')
#################################################

#################################################
#flux cal
#################################################
print 'Flux calibration...'

maps=raw_input('Was the reference field observed in all spws that transfer fields were? y or n?--> ')
dict_log.append(('fluxscale maps',maps))
if maps=='y':
	fluxscale(vis=ms_name,caltable=cal_table_prefix+'.amp.gcal',\
	fluxtable=cal_table_prefix+'.flux.cal',reference=flux_cal,\
	listfile=cal_table_prefix+'.fluxscale.txt')
else:
	maps_val=raw_input('Please enter map, e.g., 0,0,3,3--> ').split(',')
	fluxscale(vis=ms_name,caltable=cal_table_prefix+'.amp.gcal',\
		fluxtable=cal_table_prefix+'.flux.cal',reference=flux_cal,\
		refspwmap=[int(maps_val[0]),int(maps_val[1]),int(maps_val[2]),int(maps_val[3])],\
		listfile=cal_table_prefix+'.fluxscale.txt')
print 'Plotting solutions...'
plotcal(caltable=cal_table_prefix+'.flux.cal',xaxis='time',yaxis='amp',\
	antenna='',iteration='antenna',subplot=421,spw='',poln='X',\
	plotrange=[0,0,0.0,0],showgui=True)
raw_input('Please press enter when ready to continue.')
plotcal(caltable=cal_table_prefix+'.flux.cal',xaxis='time',yaxis='amp',\
	antenna='',iteration='antenna',subplot=421,spw='',poln='X',\
	plotrange=[0,0,0.0,0],showgui=True)
raw_input('Please press enter when ready to continue.')
gt_lst.remove(cal_table_prefix+'.amp.gcal')
gt_lst.append(cal_table_prefix+'.flux.cal')
#################################################

#################################################
#apply cal solutions
#################################################
print 'Apply calibration solutions...'
if do_ant=='n':
	#bandpass cal
	print 'Bandpass cals.'
	for i in range(0,len(bpf_lst)):
		applycal(vis=ms_name,field=bpf_lst[i],gaintable=gt_lst,\
			interp=['nearest','nearest','nearest'],\
			gainfield=[bpf_lst[i],bpf_lst[i],bpf_lst[i]], flagbackup=True, calwt=False)
	#flux cal if different from bandpass cal
	print 'Flux cals.'
	for i in range(0,len(flux_lst)):
		if flux_lst[i] not in bpf_lst:
			applycal(vis=ms_name,field=flux_lst[i],gaintable=gt_lst,\
				interp=['nearest','nearest','nearest'],\
				gainfield=[bpf_cal,flux_lst[i],flux_lst[i]], flagbackup=True, calwt=False)
	#primary phase cal if different from bandpass and flux cals
	print 'Primary phase cal.'
	for i in range(0,len(second_lst)):
		if second_lst[i] not in bpf_lst and second_lst[i] not in flux_lst and second_lst[i] != check_src:
			applycal(vis=ms_name,field=second_lst[i],gaintable=gt_lst,\
				interp=['nearest','nearest','nearest'],\
				gainfield=[bpf_cal,second_lst[i],second_lst[i]], flagbackup=True, calwt=False)
	#check source
	print 'Secondary phase cal/check source.'
	if check_src != '':
		applycal(vis=ms_name,field=check_src,gaintable=[cal_table_prefix+'.bandpass.bcal',\
			cal_table_prefix+'.scanphase.gcal',cal_table_prefix+'.flux.cal'],\
			interp=['nearest','linear','linear'],gainfield=[bpf_cal,prim_cal,prim_cal],\
			flagbackup=True, calwt=False)
	#target
	print 'Target.'
	for i in range(0,len(target_lst)):
		applycal(vis=ms_name,field=target_lst[i],gaintable=[cal_table_prefix+'.bandpass.bcal',\
			cal_table_prefix+'.scanphase.gcal',cal_table_prefix+'.flux.cal'],\
			interp=['nearest','linear','linear'],\
			gainfield=[bpf_cal,second_cal,second_cal], flagbackup=True, calwt=False)
elif do_ant=='y':
	#bandpass cal
	print 'Bandpass cals.'
	for i in range(0,len(bpf_lst)):
		applycal(vis=ms_name,field=bpf_lst[i],gaintable=gt_lst,\
			interp=['','nearest','nearest','nearest'],\
			gainfield=['',bpf_lst[i],bpf_lst[i],bpf_lst[i]], flagbackup=True, calwt=False)
	#flux cal if different from bandpass cal
	print 'Flux cals.'
	for i in range(0,len(flux_lst)):
		if flux_lst[i] not in bpf_lst:
			applycal(vis=ms_name,field=flux_lst[i],gaintable=gt_lst,\
				interp=['','nearest','nearest','nearest'],\
				gainfield=['',bpf_cal,flux_lst[i],flux_lst[i]], flagbackup=True, calwt=False)
	#primary phase cal if different from bandpass and flux cals
	print 'Primary phase cal.'
	for i in range(0,len(second_lst)):
		if second_lst[i] not in bpf_lst and second_lst[i] not in flux_lst and second_lst[i] != check_src:
			applycal(vis=ms_name,field=second_lst[i],gaintable=gt_lst,\
				interp=['','nearest','nearest','nearest'],\
				gainfield=['',bpf_cal,second_lst[i],second_lst[i]], flagbackup=True, calwt=False)
	#check source
	print 'Secondary phase cal/check source.'
	if check_src != '':
		applycal(vis=ms_name,field=check_src,gaintable=[cal_table_prefix+'.antpos',cal_table_prefix+'.bandpass.bcal',\
			cal_table_prefix+'.scanphase.gcal',cal_table_prefix+'.flux.cal'],\
			interp=['','nearest','linear','linear'],gainfield=['',bpf_cal,prim_cal,prim_cal],\
			flagbackup=True, calwt=False)
	#target
	print 'Target.'
	for i in range(0,len(target_lst)):
		applycal(vis=ms_name,field=target_lst[i],gaintable=[cal_table_prefix+'.antpos',cal_table_prefix+'.bandpass.bcal',\
			cal_table_prefix+'.scanphase.gcal',cal_table_prefix+'.flux.cal'],\
			interp=['','nearest','linear','linear'],\
			gainfield=['',bpf_cal,second_cal,second_cal], flagbackup=True, calwt=False)
flagmanager(vis=ms_name,mode='save',versionname=target+'_'+obsDate+'_'+band+'_applycal',comment='after applycal')
#################################################

#################################################
#check corrected data
#################################################
print 'Checking corrected data...'
print '(1) Amp vs time for all sources'
plotms(vis=ms_name,spw='',xaxis='time',yaxis='amp',field='',avgchannel='128',\
	coloraxis='field',iteraxis='spw',ydatacolumn='corrected',showgui=True)
raw_input('Please press enter when ready to continue.')
print '(2) Phase vs time for cals'
plotms(vis=ms_name,spw='',xaxis='time',yaxis='phase',field=",".join(cal_lst),avgchannel='128',
         coloraxis='field',iteraxis='spw',ydatacolumn='corrected',showgui=True)
raw_input('Please press enter when ready to continue.')

if check_src != '':
	print 'Checking phase transfer...'
	print 'Applying check src solutions to itself and compare to transfer solutions.'
	if do_ant=='n':
		applycal(vis=ms_name,field=check_src,gaintable=gt_lst,\
			interp=['nearest','nearest','nearest'],gainfield=[bpf_cal,check_src,check_src],\
			flagbackup=True, calwt=False)
	elif do_ant=='y':
		applycal(vis=ms_name,field=check_src,gaintable=gt_lst,\
			interp=['','nearest','nearest','nearest'],gainfield=['',bpf_cal,check_src,check_src],\
			flagbackup=True, calwt=False)
	print 'Plotting again...'
	print '(1) Amp vs time for all sources'
	plotms(vis=ms_name,spw='',xaxis='time',yaxis='amp',field='',avgchannel='128',\
		coloraxis='field',iteraxis='spw',ydatacolumn='corrected',showgui=True)
	raw_input('Please press enter when ready to continue.')
	print '(2) Phase vs time for cals'
	plotms(vis=ms_name,spw='',xaxis='time',yaxis='phase',field=",".join(cal_lst),\
		avgchannel='128',coloraxis='field',iteraxis='spw',ydatacolumn='corrected',showgui=True)
	raw_input('Please press enter when ready to continue.')
	print 'Reapply first solutions from primary cal to check src...'
	if do_ant=='n':
		applycal(vis=ms_name,field=check_src,gaintable=[cal_table_prefix+'.bandpass.bcal',\
			cal_table_prefix+'.scanphase.gcal',cal_table_prefix+'.flux.cal'],\
			interp=['nearest','linear','linear'],gainfield=[bpf_cal,prim_cal,prim_cal],\
			flagbackup=True, calwt=False)
	elif do_ant=='y':
		applycal(vis=ms_name,field=check_src,gaintable=[cal_table_prefix+'.antpos',cal_table_prefix+'.bandpass.bcal',\
			cal_table_prefix+'.scanphase.gcal',cal_table_prefix+'.flux.cal'],\
			interp=['','nearest','linear','linear'],gainfield=['',bpf_cal,prim_cal,prim_cal],\
			flagbackup=True, calwt=False)

print 'Checking spectral calibration...'
print 'Interate over spws.'
print 'First bandpass and flux cals...'
cals_add=[]
for cal in flux_lst:
	if cal not in bpf_lst:
		cals_add.append(cal)
if cals_add==[]:
	plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=bpf_cal,\
		avgtime='1e8',avgscan=True,coloraxis='field',ydatacolumn='corrected',\
		iteraxis='spw',xselfscale=True,showgui=True)
	raw_input('Please press enter when ready to continue.')
else:
	plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=bpf_cal+','+",".join(cals_add),\
		avgtime='1e8',avgscan=True,coloraxis='field',ydatacolumn='corrected',\
		iteraxis='spw',xselfscale=True,showgui=True)
	raw_input('Please press enter when ready to continue.')
print 'Next phase cals...'
plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=second_cal,\
	avgtime='1e8',avgscan=True,coloraxis='field',ydatacolumn='corrected',\
	iteraxis='spw',xselfscale=True,showgui=True)
raw_input('Please press enter when ready to continue.')
print 'Now the target...'
plotms(vis=ms_name,spw='',xaxis='frequency',yaxis='amp',field=target_id,\
	avgtime='1e8',avgscan=True,coloraxis='field',ydatacolumn='corrected',\
	iteraxis='spw',xselfscale=True,showgui=True)
raw_input('Please press enter when ready to continue.')
extraf=raw_input('Do you need to do additional flagging? y or n-->')
if extraf=='y':
	while extraf=='y':
		badasfextra=raw_input('Please enter bad ant,spw,field, correlation, and timerange to flag (enter if none). e.g., DV10,DA12;5:4~9;3;YY;10:52:11.0~10:53:11.0 ;5;3;-->').split(' ')
		if '' in badasfextra:
			print 'Nothing to flag.'
			extraf=raw_input('Do you need to do additional flagging? y or n-->')
		else:
			print 'Flagging selected data.'
			for i in range(0,len(badasfextra)):
				strge=badasfextra[i].split(';')
				flagdata(vis=ms_name,flagbackup=True, mode='manual', antenna=strge[0],spw=strge[1],field=strge[2],correlation=strge[3],timerange=strge[4])
			raw_input('press enter to continue.')
			extraf=raw_input('Do you need to do additional flagging? y or n-->')
else:
	print 'No extra flagging requested.'
flagmanager(vis=ms_name,mode='save',versionname=target+'_'+obsDate+'_'+band+'_addflag',comment='after additional flagging')
#################################################

#################################################
#split out target data
#################################################
print 'Splitting out target data...'
#full MS
for iii in range(0,len(target_lst)):
	print 'Full band for TID: ',target_lst[iii]
	ms_name_final=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_calibrated.ms'
	os.system('rm -rf '+ms_name_final)
	split(vis=ms_name,outputvis=ms_name_final,datacolumn='corrected',\
		field=target_lst[iii],antenna='',spw='')
	#individual spws
	for k in range(0,len(science_spw)):
		print 'SPW: ', k, ' for TID: ',target_lst[iii]
		ms_name_finalspw=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(k)+'_calibrated.ms'
		os.system('rm -rf '+ms_name_finalspw)
		split(vis=ms_name,outputvis=ms_name_finalspw,datacolumn='corrected',\
			field=target_lst[iii],antenna='',spw=str(k))
#################################################

#################################################
#Imaging
#################################################
if doImage=='T':
	print 'Starting Imaging...'
	dopscbb='n'
	dopscs='n'
	for iii in range(0,len(target_lst)):
		print 'Beginnning Imaging of target id: ', target_lst[iii]
		if 'A' in bandsIM:
			vis=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_calibrated.ms'
			print 'Imaging all base-bands together...'
			if mymask=='':
				print 'Using interactive mode so you can make a mask...'
				print 'Cleaning...'
				os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim*')
				clean(vis=vis,imagename=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim',\
					mode='mfs',imagermode='csclean',imsize=myimsize,cell=mycell,spw='',\
					gain=0.1,weighting=weighting,robust=robust,nterms=mynterms, mask='',usescratch=False,\
					interactive=True,threshold=mythreshold,niter=0,pbcor=False,minpb=0.2,multiscale=multiscale)
				raw_input('Please press enter to continue when you are done.')
			else:
				print 'Cleaning...'
				print 'This make take awhile, please go do something else for a while.'
				os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim*')
				clean(vis=vis,imagename=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim',\
					mode='mfs',imagermode='csclean',imsize=myimsize,cell=mycell,spw='',\
					gain=0.1,weighting=weighting,robust=robust,nterms=mynterms, mask=mymask,usescratch=False,\
					interactive=False,threshold=mythreshold,niter=0,pbcor=False,minpb=0.2,multiscale=multiscale)
				raw_input('Please press enter to continue when you are done.')
			os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.pbcor')
			os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.pbcor.fits')
			if int(mynterms) > 1:
				os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image')
				os.system('sudo mv '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.tt0 '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image')
			immath(imagename=[my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+\
				'fullbandim.image',my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.flux'],\
				expr='IM0/IM1',outfile=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.pbcor')
			exportfits(imagename=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.pbcor',\
				fitsimage=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.pbcor.fits')
			imagenbb=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'fullbandim.image.pbcor'
			fluxbb,errbb,unitbb,freqbb,errbb_real=imfit_point(imagenbb,my_dir,'I')
			print 'Combined base-band flux density of ',fluxbb,' +/- ',errbb, unitbb
			print 'Local RMS in Combined base-band image is: ',errbb_real,' Jy'
			dopscbb=raw_input('Do you want to do phase selfcal?y or n-->')
			dict_log.append(('selfcal full bb',dopscbb))
			if dopscbb=='y':
				selfcal_both,scim_both=phselfcal(vis,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
					obsDate,band,'n',outlierfile,multiscale,robust,weighting)
				fluxbb_sc,errbb_sc,unitbb_sc,freqbb_sc,errbb_real_sc=imfit_point(scim_both,my_dir)
				print 'Combined base-band flux density post selfcal of ',fluxbb_sc,' +/- ',errbb_sc, unitbb_sc
				print 'Local RMS in Combined base-band image post selfcal is: ',errbb_real_sc,' Jy'
		if 'S' in bandsIM:
			spw_presc=[]
			spw_postsc=[]
			print 'Imaging base-band spws separetly...'
			for j in range(0,len(science_spw)):
				print 'Imaging spw: ',j
				vis=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'_calibrated.ms'
				if mymask=='':
					print 'Using interactive mode so you can make a mask...'
					print 'Cleaning...'
					os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im*')
					clean(vis=vis,imagename=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im',\
						mode='mfs',imagermode='csclean',imsize=myimsize,cell=mycell,spw='',\
						gain=0.1,weighting=weighting,robust=robust,nterms=mynterms, mask='',usescratch=False,\
						interactive=True,threshold=mythreshold,niter=0,pbcor=False,minpb=0.2,multiscale=multiscale)
					raw_input('Please press enter to continue when you are done.')
				else:
					print 'Cleaning...'
					print 'This make take awhile, please go do something else for a while.'
					os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im*')
					clean(vis=vis,imagename=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im',\
						mode='mfs',imagermode='csclean',imsize=myimsize,cell=mycell,spw='',\
						gain=0.1,weighting=weighting,robust=robust,nterms=mynterms, mask=mymask,usescratch=False,\
						interactive=False,threshold=mythreshold,niter=0,pbcor=False,minpb=0.2,multiscale=multiscale)
					raw_input('Please press enter to continue when you are done.')
				os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.pbcor')
				os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.pbcor.fits')
				if int(mynterms) > 1:
					os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image')
					os.system('sudo mv '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.tt0 '+my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image')
				immath(imagename=[my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+\
					'_spw'+str(j)+'im.image',my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.flux'],\
					expr='IM0/IM1',outfile=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.pbcor')
				exportfits(imagename=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.pbcor',\
					fitsimage=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.pbcor.fits')
				imagens=my_dir+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_spw'+str(j)+'im.image.pbcor'
				fluxs,errs,units,freqs,errs_real=imfit_point(imagens,my_dir,'I')
				print 'Spw ', j,' flux density of ',fluxs,' +/- ',errs, units
				print 'Local RMS in image is: ',errs_real,' Jy'
				spw_presc.append([freqs,fluxs,errs,units,errs_real])
				dopscs=raw_input('Do you want to do phase selfcal?y or n-->')
				dict_log.append(('selfcal spws',dopscs))
				if dopscs=='y':
					selfcal_s,scim_s=phselfcal(vis,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
						obsDate,band,'n',outlierfile,multiscale,robust,weighting)
					fluxs_sc,errs_sc,units_sc,freqs_sc,errs_real_sc=imfit_point(scim_s,my_dir)
					print 'Combined base-band flux density post selfcal of ',fluxs_sc,' +/- ',errs_sc, units_sc
					print 'Local RMS in Combined base-band image post selfcal is: ',errs_real_sc,' Jy'
					spw_postsc.append([freqs_sc,fluxs_sc,errs_sc,units_sc,errs_real_sc])
		#writing imfit result to file
		print 'Writing imfit results to file...'
		if iii==0:
			resul_file=open(my_dir+'imfit_results'+'_'+band+'.txt','w')
			resul_file.write('#freq flux err unit err_real\n')
		else:
			resul_file=open(my_dir+'imfit_results.txt','a')
		if 'A' in bandsIM:
			if dopscbb=='n':
				resul_file.write('Target'+str(target_id[iii])+':\n')
				resul_file.write('{0} {1} {2} {3} {4}\n'.format(freqbb,fluxbb,errbb,unitbb,errbb_real))
			else:
				resul_file.write('Target'+str(target_id[iii])+':\n')
				resul_file.write('Pre-selfcal:\n')
				resul_file.write('{0} {1} {2} {3} {4}\n'.format(freqbb,fluxbb,errbb,unitbb,errbb_real))
				resul_file.write('Post selfcal:\n')
				resul_file.write('{0} {1} {2} {3} {4}\n'.format(freqbb_sc,fluxbb_sc,errbb_sc,unitbb_sc,errbb_real_sc))
		if 'S' in bandsIM:
			for kk in range(0,len(science_spw)):
				if dopscs=='n':
					vals=spw_presc[kk]
					resul_file.write('Target'+str(target_id[iii])+' SPW'+str(kk)+':\n')
					resul_file.write('{0} {1} {2} {3} {4}\n'.format(vals[0],vals[1],vals[2],vals[3],vals[4]))
				else:
					resul_file.write('Target'+str(target_id[iii])+' SPW'+str(kk)+':\n')
					resul_file.write('Pre-selfcal:\n')
					vals=spw_presc[kk]
					resul_file.write('{0} {1} {2} {3} {4}\n'.format(vals[0],vals[1],vals[2],vals[3],vals[4]))
					resul_file.write('Post selfcal:\n')
					vals2=spw_postsc[kk]
					resul_file.write('{0} {1} {2} {3} {4}\n'.format(vals2[0],vals2[1],vals2[2],vals2[3],vals2[4]))
		resul_file.close()
###########################################

###########################################
#UVfitting
###########################################
if uv_fit=='T':
	print 'UV fitting...'
	for iii in range(0,len(target_lst)):
		print 'Beginnning UV fitting of target id: ', target_lst[iii]
		comp_uv='delta'
		stokes_param=mystokes
		vis=my_dir+'data_'+target+'_'+obsDate+'_'+band+'_TID'+target_lst[iii]+'_calibrated.ms'
		print 'Combined base-band...'
		combboth=vis.strip('.ms')
		mstransform(vis=vis, outputvis=combboth+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
		fitfulluv=uvm.uvmultifit(vis=combboth+'_mstrans.ms', spw='', column = "data",\
			uniform=False, model=[comp_uv],stokes = stokes_param,\
			outfile=my_dir+'combmodelfit.dat',var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,\
			cov_return=False, finetune=False, method="levenberg")
		src_uv_init=fitfulluv.result['Parameters'][2]
		src_uv_err=fitfulluv.result['Uncertainties'][2]
		print 'Writing uvfit results to file...'
		if iii==0:
			resuluv_file=open(my_dir+'uvfit_results'+'_'+band+'.txt','w')
		else:
			resuluv_file=open(my_dir+'uvfit_results'+'_'+band+'.txt','a')
		resuluv_file.write('{0} {1} {2} {3}\n'.format('combined',src_uv_init,src_uv_err))
		resuluv_file.close()
###########################################

###########################################
#Obs info
###########################################
print 'Observation INFO:'

#obs dates in MJD
start_d=au.getObservationStartDate(ms_name_final)
endd_d=au.getObservationStopDate(ms_name_final)
start_m=au.dateStringToMJD(start_d)
end_m=au.dateStringToMJD(endd_d)
mjd_err=((end_m-start_m)/2.)
mjd_mid=mjd_err+start_m
dict_log.append(('mjd_mid',mjd_mid))
dict_log.append(('mjd_err',mjd_err))
print 'OBS MJD= ',mjd_mid, '+/- ',mjd_err 

#time on src
min_onsrc=au.timeOnSource(ms_name_final)
print 'Time on Target Source: ',min_onsrc['minutes_on_science'],'minutes'
dict_log.append(('min_on_src',min_onsrc['minutes_on_science']))

#ALMA cycle
alma_cy=au.surmiseCycle(ms_name_final)
print 'Data is from ALMA Cycle ',alma_cy
dict_log.append(('alma_cycle',alma_cy))

#ALMA config
alma_config=au.surmiseConfiguration(ms_name_final)
print 'Data is taken in ALMA configuration: ',alma_config
dict_log.append(('alma_config',alma_config))

#spectral resolution
spec_res=au.effectiveResolution(ms_name_final,0,kms=True)
spec_res2=au.effectiveResolution(ms_name_final,0),kms=False)
print 'Spectral resolution is: ',spec_res,'km/s or ',spec_res2, 'Hz'
dict_log.append(('spec_res_kms',spec_res))
dict_log.append(('spec_res_Hz',spec_res2))

#median PWV
pwv=au.getMedianPWV(ms_name_final)
print 'Median PWV in data set is: ',pwv[0],'+/-',pwv[1],'mm'
dict_log.append(('med_pwv',pwv[0]))
dict_log.append(('med_pwv_err',pwv[1]))

#plot weather cond.
print 'Plotting weather conditions...'
au.plotWeather(ms_name_final,figfile=my_dir+'weather_cond.png')
raw_input('Please press enter when ready to continue.')

#plot bands
print 'Plotting bands...'
au.plotspws(ms_name_final,intents=['*TARGET*'],plotfile=my_dir+'bands.png')
raw_input('Please press enter when ready to continue.')
###########################################

print 'Cleaning up...'
os.system('rm -rf casa*.log')
os.system('rm -rf *.last')
os.system('rm -rf *.png')
print 'Writing user_input log file...'
writeDict(dict_log, my_dir+'user_input_'+obsDate+'_full.logg',str(datetime.datetime.now()))
print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'		