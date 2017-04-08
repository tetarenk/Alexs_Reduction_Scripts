#################################################
#VLA CASA Continuum Reduction Script
#################################################
'''CASA script to be used for the flagging, calibration and imaging of VLA Continuum Data
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) Calibrated Split MS for each band (and base-band) -- [target]_[obsDate]_[band]_[baseband].ms
        (2) Continuum images in each band (and base-band) -- [target]_[obsDate]_[band]_[baseband]_clean1.image(.tt0).pbcor
        (3) File of flux densities from image/UV plane fitting -- imfit_results.txt/uvfit_results.txt
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
       - All output images are also converted to fits format (just append .fits to end of images above)
       - This script is intended to be used with raw data.
       - All input logged in user_input.log.
       - If autoflag is used summary presented in autoflag_log.txt
Written by: Alex J. Tetarenko
Last Updated: Jan 3 2017'''
#Still need to add pol cal!!!!

print '##################################################'
print 'Welcome to Alexs VLA Continuum Reduction Script'
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

#define output directory
my_dir='/mnt/bigdata/tetarenk/test_scripts/vla/'
if not os.path.isdir(my_dir):
	os.system('sudo mkdir '+my_dir)
	os.system('sudo chown ubuntu '+my_dir)
	os.system('sudo chmod -R u+r '+my_dir) 
	os.system('sudo chmod -R u+w '+my_dir)
	os.system('sudo chmod -R u+x '+my_dir)
print 'You have set your output directory to ', my_dir
print 'All output images & intermediate data products are put in this directory.\n'

#param file location
param_dir_file='/home/ubuntu/CASA_reduction_scripts/params_vla.txt'
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
ms_name=data_params.ms_name
ms_name_prefix=ms_name.strip('.ms')
obsDate=data_params.obsDate
target=data_params.target
bands=data_params.bands
spw_bands=data_params.spw_bands
scans=data_params.scans
bitdata=data_params.bitdata
remakems=data_params.remakems
doImage=data_params.doImage
#general image params
use_auto=data_params.use_auto
mythreshold=data_params.mythreshold
myimsize=data_params.myimsize
mycell=data_params.mycell
myniter=data_params.myniter
mynterms=data_params.mynterms
mystokes=data_params.mystokes
outlierfile=data_params.outlierfile
multiscale=data_params.multiscale
robust=data_params.robust
weighting=data_params.weighting
#mask options
mymask=data_params.mymask
#uv fitting
uv_fit=data_params.uv_fit

#write all variables to log dictionary
dict_log.extend([('ms_name_input',ms_name),('ms_name_prefix_input',ms_name_prefix),('obsDate', obsDate),('target',target),\
	('bands_input',bands),('spw_bands_input',spw_bands),('scans_input',scans),('bitdata',bitdata),('remakeMS',remakems),\
	('mythreshold',mythreshold),('myimsize',myimsize),('mycell',mycell),('myniter',myniter),('mynterms',mynterms),\
	('mystokes',mystokes),('mymask',mymask),('uv_fit',uv_fit),('doImage',doImage),('use_auto',use_auto)])
#################################################


#################################################
#Examining and Splitting Data Section
#################################################
#listobs first
if not os.path.isfile(my_dir+obsDate+'_listfile.txt'):
	print 'Making listobs file...'
	listobs(ms_name,listfile=my_dir+obsDate+'_listfile.txt')
else:
	print 'Listobs file already exists'
seelo=raw_input('Do you want to see the listobs? y or n-->')
if seelo=='y':
	os.system('pluma '+my_dir+obsDate+'_listfile.txt &')
	raw_input('Please press enter to continue when you are done.')
else:
	print 'Okay. Moving on.'
seeurl=raw_input('Do you want to open the log file in your web browser?y or n-->')
if seeurl =='y':
	webbrowser.open('http://www.vla.nrao.edu/cgi-bin/oplogs.cgi', new=0, autoraise=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Moving on to making split data sets.'
else:
	print 'Okay. Moving on to making split data sets.'

print '****************************************************************'
print 'You have selected the ',bands,' bands to reduce.'
print '****************************************************************'
print 'We need to split the data into these individual bands to reduce seperately.'

ms_name_list=[]
for i in range(0,len(bands)):
	#make split data sets of different bands if they dont exist
	if not os.path.isdir(my_dir+'data_'+obsDate+'_'+bands[i]+'.ms'):
		print 'Making split MS for ', bands[i],' band...'
		split(vis=ms_name,outputvis=my_dir+'data_'+obsDate+'_'+bands[i]+'.ms',\
			spw=spw_bands[i],datacolumn='data',scan=scans[i])
		ms_name_list.append(my_dir+'data_'+obsDate+'_'+bands[i]+'.ms')
	elif remakems=='T':
		print 'Remaking MS for ',bands[i],' band...'
		os.system('rm -rf '+my_dir+'data_'+obsDate+'_'+bands[i]+'.ms')
		os.system('rm -rf '+my_dir+'data_'+obsDate+'_'+bands[i]+'.ms.flagversions')
		split(vis=ms_name,outputvis=my_dir+'data_'+obsDate+'_'+bands[i]+'.ms',\
			spw=spw_bands[i],datacolumn='data',scan=scans[i])
		ms_name_list.append(my_dir+'data_'+obsDate+'_'+bands[i]+'.ms')
	else:
		print 'Split MS for ',bands[i],' already exists.'
		ms_name_list.append(my_dir+'data_'+obsDate+'_'+bands[i]+'.ms')
dict_log.append(('ms_name_list',ms_name_list))
#################################################

print 'Beginning reduction process...'
#will loop through reduction procedure for all split data sets 
#first define sub-bands
bandL_dict=dict([("L", 'l1'),("S", 's2'),("C", 'c5'),("X", 'x9'),("Ku", 'ku14'),\
	("K", 'k21'),("Ka", 'ka31'),("Q", 'q40')])
bandU_dict=dict([("L", 'l2'),("S", 's3'),("C", 'c7'),("X", 'x11'),("Ku", 'ku17'),\
	("K", 'k26'),("Ka", 'ka36'),("Q", 'q45')])
#will loop through all split data sets now
for kk in range(0,len(ms_name_list)):
	if kk==0:
		firstf='y'
	else:
		firstf='n'
	print 'Starting with ',ms_name_list[kk],' data set.'
	print 'Running listobs. Please examine file for fields ids & spw ranges.'
	os.system('rm -rf '+my_dir+obsDate+'_'+bands[kk]+'_listfile.txt')
	listobs(ms_name_list[i],listfile=my_dir+obsDate+'_'+bands[kk]+'_listfile.txt')
	os.system('pluma '+my_dir+obsDate+'_'+bands[kk]+'_listfile.txt &')
	raw_input('Please press enter to continue when you are done.')
	#define variables specific to each band
	print 'Please enter the following data set specifics:'
	bpf_cal=raw_input('Please enter field id for bandpass cal (if >1 seperate by comma)-->')
	bpf_lst=bpf_cal.split(',')
	second_cal=raw_input('Please enter field id for secondary cal (if >1 seperate by comma)-->')
	second_lst=second_cal.split(',')
	polleak_cal=raw_input('Please enter field id for pol cal (if >1 seperate by comma, if none press enter)-->')
	if polleak_cal=='':
		polleak_lst=[]
	else:
		polleak_lst=polleak_cal.split(',')
	target_id=raw_input('Please enter field id for target (if >1 seperate by comma)-->')
	target_lst=target_id.split(',')
	timerbp=raw_input('Please enter the time range of the bandpass cal scan(s). e.g., 11:43:30.0~11:47:27.0-->')
	field_lst=[]
	[field_lst.append(i) for i in bpf_lst]
	[field_lst.append(i) for i in second_lst]
	[field_lst.append(i) for i in polleak_lst]
	spw_low=raw_input('Please enter lower base-band spw range. e.g., 0~15-->')
	spw_high=raw_input('Please enter upper base-band spw range. e.g., 0~15-->')
	spw_full=spw_low.split('~')[0]+'~'+spw_high.split('~')[1]
	band=bands[kk]
	ms_name=ms_name_list[kk]
	ms_name_prefix=ms_name.strip('.ms')
	band_low=bandL_dict[band]
	band_high=bandU_dict[band]
	date=obsDate+'_'+band
	ant_plot=my_dir+'antennas_'+date+'.png'
	ant_cal=my_dir+date+'.antpos'
	opac_cal=my_dir+date+'.opac'
	gc_cal=my_dir+date+'.gaincurve'
	req_cal=my_dir+date+'.req'
	cal_table_prefix=my_dir+target+'_'+date
	split_low=cal_table_prefix+'_'+band_low+'.ms'
	split_high=cal_table_prefix+'_'+band_high+'.ms'
	split_full=cal_table_prefix+'_comb.ms'
	#add to log dictionary
	dict_log.append((ms_name_prefix+'_bpf_cal',bpf_cal))
	dict_log.append((ms_name_prefix+'_second_cal',second_cal))
	dict_log.append((ms_name_prefix+'_polleak_cal',polleak_cal))
	dict_log.append((ms_name_prefix+'_target_id',target_id))
	dict_log.append((ms_name_prefix+'_time_range_bp',timerbp))
	dict_log.append((ms_name_prefix+'_spw_low',spw_low))
	dict_log.append((ms_name_prefix+'_spw_high',spw_high))
	dict_log.append((ms_name_prefix+'_spw_full',spw_full))
	dict_log.append((ms_name_prefix+'_band',band))
	dict_log.append((ms_name_prefix+'_band_low',band_low))
	dict_log.append((ms_name_prefix+'_band_high',band_high))
	#########################################
	#Reference antenna and initial flagging
	#########################################
	#plot antenna positions to determne reference antenna near centre of array that 
	#is present throughout whole observation and has no problems
	print 'Plotting antenna positions. Please choose a reference antenna near centre of array.'
	os.system('rm -rf '+ant_plot)
	plotants(vis=ms_name,figfile=ant_plot)
	ref_ant=raw_input('Please enter reference antenna. e.g., ea02-->')
	dict_log.append((ms_name_prefix+'_ref_ant',ref_ant))
	flag_done=raw_input('Is the data set already flagged?y or n-->')
	dict_log.append((ms_name_prefix+'_flag_done',flag_done))
	if flag_done=='n':
		#initial flagging
		print 'Starting Initial Flagging...'
		flag_dummy=raw_input('Do you want to flag the dummy scans?y or n-->')
		dict_log.append((ms_name_prefix+'_flag_dummy',flag_dummy))
		if flag_dummy=='y':
			dum_scan=raw_input('Please enter dummy scan ids. e.g., 1,2-->')
			print 'Flagging Dummy scan(s) ', dum_scan,' ...'
			dict_log.append((ms_name_prefix+'_dum_scan',dum_scan))
			flagdata(vis=ms_name, flagbackup=T, mode='manual', scan=dum_scan)
		else:
			print 'Dummy scans not flagged.'
		flag_shad=raw_input('Do you want to flag the shadow data?y or n-->')
		dict_log.append((ms_name_prefix+'_flag_shad',flag_shad))
		if flag_shad=='y':
			print 'Flagging Shadow data...'
			flagdata(vis=ms_name, mode='shadow')
		else:
			print 'Shadow data not flagged.'
		flag_zero=raw_input('Do you want to flag the zero amp data?y or n-->')
		dict_log.append((ms_name_prefix+'_flag_zero',flag_zero))
		if flag_zero=='y':
			print 'Flagging zero amp data...'
			flagdata(vis=ms_name, autocorr=T)
			flagdata(vis=ms_name,mode='clip',clipzeros=True)
		else:
			print 'Zero amp data not flagged.'
		flag_end=raw_input('Do you want to flag the end channels?y or n-->')
		dict_log.append((ms_name_prefix+'_flag_end_chan',flag_end))
		if flag_end=='y':
			beg=raw_input('Please enter beginning channels to flag. e.g., 0~3-->')
			endd=raw_input('Please enter end channels to flag. e.g., 60~63-->')
			print 'Flagging beginning channels ',beg,' and end channels ',endd,' ...'
			dict_log.append((ms_name_prefix+'_flag_beg_chan',beg))
			dict_log.append((ms_name_prefix+'_flag_end_chan',endd))
			flagdata(vis=ms_name,spw=spw_full+':'+beg,field='',antenna='')
			flagdata(vis=ms_name,spw=spw_full+':'+endd,field='',antenna='')
		else:
			print 'End channels not flagged.'
		flag_badant=raw_input('Do you want to flag the bad antennas from log file?y or n-->')
		dict_log.append((ms_name_prefix+'_flag_bad_ant',flag_badant))
		if flag_badant=='y':
			bada=raw_input('Please enter bad ants to flag. e.g. ea01,ea02-->')
			print 'Flagging bad ants ', bada,' ...'
			dict_log.append((ms_name_prefix+'_badant',bada))
			flagdata(vis=ms_name, flagbackup=T, mode='manual', antenna=bada)
		else:
			print 'No bad ants flagged.'
		print 'Flagging samples at start of each scan...'
		quakint=raw_input('Please enter integration time in seconds. e.g., 3-->')
		dict_log.append((ms_name_prefix+'_flag_quack_int_time',quakint))
		flagdata(vis=ms_name, mode='quack', quackinterval=float(quakint), quackmode='beg',quackincrement=True)
		########################################

		########################################
		#detailed flagging
		########################################
		print 'Examining Data for detailed flagging...'
		print 'Look for obvious bad data.'
		print 'Keep track of ants/spws/channels to flag. You will be prompted for their values after plotting.'
		print '(1) Amp vs time...'
		plotms(vis=ms_name,xaxis="time",yaxis="amp",coloraxis="field",iteraxis="antenna",\
		 avgtime='60s')
		raw_input('Please press enter when ready to continue.')
		autoflag=raw_input('Is the RFI bad enough that you need to do autoflag?y or n-->')
		dict_log.append((ms_name_prefix+'_flag_autoflag',autoflag))
		########################################
		#Auto-flagging (if needed)
		#########################################
		if autoflag=='y':
			print 'Preparing data for autoflag procedure...'
			print 'Solving for baseline solutions...'
			#baseline solutions
			os.system('rm -rf '+ ant_cal)
			os.system('rm -rf '+ gc_cal)
			os.system('rm -rf '+ opac_cal)
			os.system('rm -rf '+ req_cal)
			gencal(vis=ms_name, caltable=ant_cal,caltype='antpos')
			if band in ['Ku','K','Ka','Q']:
				print 'Solving for gain-elevation curves...'
				#gain-elev curves
				gencal(ms_name,caltable=gc_cal,caltype='gc')
				print 'Solving for opacity corrections...'
				#opacity corrections
				clearstat()
				myTau = plotweather(vis=ms_name, doPlot=T)
				gencal(vis=ms_name,caltable=opac_cal, caltype='opac',spw=spw_full,parameter=myTau)
			if bitdata==3:
				print 'Solving for requantizer gains...'
				#requantizer gains
				gencal(ms_name,caltable=req_cal,caltype='rq')
			print 'Hanning smooth of data...'
			#hanning smooth data
			os.system('rm -rf '+ ms_name_prefix+'_hs.ms')
			hanningsmooth(vis=ms_name,datacolumn='all', outputvis=ms_name_prefix+'_hs.ms')
			ms_name=ms_name_prefix+'_hs.ms'
			ms_name_prefix=ms_name.strip('.ms')
			#plot hanning smoothed data
			print 'Plotting Hanning smoothed data of target...'
			plotms(vis=ms_name, field=target_id, antenna=ref_ant, spw=spw_full,\
				xaxis='freq', yaxis='amp', coloraxis='spw', symbolshape = 'circle',\
				 correlation='RR,LL')
			raw_input('Please press enter when ready to continue.')
			print 'Using secondary cal to do preliminary bandpass...'
			print 'Plotting Amp vs channel...'
			print 'Choose 3 channel range for each spw free of RFI.'
			plotms(vis=ms_name, field=second_cal,antenna=ref_ant, \
				xaxis='channel', yaxis='amp', iteraxis='spw',yselfscale=True, \
				correlation='RR,LL', symbolshape='circle')
			raw_input('Please press enter when ready to continue.')
			spw_string_rfifree=raw_input('Enter spw string for RFI free channels. e.g., 0:30~33,1:30~33-->')
			dict_log.append((ms_name_prefix+'_autoflag_rfi_free_chan',spw_string_rfifree))
			print 'Performing gaincal for initial phase solutions...'
			os.system('rm -rf '+cal_table_prefix+'.initPhautoF')
			gaincal(vis=ms_name, caltable=cal_table_prefix+'.initPhautoF',\
			 field=second_cal, solint='int', spw=spw_string_rfifree, refant=ref_ant, \
			 minblperant=3, minsnr=3.0, calmode='p')
			print 'Doing bandpass...'
			os.system('rm -rf '+cal_table_prefix+'.initBPautoF')
			bandpass(vis=ms_name, caltable=cal_table_prefix+'.initBPautoF', \
				field=second_cal, solint='inf', combine='scan', refant=ref_ant,\
				 minblperant=3, minsnr=10.0, gaintable=[cal_table_prefix+'.initPhautoF'],\
				  interp=['', 'nearest'], solnorm=False)
			print 'Inspecting solutions...'
			print 'Look for spws which have so much RFI they are a lost cause.'
			print 'Amp vs freq...'
			plotcal(caltable=cal_table_prefix+'.initBPautoF', xaxis='freq', yaxis='amp',\
			 iteration='antenna', subplot=331)
			raw_input('Please press enter when ready to continue.')
			print 'Phase vs frequency...'
			plotcal(caltable=cal_table_prefix+'.initBPautoF', xaxis='freq', yaxis='phase',\
			 iteration='antenna', subplot=331)
			raw_input('Please press enter when ready to continue.')
			spw_string_rfiflag=raw_input('Enter spw/channels string to flag (enter if none). e.g., 8,9:20~23,10:32-->')
			dict_log.append((ms_name_prefix+'_autoflag_rfi_flag',spw_string_rfiflag))
			if spw_string_rfiflag=='':
				print 'No flagging done.'
			else:
				flagdata(vis=ms_name, spw=spw_string_rfiflag)
			print 'Applying calibration...'
			applycal(vis=ms_name, gaintable=[cal_table_prefix+'.initBPautoF'], calwt=False)
			print 'Beginning autoflag...'
			print 'Will test first then apply after.'
			cont_af='y'
			while cont_af=='y':
				fdev,tdev=raw_input('Please enter deviations of amp (sigma) in freq and time to flag. e.g., 4 4-->').split(' ')
				flagdata(vis=ms_name, mode='rflag', field=second_cal,spw='', \
					datacolumn='corrected', freqdevscale=int(fdev), timedevscale=int(tdev), \
					action='calculate', display='both', flagbackup=False)
				raw_input('Please press enter when ready to continue.')
				print 'Extending flags in freq and time...'
				print 'Growtime means a channel is flagged if > X% of times already flagged.'
				print 'Growfreq means entire spectrum for an integration is flagged if >Y% of channels already flagged'
				growf,growt=raw_input('Please enter growfreq and growtime params. e.g., 90.0 50.0-->').split(' ')
				flagdata(vis=ms_name, mode='extend', spw='', growtime=float(growt), \
					growfreq=float(growf), action='calculate', display='data', flagbackup=False)
				raw_input('Please press enter when ready to continue.')
				cont_af=raw_input('Do you want to try different autoflag param values?y or n-->')
			print 'Now applying autoflag selections...'
			fdev,tdev,growf,growt=raw_input('Please enter final fdev,tdev,growf,growt. e.g., 4 4 90.0 50.0-->').split(' ')
			dict_log.append((ms_name_prefix+'_autoflag_fdev',fdev))
			dict_log.append((ms_name_prefix+'_autoflag_tdev',tdev))
			dict_log.append((ms_name_prefix+'_autoflag_growf',growf))
			dict_log.append((ms_name_prefix+'_autoflag_growt',growt))
			flagdata(vis=ms_name, mode='rflag', field='',spw='', datacolumn='corrected',\
			 freqdevscale=int(fdev), timedevscale=int(tdev), action='apply')
			flagdata(vis=ms_name, mode='extend', spw='', growtime=float(growt), growfreq=float(growf), action='apply')
			print 'Writing autoflag log...'
			flagInfo = flagdata(vis=ms_name, mode='summary')
			openf=open(my_dir+'autoflag_log.txt','w')
			openf.write('Autoflag stats for '+str(band)+' band:')
			for i in flagInfo['field'].keys():
				print("\n %2.1f%% of %s are flagged.\n" % (100.0*flagInfo['field'][i]['flagged']/flagInfo['field'][i]['total'],i))
				openf.write(("\n %2.1f%% of %s are flagged.\n" % (100.0*flagInfo['field'][i]['flagged']/flagInfo['field'][i]['total'],i)))
			print("Spectral windows are flagged as follows:")
			for spw in range(0,int(spw_full.split('~')[1])+1):
				print("SPW %s: %2.1f%%" % (spw, 100.0 * flagInfo['spw'][str(spw)]['flagged'] / flagInfo['spw']\
					[str(spw)]['total']))
				openf.write(("\nSPW %s: %2.1f%%" % (spw, 100.0 * flagInfo['spw'][str(spw)]['flagged'] / flagInfo['spw']\
					[str(spw)]['total'])))
			openf.close()
			print 'Autoflag finished. Back to normal flagging procedure...'
			#########################################
		else:
			print 'Autoflag not selected.'
		print '(2) Baselines versus antenna to look for bad antennas,bad channels/spws...'
		lastf,fant=raw_input('Please enter last field id and first antenna. e.g., 2 ea01-->').split(' ')
		dict_log.append((ms_name_prefix+'_flag_lastfield',lastf))
		dict_log.append((ms_name_prefix+'_flag_firstant',fant))
		plotms(vis=ms_name,field=lastf,spw='', antenna=fant,correlation='RR,LL',xaxis='antenna2',yaxis='amp')
		raw_input('Please press enter when ready to continue.')
		print '(3) Bandpass for reference antenna (second_cal)...'
		print 'Look for bad spws and RFI spikes'
		print 'Amp first...'
		plotms(vis=ms_name,field=second_cal, spw='',antenna=ref_ant+'&*',correlation='RR,LL',\
			xaxis='frequency',yaxis='amp',iteraxis='baseline')
		raw_input('Please press enter when ready to continue.')
		print 'Phase next...'
		print 'RR correlation...'
		plotms(vis=ms_name,field=second_cal, spw='',antenna=ref_ant+'&*',correlation='RR',\
			xaxis='frequency',yaxis='phase',iteraxis='baseline')
		raw_input('Please press enter when ready to continue.')
		print 'LL correlation...'
		plotms(vis=ms_name,field=second_cal, spw='',antenna=ref_ant+'&*',correlation='LL',\
			xaxis='frequency',yaxis='phase',iteraxis='baseline')
		raw_input('Please press enter when ready to continue.')
		print '(4A) Baselines (all not just to ref antenna) vs amp for target...'
		plotms(vis=ms_name,field=target_id, xaxis='baseline', yaxis='amp', spw='', iteraxis='spw',\
			correlation='RR,LL', coloraxis='antenna1', symbolshape = 'circle')
		raw_input('Please press enter when ready to continue.')
		print '(4B) Amp vs uvdist for target...'
		plotms(vis=ms_name, field=target_id, xaxis='uvdist', yaxis='amp', spw='', iteraxis='spw',\
			correlation='RR,LL', coloraxis='antenna1', symbolshape = 'circle')
		raw_input('Please press enter when ready to continue.')
		print '(5) Data stream plot...'
		plotms(vis=ms_name,field='',correlation='RR,LL',timerange='',\
			antenna=fant,spw='',xaxis='time',yaxis='antenna2',\
			plotrange=[-1,-1,0,28],coloraxis='field')
		raw_input('Please press enter when ready to continue.')
		print 'Flagging...'
		badasf=raw_input('Please enter bad ant,spw,and field to flag (enter if none). e.g., ea10,ea12;5:4~9;3 ;5;3-->').split(' ')
		dict_log.append((ms_name_prefix+'_flag_antspwfield',badasf))
		if badasf=='':
			print 'Nothing to flag.'
		else:
			print 'Flagging selected data.'
			for i in range(0,len(badasf)):
				strg=badasf[i].split(';')
				flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=strg[0],spw=strg[1],field=strg[2])
		print 'Final check of flagged data...'
		plotms(vis=ms_name,field=second_cal,spw='', antenna=ref_ant,correlation='RR,LL',xaxis='frequency',yaxis='amp')
		raw_input('Please press enter when ready to continue.')
		flag_again=raw_input('Do you need to do more flagging? y or n-->')
		while flag_again=='y':
			count_f=1
			badasf2=raw_input('Please enter bad ant,spw,and field to flag (enter if none). e.g., ea10,ea12;5:4~9;3 ;5;3-->').split(' ')
			dict_log.append((ms_name_prefix+'_flag_antspwfield_'+str(count_f),badasf2))
			if badasf2=='':
				print 'Nothing to flag.'
			else:
				print 'Flagging selected data.'
				for i in range(0,len(badasf2)):
					strg2=badasf2[i].split(';')
					flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=strg2[0],spw=strg2[1],field=strg2[2])
			print 'Plotting...'
			plotms(vis=ms_name,field=second_cal,spw='', antenna=ref_ant,correlation='RR,LL',xaxis='frequency',yaxis='amp')
			raw_input('Please press enter when ready to continue.')
			count_f=count_f+1
			flag_again=raw_input('Do you need to do more flagging? y or n-->')
	else:
		autoflag='n'
		dict_log.append((ms_name_prefix+'_flag_autoflag',autoflag))
	intera=raw_input('Flagging is finished. Do you want to do interactive calibration?y or n-->')
	dict_log.append((ms_name_prefix+'_interactive',intera))
	if intera=='n':
		print 'You have chosen to not do interactive calibration.'
		print 'No plots will be made and no additional flagging is required.'
		print 'Please go do something else for a while.'
		print 'You will not be prompted until imaging.'
		print '*********************************************************************'
	else:
		print 'The rest of the script is interactive. Please stay by the computer.'
		print '*********************************************************************'
	########################################

	########################################
	##baseline,gain curve, opacity solutions
	#(if no autoflag)
	########################################
	if autoflag=='n':
		print 'Solving for baseline solutions...'
		#baseline solutions
		os.system('rm -rf '+ant_cal)
		os.system('rm -rf '+gc_cal)
		os.system('rm -rf '+opac_cal)
		os.system('rm -rf '+req_cal)
		gencal(vis=ms_name, caltable=ant_cal,caltype='antpos')
		if band in ['Ku','K','Ka','Q']:
			print 'Solving for gain-elevation curves...'
			#gain-elev curves
			gencal(ms_name,caltable=gc_cal,caltype='gc')
			print 'Solving for opacity corrections...'
			#opacity corrections
			clearstat()
			myTau = plotweather(vis=ms_name, doPlot=T)
			gencal(vis=ms_name,caltable=opac_cal, caltype='opac',spw=spw_full,parameter=myTau)
		if bitdata==3:
			print 'Solving for requantizer gains...'
			#requantizer gains
			gencal(ms_name,caltable=req_cal,caltype='rq')
	#####################################

	#make list of pre gain tables to be applied in tandem below
	gt_lst_low=[]
	gf_lst_low=[]
	gi_lst_low=[]
	gt_lst_high=[]
	gf_lst_high=[]
	gi_lst_high=[]
	if os.path.isdir(ant_cal):
		gt_lst_low.append(ant_cal)
		gf_lst_low.append('')
		gi_lst_low.append('')
		gt_lst_high.append(ant_cal)
		gf_lst_high.append('')
		gi_lst_high.append('')
	if os.path.isdir(gc_cal):
		gt_lst_low.append(gc_cal)
		gf_lst_low.append('')
		gi_lst_low.append('')
		gt_lst_high.append(gc_cal)
		gf_lst_high.append('')
		gi_lst_high.append('')
	if os.path.isdir(opac_cal):
		gt_lst_low.append(opac_cal)
		gf_lst_low.append('')
		gi_lst_low.append('')
		gt_lst_high.append(opac_cal)
		gf_lst_high.append('')
		gi_lst_high.append('')
	if os.path.isdir(req_cal):
		gt_lst_low.append(req_cal)
		gf_lst_low.append('')
		gi_lst_low.append('')
		gt_lst_high.append(req_cal)
		gf_lst_high.append('')
		gi_lst_high.append('')

	#####################################
	#fluxscale
	#####################################
	#list models available
	print 'Setting flux scale...'
	print 'Listing flux models available.'
	setjy(vis=ms_name, listmodels=T)
	print 'Set fluxscale for lower base-band...'
	flux_mod_low=raw_input('Please enter flux model for lower base-band. e.g., 3C48_C.im-->')
	dict_log.append((ms_name_prefix+'_fluxmodlsb',flux_mod_low))
	setjy(vis=ms_name,field=bpf_cal,standard='Perley-Butler 2013',
      model=flux_mod_low,usescratch=False,scalebychan=True,spw=spw_low)
	print 'Set fluxscale for upper base-band...'
	flux_mod_high=raw_input('Please enter flux model for upper base-band. e.g., 3C48_C.im-->')
	dict_log.append((ms_name_prefix+'_fluxmodusb',flux_mod_high))
	setjy(vis=ms_name,field=bpf_cal,standard='Perley-Butler 2013',
      model=flux_mod_high,usescratch=False,scalebychan=True,spw=spw_low)
	########################################

	########################################
	##initial phase cal before bandpass
	########################################
	os.system('rm -rf '+cal_table_prefix+'_phase_int_all.cal')
	os.system('rm -rf '+cal_table_prefix+'_phase_int_all2.cal')
	os.system('rm -rf '+cal_table_prefix+'_phase_int_bp.cal')
	os.system('rm -rf '+cal_table_prefix+'_phase_int_bp2.cal')
	print 'Initial phase cal before bandpass...'
	if intera=='y':
		print 'First look for RFI free channels in lower base-band...'
		plotms(vis=ms_name, spw=spw_low, antenna=ref_ant, xaxis='freq', yaxis='amp',iteraxis='field',\
			correlation='RR,LL', coloraxis='spw', symbolshape = 'circle')
		raw_input('Please press enter when ready to continue.')
		rfifree_low=raw_input('Please enter rfi free channels for lower base-band. e.g. 0~7:27~36-->')
		dict_log.append((ms_name_prefix+'_initphase_rfifreelsb',rfifree_low))
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_all.cal',field=",".join(field_lst),\
		 refant=ref_ant, spw=rfifree_low,gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_low)
		print 'Plotting solutions...'
		plotcal(caltable=cal_table_prefix+'_phase_int_all.cal',xaxis='time',yaxis='phase',\
			poln='R',iteration='antenna',plotrange=[-1,-1,-180,180],field='',antenna='')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_phase_int_all.cal',xaxis='time',yaxis='phase',\
			poln='L',iteration='antenna',plotrange=[-1,-1,-180,180],field='')
		raw_input('Please press enter when ready to continue.')
		badantsf3a=raw_input('Please enter bad ants to flag (enter if none). e.g., ea10,ea16-->')
		dict_log.append((ms_name_prefix+'_initphase_flag_badantlsb',badantsf3a))
		if badantsf3a=='':
				print 'No antennas to flag.'
		else:
			print 'Flagging selected ants.'
			flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=badantsf3a,spw=spw_low)
		print 'Look for RFI free channels in upper base-band...'
		plotms(vis=ms_name, spw=spw_high, antenna=ref_ant, xaxis='freq', yaxis='amp',iteraxis='field',\
			correlation='RR,LL', coloraxis='spw', symbolshape = 'circle')
		raw_input('Please press enter when ready to continue.')
		rfifree_high=raw_input('Please enter rfi free channels for upper base-band. e.g. 8~15:27~36-->')
		dict_log.append((ms_name_prefix+'_initphase_rfifreeusb',rfifree_high))
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_all2.cal', field=",".join(field_lst),\
		 refant=ref_ant, spw=rfifree_high,gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_high)
		print 'Plotting solutions...'
		plotcal(caltable=cal_table_prefix+'_phase_int_all2.cal',xaxis='time',yaxis='phase',\
			poln='R',iteration='antenna',plotrange=[-1,-1,-180,180],field='',antenna='')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_phase_int_all2.cal',xaxis='time',yaxis='phase',\
			poln='L',iteration='antenna',plotrange=[-1,-1,-180,180],field='')
		raw_input('Please press enter when ready to continue.')
		badantsf3b=raw_input('Please enter bad ants to flag (enter if none). e.g., ea10,ea16-->')
		dict_log.append((ms_name_prefix+'_initphase_flag_badantusb',badantsf3b))
		if badantsf3b=='':
				print 'No antennas to flag.'
		else:
			print 'Flagging selected ants.'
			flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=badantsf3b,spw=spw_high)
		print 'Doing only Bandpass cal solution...'
		print 'Lower base-band...'
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_bp.cal',field=bpf_cal, refant=ref_ant,\
		 spw=rfifree_low,gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_low)
		print 'Plotting solutions...'
		plotcal(caltable=cal_table_prefix+'_phase_int_bp.cal',xaxis='time',yaxis='phase',\
			poln='R',iteration='antenna',plotrange=[-1,-1,-180,180],field=bpf_cal,antenna='',\
			timerange=timerbp)
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_phase_int_bp.cal',xaxis='time',yaxis='phase',\
			poln='L',iteration='antenna',plotrange=[-1,-1,-180,180],field=bpf_cal,antenna='',\
			timerange=timerbp)
		raw_input('Please press enter when ready to continue.')
		print 'Upper base-band...'
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_bp2.cal', field=bpf_cal, refant=ref_ant,\
		 spw=rfifree_high,gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_high)
		print 'Plotting solutions...'
		plotcal(caltable=cal_table_prefix+'_phase_int_bp2.cal',xaxis='time',yaxis='phase',\
			poln='R',iteration='antenna',plotrange=[-1,-1,-180,180],field=bpf_cal,antenna='',\
			timerange=timerbp)
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_phase_int_bp2.cal',xaxis='time',yaxis='phase',\
			poln='L',iteration='antenna',plotrange=[-1,-1,-180,180],field=bpf_cal,antenna='',\
			timerange=timerbp)
		raw_input('Please press enter when ready to continue.')
	else:
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_all.cal',field=",".join(field_lst),\
		 refant=ref_ant, spw=spw_low+':27~36',gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_low)
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_all2.cal', field=",".join(field_lst),\
		 refant=ref_ant, spw=spw_high+':27~36',gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_high)
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_bp.cal',field=bpf_cal, refant=ref_ant,\
		 spw=spw_low+':27~36',gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_low)
		gaincal(vis=ms_name, caltable=cal_table_prefix+'_phase_int_bp2.cal', field=bpf_cal, refant=ref_ant,\
		 spw=spw_high+':27~36',gaintype='G',calmode='p', solint='int', minsnr=2.0,gaintable=gt_lst_high)

	
	gt_lst_high.append(cal_table_prefix+'_phase_int_bp2.cal')
	gt_lst_low.append(cal_table_prefix+'_phase_int_bp.cal')
	########################################

	########################################
	#Delay calibration
	########################################
	os.system('rm -rf '+cal_table_prefix+'.K0')
	os.system('rm -rf '+cal_table_prefix+'_2.K0')
	print 'Delay calibration...'
	if intera=='y':
		print 'Lower base-band...'
		gaincal(vis=ms_name,caltable=cal_table_prefix+'.K0', field=bpf_cal,\
			refant=ref_ant,spw=spw_low,gaintype='K', solint='inf',combine='scan',\
			minsnr=2, gaintable=gt_lst_low)
		print 'Plotting solutions...'
		plotcal(caltable=cal_table_prefix+'.K0',xaxis='antenna',yaxis='delay',figfile='plotcal_K0-delay.png')
		raw_input('Please press enter when ready to continue.')
		badantsf4a=raw_input('Please enter bad ants to flag (enter if none). e.g., ea10,ea16-->')
		dict_log.append((ms_name_prefix+'_delay_flag_badantlsb',badantsf4a))
		if badantsf4a=='':
				print 'No antennas to flag.'
		else:
			print 'Flagging selected ants.'
			flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=badantsf4a,spw=spw_low)
		print 'Upper base-band...'
		gaincal(vis=ms_name,caltable=cal_table_prefix+'_2.K0', field=bpf_cal,\
			refant=ref_ant,spw=spw_high,gaintype='K', solint='inf',combine='scan',\
			minsnr=2, gaintable=gt_lst_high)
		print 'Plotting solutions...'
		plotcal(caltable=cal_table_prefix+'_2.K0',xaxis='antenna',yaxis='delay',figfile='plotcal_K02-delay.png')
		raw_input('Please press enter when ready to continue.')
		badantsf4b=raw_input('Please enter bad ants to flag (enter if none). e.g., ea10,ea16-->')
		dict_log.append((ms_name_prefix+'_delay_flag_badantusb',badantsf4b))
		if badantsf4b=='':
				print 'No antennas to flag.'
		else:
			print 'Flagging selected ants.'
			flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=badantsf4b,spw=spw_high)
	else:
		gaincal(vis=ms_name,caltable=cal_table_prefix+'.K0', field=bpf_cal,\
			refant=ref_ant,spw=spw_low,gaintype='K', solint='inf',combine='scan',\
			minsnr=2, gaintable=gt_lst_low)
		gaincal(vis=ms_name,caltable=cal_table_prefix+'_2.K0', field=bpf_cal,\
			refant=ref_ant,spw=spw_high,gaintype='K', solint='inf',combine='scan',\
			minsnr=2, gaintable=gt_lst_high)

	gt_lst_low.append(cal_table_prefix+'.K0')
	gt_lst_high.append(cal_table_prefix+'_2.K0')
	gf_lst_low.append('')
	gi_lst_low.append('')
	gf_lst_high.append('')
	gi_lst_high.append('')
	########################################

	########################################
	#Bandpass calibration
	########################################
	os.system('rm -rf '+cal_table_prefix+'.B0')
	os.system('rm -rf '+cal_table_prefix+'_2.B0')
	print 'Bandpass calibration...'
	if intera=='y':
		print 'Lower base-band...'
		bandpass(vis=ms_name,caltable=cal_table_prefix+'.B0',field=bpf_cal,spw=spw_low,\
			refant=ref_ant,solnorm=True,combine='scan', solint='inf',bandtype='B',\
			gaintable=gt_lst_low)
		print 'Plotting solutions...'
		plotcal(caltable= cal_table_prefix+'.B0',poln='R', xaxis='chan',yaxis='amp',field= bpf_cal,\
			subplot=221, iteration='antenna',figfile='plotcal_BP-B0-R-amp.png')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable= cal_table_prefix+'.B0',poln='L', xaxis='chan',yaxis='amp',field= bpf_cal,\
			subplot=221, iteration='antenna',figfile='plotcal_BP-B0-L-amp.png')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable= cal_table_prefix+'.B0',poln='R', xaxis='chan',yaxis='phase',field= bpf_cal,\
			subplot=221, iteration='antenna',plotrange=[-1,-1,-180,180],figfile='plotcal_BP-B0-R-phase.png')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable= cal_table_prefix+'.B0',poln='L', xaxis='chan',yaxis='phase',field= bpf_cal,\
			subplot=221, iteration='antenna',plotrange=[-1,-1,-180,180],figfile='plotcal_BP-B0-L-phase.png')
		raw_input('Please press enter when ready to continue.')
		print 'Upper base-band...'
		bandpass(vis=ms_name,caltable=cal_table_prefix+'_2.B0',field=bpf_cal,spw=spw_high,\
			refant=ref_ant,solnorm=True,combine='scan', solint='inf',bandtype='B',\
			gaintable=gt_lst_high)
		print 'Plotting solutions...'
		plotcal(caltable= cal_table_prefix+'_2.B0',poln='R', xaxis='chan',yaxis='amp',field= bpf_cal,\
			subplot=221, iteration='antenna',figfile='plotcal_BP-B02-R-amp.png')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable= cal_table_prefix+'_2.B0',poln='L', xaxis='chan',yaxis='amp',field= bpf_cal,\
			subplot=221, iteration='antenna',figfile='plotcal_BP-B02-L-amp.png')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable= cal_table_prefix+'_2.B0',poln='R', xaxis='chan',yaxis='phase',field= bpf_cal,\
			subplot=221, iteration='antenna',plotrange=[-1,-1,-180,180],figfile='plotcal_BP-B02-R-phase.png')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable= cal_table_prefix+'_2.B0',poln='L', xaxis='chan',yaxis='phase',field= bpf_cal,\
			subplot=221, iteration='antenna',plotrange=[-1,-1,-180,180],figfile='plotcal_BP-B02-L-phase.png')
		raw_input('Please press enter when ready to continue.')
	else:
		bandpass(vis=ms_name,caltable=cal_table_prefix+'.B0',field=bpf_cal,spw=spw_low,\
			refant=ref_ant,solnorm=True,combine='scan', solint='inf',bandtype='B',\
			gaintable=gt_lst_low)
		bandpass(vis=ms_name,caltable=cal_table_prefix+'_2.B0',field=bpf_cal,spw=spw_high,\
			refant=ref_ant,solnorm=True,combine='scan', solint='inf',bandtype='B',\
			gaintable=gt_lst_high)

	gt_lst_low.append(cal_table_prefix+'.B0')
	gt_lst_high.append(cal_table_prefix+'_2.B0')
	gf_lst_low.append('')
	gi_lst_low.append('')
	gf_lst_high.append('')
	gi_lst_high.append('')
	########################################


	########################################
	#gain cal
	########################################
	os.system('rm -rf '+cal_table_prefix+'.G1')
	os.system('rm -rf '+cal_table_prefix+'_2.G1')
	gt_lst_low.remove(cal_table_prefix+'_phase_int_bp.cal')
	gt_lst_high.remove(cal_table_prefix+'_phase_int_bp2.cal')
	print 'Gain cal...'
	print 'First BP/flux cal for both base-bands...'
	for i in range(0,len(bpf_lst)):
		gaincal(vis=ms_name,caltable=cal_table_prefix+'.G1',field=bpf_lst[i],spw=spw_low,\
			solint='inf',refant=ref_ant,gaintype='G',calmode='ap',solnorm=F, gaintable=gt_lst_low)
		gaincal(vis=ms_name,caltable=cal_table_prefix+'_2.G1',field=bpf_lst[i],spw=spw_high,\
			solint='inf',refant=ref_ant,gaintype='G',calmode='ap',solnorm=F, gaintable=gt_lst_high)	
	print 'Append second cals to same table...'
	for i in range(0,len(second_lst)):
		gaincal(vis=ms_name,caltable=cal_table_prefix+'.G1',field=second_lst[i],spw=spw_low,\
			solint='inf',refant=ref_ant,gaintype='G',calmode='ap', gaintable=gt_lst_low,append=True)
		gaincal(vis=ms_name,caltable=cal_table_prefix+'_2.G1',field=second_lst[i],spw=spw_high,\
			solint='inf',refant=ref_ant,gaintype='G',calmode='ap', gaintable=gt_lst_high,append=True)
	if len(polleak_lst)>0:
		print 'Append polleak cals to same table...'
		for i in range(0,len(polleak_lst)):
			gaincal(vis=ms_name,caltable=cal_table_prefix+'.G1',field=polleak_lst[i],spw=spw_low,\
				solint='inf',refant=ref_ant,gaintype='G',calmode='ap', gaintable=gt_lst_low,append=True)
			gaincal(vis=ms_name,caltable=cal_table_prefix+'_2.G1',field=polleak_lst[i],spw=spw_high,\
				solint='inf',refant=ref_ant,gaintype='G',calmode='ap', gaintable=gt_lst_high,append=True)
	if intera=='y':
		print 'Plotting solutions...'
		print 'Lower base-band phases...'
		plotcal(caltable=cal_table_prefix+'.G1',xaxis='time',yaxis='phase',poln='R',\
			plotrange=[-1,-1,-180,180],iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'.G1',xaxis='time',yaxis='phase',poln='L',\
			plotrange=[-1,-1,-180,180],iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		print 'Lower base-band amps...'
		plotcal(caltable=cal_table_prefix+'.G1',xaxis='time',yaxis='amp',poln='R',iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'.G1',xaxis='time',yaxis='amp',poln='L',iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		print 'Upper base-band phases...'
		plotcal(caltable=cal_table_prefix+'_2.G1',xaxis='time',yaxis='phase',poln='R',\
			plotrange=[-1,-1,-180,180],iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_2.G1',xaxis='time',yaxis='phase',poln='L',\
			plotrange=[-1,-1,-180,180],iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		print 'Upper base-band amps...'
		plotcal(caltable=cal_table_prefix+'_2.G1',xaxis='time',yaxis='amp',poln='R',iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_2.G1',xaxis='time',yaxis='amp',poln='L',iteration='antenna')
		raw_input('Please press enter when ready to continue.')
		print 'Check phase stability of reference antenna (diff between R and L stable over time)...'
		print 'Lower base-band...'
		plotcal(caltable=cal_table_prefix+'.G1', xaxis='time', yaxis='phase',poln='/',\
			plotrange=[-1,-1,-180,180],iteration='spw')
		raw_input('Please press enter when ready to continue.')
		print 'Upper base-band...'
		plotcal(caltable=cal_table_prefix+'_2.G1', xaxis='time', yaxis='phase',poln='/',\
			plotrange=[-1,-1,-180,180],iteration='spw')
		raw_input('Please press enter when ready to continue.')
	#gt_lst_low.append(cal_table_prefix+'.G1')
	#gt_lst_high.append(cal_table_prefix+'_2.G1')
	########################################

	########################################
	#Polarization calibration (TBD)
	########################################
	pol_calib='n'#raw_input('Do you wish to do poalrization calibration? y or n-->')
	dict_log.append((ms_name_prefix+'_pol_calib',pol_calib))
	if pol_calib=='y':
		print 'Polarization Calibration...'
		foo=2
	else:
		print 'No polarization calibration done.'
	########################################

	########################################
	#Scaling amp gains
	########################################
	os.system('rm -rf '+cal_table_prefix+'.fluxscale1')
	os.system('rm -rf '+cal_table_prefix+'_2.fluxscale1')
	print 'Scaling amp gains for lower and upper base-band...'
	myscale = fluxscale(vis=ms_name,caltable=cal_table_prefix+'.G1', fluxtable=cal_table_prefix+'.fluxscale1',\
		reference=bpf_lst,transfer=second_lst+polleak_lst,incremental=False)
	myscale2 = fluxscale(vis=ms_name,caltable=cal_table_prefix+'_2.G1', fluxtable=cal_table_prefix+'_2.fluxscale1',\
		reference=bpf_lst,transfer=second_lst+polleak_lst,incremental=False)
	if intera=='y':
		print 'Plotting solutions...'
		print 'Lower Base-band...'
		plotcal(caltable=cal_table_prefix+'.fluxscale1',xaxis='time',yaxis='amp',poln='R',iteration='spw')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'.fluxscale1',xaxis='time',yaxis='amp',poln='L',iteration='spw')
		raw_input('Please press enter when ready to continue.')
		print 'Upper Base-band...'
		plotcal(caltable=cal_table_prefix+'_2.fluxscale1',xaxis='time',yaxis='amp',poln='R',iteration='spw')
		raw_input('Please press enter when ready to continue.')
		plotcal(caltable=cal_table_prefix+'_2.fluxscale1',xaxis='time',yaxis='amp',poln='L',iteration='spw')
		raw_input('Please press enter when ready to continue.')
	gt_lst_low.append(cal_table_prefix+'.fluxscale1')
	gt_lst_high.append(cal_table_prefix+'_2.fluxscale1')
	########################################

	########################################
	##Applying the calibration
	########################################
	print 'Applying calibration...'
	print 'First the cals...'
	print 'BP cal(s)...'
	for i in range(0,len(bpf_lst)):
		applycal(vis=ms_name,field=bpf_lst[i],gaintable=gt_lst_low,\
			gainfield=gf_lst_low+[bpf_lst[i]],interp=gi_lst_low+['nearest'],\
			calwt=[False],parang=False,spw=spw_low)
		applycal(vis=ms_name,field=bpf_lst[i],gaintable=gt_lst_high,\
			gainfield=gf_lst_high+[bpf_lst[i]],interp=gi_lst_high+['nearest'],\
			calwt=[False],parang=False,spw=spw_high)
	print 'Second cal(s)...'
	for i in range(0,len(second_lst)):
		applycal(vis=ms_name,field=second_lst[i],gaintable=gt_lst_low,\
			gainfield=gf_lst_low+[second_lst[i]],interp=gi_lst_low+['nearest'],\
			calwt=[False],parang=False,spw=spw_low)
		applycal(vis=ms_name,field=second_lst[i],gaintable=gt_lst_high,\
			gainfield=gf_lst_high+[second_lst[i]],interp=gi_lst_high+['nearest'],\
			calwt=[False],parang=False,spw=spw_high)
	if len(polleak_lst)>0:
		print 'Pol cal(s)...'
		for i in range(0,len(polleak_lst)):
			applycal(vis=ms_name,field=polleak_lst[i],gaintable=gt_lst_low,\
				gainfield=gf_lst_low+[polleak_lst[i]],interp=gi_lst_low+['nearest'],calwt=[False],\
				parang=False,spw=spw_low)
			applycal(vis=ms_name,field=polleak_lst[i],gaintable=gt_lst_high,\
				gainfield=gf_lst_high+[polleak_lst[i]],interp=gi_lst_high+['nearest'],calwt=[False],\
				parang=False,spw=spw_high)
	print 'Now the target...'
	for i in range(0,len(target_lst)):
		applycal(vis=ms_name,field=target_lst[i],\
			gaintable=gt_lst_low,gainfield=gf_lst_low+[second_cal],interp=gi_lst_low+['linear'],\
			calwt=[False],parang=False,spw=spw_low)
		applycal(vis=ms_name,field=target_lst[i],\
			gaintable=gt_lst_high,gainfield=gf_lst_high+[second_cal],interp=gi_lst_high+['linear'],\
			calwt=[False],parang=False,spw=spw_high)
	########################################

	########################################
	#check calibrated data
	##########################################
	if intera=='y':
		print 'Check the calibrated data...'
		print 'Lower base-band...'
		print '(1) BP spectra amp'
		plotms(vis=ms_name,field=bpf_cal,correlation='',timerange=timerbp,antenna='',avgtime='60s',\
			xaxis='channel',yaxis='amp',ydatacolumn='corrected',spw=spw_low)
		raw_input('Please press enter when ready to continue.')
		print '(2) Second cal amp spectra'
		for i in range(0,len(second_lst)):
			plotms(vis=ms_name,field=second_lst[i],correlation='RR,LL',timerange='',antenna='',avgtime='60s',\
				xaxis='channel',yaxis='amp',ydatacolumn='corrected',spw=spw_low)
			raw_input('Please press enter when ready to continue.')
		print '(3) Second cal phase spectra'
		for i in range(0,len(second_lst)):
			plotms(vis=ms_name,field=second_lst[i],correlation='RR,LL',timerange='',antenna='',avgtime='60s',\
				xaxis='channel',yaxis='phase',ydatacolumn='corrected',plotrange=[-1,-1,-180,180],\
				coloraxis='corr',spw=spw_low)
			raw_input('Please press enter when ready to continue.')
		print '(4) Second cal amp vs phase'
		for i in range(0,len(second_lst)):
			plotms(vis=ms_name,field=second_lst[i],correlation='RR,LL',timerange='',antenna='',avgtime='60s',\
				xaxis='phase',xdatacolumn='corrected',yaxis='amp',ydatacolumn='corrected',\
				plotrange=[-180,180,0,3],coloraxis='field',spw=spw_low)
			raw_input('Please press enter when ready to continue.')
		print '(5) Target amp vs uvdist'
		for i in range(0,len(target_lst)):
			plotms(vis=ms_name,field=target_lst[i],correlation='',timerange='',antenna='',avgtime='60s',\
				xaxis='uvdist',yaxis='amp',ydatacolumn='corrected',spw=spw_low)
			raw_input('Please press enter when ready to continue.')
		print '(6) Target amp vs time'
		for i in range(0,len(target_lst)):
			plotms(vis=ms_name,field=target_lst[i],correlation='',timerange='',antenna='',avgtime='60s',\
				xaxis='time',yaxis='amp',ydatacolumn='corrected',spw=spw_low)
			raw_input('Please press enter when ready to continue.')
		print 'Upper base-band...'
		print '(1) BP spectra amp'
		plotms(vis=ms_name,field=bpf_cal,correlation='',timerange=timerbp,antenna='',avgtime='60s',\
			xaxis='channel',yaxis='amp',ydatacolumn='corrected',spw=spw_high)
		raw_input('Please press enter when ready to continue.')
		print '(2) Second cal amp spectra'
		for i in range(0,len(second_lst)):
			plotms(vis=ms_name,field=second_lst[i],correlation='RR,LL',timerange='',antenna='',avgtime='60s',\
				xaxis='channel',yaxis='amp',ydatacolumn='corrected',spw=spw_high)
			raw_input('Please press enter when ready to continue.')
		print '(3) Second cal phase spectra'
		for i in range(0,len(second_lst)):
			plotms(vis=ms_name,field=second_lst[i],correlation='RR,LL',timerange='',antenna='',avgtime='60s',\
				xaxis='channel',yaxis='phase',ydatacolumn='corrected',plotrange=[-1,-1,-180,180],\
				coloraxis='corr',spw=spw_high)
			raw_input('Please press enter when ready to continue.')
		print '(4) Second cal amp vs phase'
		for i in range(0,len(second_lst)):
			plotms(vis=ms_name,field=second_lst[i],correlation='RR,LL',timerange='',antenna='',avgtime='60s',\
				xaxis='phase',xdatacolumn='corrected',yaxis='amp',ydatacolumn='corrected',\
				plotrange=[-180,180,0,3],coloraxis='field',spw=spw_high)
			raw_input('Please press enter when ready to continue.')
		print '(5) Target amp vs uvdist'
		for i in range(0,len(target_lst)):
			plotms(vis=ms_name,field=target_lst[i],correlation='',timerange='',antenna='',avgtime='60s',\
				xaxis='uvdist',yaxis='amp',ydatacolumn='corrected',spw=spw_high)
			raw_input('Please press enter when ready to continue.')
		print '(6) Target amp vs time'
		for i in range(0,len(target_lst)):
			plotms(vis=ms_name,field=target_lst[i],correlation='',timerange='',antenna='',avgtime='60s',\
				xaxis='time',yaxis='amp',ydatacolumn='corrected',spw=spw_high)
			raw_input('Please press enter when ready to continue.')
		extraf=raw_input('Do you need to do additional flagging? y or n-->')
		dict_log.append((ms_name_prefix+'_check_flag',extraf))
		while extraf=='y':
			badasfextra=raw_input('Please enter bad ant,spw,and field to flag (enter if none). e.g., ea10,ea12;5:4~9;3 ;5;3-->').split(' ')
			dict_log.append((ms_name_prefix+'_check_flag_antspwfield',badasfextra))
			if badasfextra=='':
				print 'Nothing to flag.'
				extraf=raw_input('Do you need to do additional flagging? y or n-->')
			else:
				print 'Flagging selected data.'
				for i in range(0,len(badasfextra)):
					strge=badasfextra[i].split(';')
					flagdata(vis=ms_name,flagbackup=T, mode='manual', antenna=strge[0],spw=strge[1],field=strge[2])
				extraf=raw_input('Do you need to do additional flagging? y or n-->')
		else:
			print 'No extra flagging requested.'
	##########################################

	##########################################
	##split out target
	##########################################
	for iii in range(0,len(target_lst)):
		os.system('rm -rf '+split_low.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms')
		os.system('rm -rf '+split_high.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms')
		os.system('rm -rf '+split_full.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms')
		print 'Splitting out target data for both base-bands...'
		split(vis=ms_name,outputvis=split_low.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms',\
			datacolumn='corrected',field=target_lst[iii],antenna='',spw=spw_low)
		split(vis=ms_name,outputvis=split_high.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms',\
			datacolumn='corrected',field=target_lst[iii],antenna='',spw=spw_high)
		split(vis=ms_name,outputvis=split_full.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms',\
			datacolumn='corrected',field=target_lst[iii],antenna='',spw=spw_full)
		split_low=split_low.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms'
		split_high=split_high.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms'
		split_full=split_full.strip('.ms')+'_TID'+str(target_lst[iii])+'.ms'
		##########################################
		if doImage=='T':
			###########################################
			##Imaging
			###########################################
			if intera=='n':
				print '**********************************************'
				print 'Entering interactive part of script again...'
				print '**********************************************'
			print 'Imaging...'
			if use_auto=='T':
				myimsize=set_imagesize(split_low,0,'0')
				mycell=set_cellsize(split_low,0)
			print 'Lower base-band...'
			if mymask=='':
				os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_low+'_clean1*')
				print 'Using interactive mode so you can make a mask...'
				print 'Cleaning...'
				clean(vis=split_low, imagename=my_dir+target+'_'+date+'_'+band_low+'_clean1',field='',spw='',interactive=T,\
					cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
					mode='mfs',niter=0,nterms=mynterms,stokes=mystokes,multiscale=multiscale,robust=robust,outlierfile=outlierfile)
			else:
				os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_low+'_clean1*')
				print 'Cleaning...'
				clean(vis=split_low, imagename=my_dir+target+'_'+date+'_'+band_low+'_clean1',field='',mask=mymask,spw='',interactive=F,\
					cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
					mode='mfs',niter=myniter,nterms=mynterms,stokes=mystokes,multiscale=multiscale,robust=robust,outlierfile=outlierfile)
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
			fluxl,errl,unitl,freql,errl_real=imfit_point(imagenl,my_dir)
			print 'Lower base-band flux density of ',fluxl,' +/- ',errl, unitl
			print 'Local RMS in Lower base-band image is: ',errl_real,' Jy'
			dopscl=raw_input('Do you want to do phase selfcal?y or n-->')
			dict_log.append((ms_name_prefix+'_phself_lsb',dopscl))
			if dopscl=='y':
				selfcal_low,scim_low=phselfcal(split_low,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
			date,band_low,'n',spw_low,outlierfile,multiscale,robust,weighting)
				fluxl_sc,errl_sc,unitl_sc,freql_sc,errl_real_sc=imfit_point(scim_low,my_dir)
			
			print 'Upper base-band...'
			if use_auto=='T':
				myimsize=set_imagesize(split_high,0,'0')
				mycell=set_cellsize(split_high,0)
			if mymask=='':
				os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_high+'_clean1*')
				print 'Using interactive mode so you can make a mask...'
				print 'Cleaning...'
				clean(vis=split_high, imagename=my_dir+target+'_'+date+'_'+band_high+'_clean1',field='',spw='',interactive=T,\
					cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
					mode='mfs',niter=0,nterms=mynterms,stokes=mystokes,multiscale=multiscale,robust=robustoutlierfile=outlierfile)
			else:
				os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_high+'_clean1*')
				print 'Cleaning...'
				clean(vis=split_high, imagename=my_dir+target+'_'+date+'_'+band_high+'_clean1',field='',mask=mymask,spw='',interactive=F,\
					cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
					mode='mfs',niter=myniter,nterms=mynterms,stokes=mystokes,multiscale=multiscale,robust=robustoutlierfile=outlierfile)
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
			fluxu,erru,unitu,frequ,erru_real=imfit_point(imagenu,my_dir)
			print 'Upper base-band flux density of ',fluxu,' +/- ',erru, unitu
			print 'Local RMS in Upper base-band image is: ',erru_real,' Jy'
			dopscu=raw_input('Do you want to do phase selfcal?y or n-->')
			dict_log.append((ms_name_prefix+'_phself_usb',dopscu))
			if dopscu=='y':
				selfcal_high,scim_high=phselfcal(split_high,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
			date,band_low,'n',spw_high,outlierfile,multiscale,robust,weighting)
				fluxu_sc,erru_sc,unitu_sc,frequ_sc,erru_real_sc=imfit_point(scim_high,my_dir)

			print 'Combined base-band...'
			if use_auto=='T':
				myimsize=set_imagesize(split_full,0,'0')
				mycell=set_cellsize(split_full,0)
			if mymask=='':
				os.system('rm -rf '+my_dir+target+'_'+date+'_both_clean1*')
				print 'Using interactive mode so you can make a mask...'
				print 'Cleaning...'
				clean(vis=[split_low,split_high], imagename=my_dir+target+'_'+date+'_both_clean1',field='',spw='',interactive=T,\
					cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
					mode='mfs',niter=0,nterms=mynterms,stokes=mystokes,multiscale=multiscale,robust=robust,outlierfile=outlierfile)
			else:
				os.system('rm -rf '+my_dir+target+'_'+date+'_both_clean1*')
				print 'Cleaning...'
				clean(vis=[split_low,split_high], imagename=my_dir+target+'_'+date+'_both_clean1',field='',mask=mymask,spw='',interactive=F,\
					cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
					mode='mfs',niter=myniter,nterms=mynterms,stokes=mystokes,multiscale=multiscale,robust=robust,outlierfile=outlierfile)
			if mynterms>1:
				imagenb=my_dir+target+'_'+date+'_both_clean1.image.tt0'
			else:
				imagenb=my_dir+target+'_'+date+'_both_clean1.image'
			print 'Correcting for PB...'
			os.system('rm -rf '+imagenb+'.pbcor')
			os.system('rm -rf '+imagenb+'.pbcor.fits')
			immath(imagename=[imagenb,my_dir+target+'_'+date+'_both_clean1.flux'],\
				expr='IM0/IM1',outfile=imagenb+'.pbcor')
			print 'Making fits image...'
			exportfits(imagename=imagenb+'.pbcor',fitsimage=imagenb+'.pbcor.fits')
			imagenb=imagenb+'.pbcor'
			fluxb,errb,unitb,freqb,errb_real=imfit_point(imagenb,my_dir)
			print 'Combined base-band flux density of ',fluxb,' +/- ',errb, unitb
			print 'Local RMS in Combined base-band image is: ',errb_real,' Jy'
			dopscb=raw_input('Do you want to do phase selfcal?y or n-->')
			dict_log.append((ms_name_prefix+'_phself_both',dopscb))
			if dopscb=='y':
				selfcal_both,scim_both=phselfcal(split_full,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
			date,band_low,'n',spw_full,outlierfile,multiscale,robust,weighting)
				fluxb_sc,errb_sc,unitb_sc,freqb_sc,errb_real_sc=imfit_point(scim_both,my_dir)

			#writing imfit result to file
			print 'Writing imfit results to file...'
			if firstf=='y':
				resul_file=open(my_dir+'imfit_results.txt','w')
			else:
				resul_file=open(my_dir+'imfit_results.txt','a')
			if dopscl=='n' and dopscu=='n' and dopscb=='n':
				resul_file.write('Target'+str(target_id[iii])+':\n')
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freql,fluxl,errl,unitl,errl_real))
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,frequ,fluxu,erru,unitu,erru_real))
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freqb,fluxb,errb,unitb,errb_real))
			else:
				resul_file.write('Target'+str(target_id[iii])+':\n')
				resul_file.write('Pre-selfcal:\n')
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freql,fluxl,errl,unitl,errl_real))
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,frequ,fluxu,erru,unitu,erru_real))
				resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freqb,fluxb,errb,unitb,errb_real))
				resul_file.write('Post selfcal:\n')
				if dopscl=='y':
					resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freql_sc,fluxl_sc,errl,unitl_sc,errl_real_sc))
				if dopscu=='y':
					resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,frequ_sc,fluxu_sc,erru_sc,unitu_sc,erru_real_sc))
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
				stokes_param=mystokes
				print 'Lower base-band...'
				comblsb=split_low.strip('.ms')
				mstransform(vis=split_low, outputvis=comblsb+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
				fitfulluv_low=uvm.uvmultifit(vis=comblsb+'_mstrans.ms', spw='', column = "data", \
					uniform=False, model=[comp_uv],stokes = stokes_param,outfile=my_dir+'lsbmodelfit.dat',\
					var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
				src_uv_init_low=fitfulluv_low.result['Parameters'][2]
				src_uv_err_low=fitfulluv_low.result['Uncertainties'][2]
				print 'Upper base-band...'
				combusb=split_high.strip('.ms')
				mstransform(vis=split_high, outputvis=combusb+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
				fitfulluv_high=uvm.uvmultifit(vis=combusb+'_mstrans.ms', spw='', column = "data", \
					uniform=False, model=[comp_uv],stokes = stokes_param,outfile=my_dir+'usbmodelfit.dat',\
					var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
				src_uv_init_high=fitfulluv_high.result['Parameters'][2]
				src_uv_err_high=fitfulluv_high.result['Uncertainties'][2]
				print 'Combined base-band...'
				combboth=split_full.strip('.ms')
				mstransform(vis=split_full, outputvis=combboth+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
				fitfulluv=uvm.uvmultifit(vis=combboth+'_mstrans.ms', spw='', column = "data", \
					uniform=False, model=[comp_uv],stokes = stokes_param,outfile=my_dir+'combmodelfit.dat',\
					var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
				src_uv_init=fitfulluv.result['Parameters'][2]
				src_uv_err=fitfulluv.result['Uncertainties'][2]
				print 'Writing uvfit results to file...'
				if firstf=='y':
					resuluv_file=open(my_dir+'uvfit_results.txt','w')
				else:
					resuluv_file=open(my_dir+'uvfit_results.txt','a')
				resuluv_file.write('{0} {1} {2} {3}\n'.format(band,'lower',src_uv_init_low,src_uv_err_low))
				resuluv_file.write('{0} {1} {2} {3}\n'.format(band,'upper',src_uv_init_high,src_uv_err_high))
				resuluv_file.write('{0} {1} {2} {3}\n'.format(band,'combined',src_uv_init,src_uv_err))
				resuluv_file.close()

	if ms_name_list[kk]==ms_name_list[-1]:
		print 'Finished ', ms_name,'.'
		if doImage=='T':
			print 'All split data sets have been reduced and imaged.'
		else:
			print 'All split data sets have been reduced. No imaging was performed.'
	else:
		print 'Finished ', ms_name,'. Moving on to next split data set.'
	###########################################

print 'Cleaning up...'
os.system('rm -rf *.log')
os.system('rm -rf *.last')
print 'Writing user_input log file...'
writeDict(dict_log, my_dir+'user_input.log',str(datetime.datetime.now()))
print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'
#need to add pol cal!!!
