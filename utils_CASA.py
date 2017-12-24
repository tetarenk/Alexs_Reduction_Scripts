import numpy as np
import os
import analysisUtils as au
import casac
from tasks import *
from taskinit import *
from collections import OrderedDict
import linecache
from datetime import datetime
import glob
from astropy.io import ascii

def phselfcal(visi='',mycell='',mynterms='',myimsize='',mythreshold='',ref_ant='',my_dir='',target='',\
	date='',bband='',combi='',outlierf='',multiscale='',robust='',weighting='',help='F'):
	'''Selfcal: can use gaintype='T' to combine polarizations and combine='spw' or 'scan' or both if low S/N '''
	if help=='T':
		print 'arguments for selfcal function are: msname,cellsize,nterms,imsize,threashold,ref_ant,my_dir,target,date,bband,combi,outlierf,multiscale,robust,weighting'
		return
	cal_table_prefix=my_dir+target+'_'+date+'_'+bband
	cont0='y'
	attemptnum=1
	scname=visi.strip('.ms')
	os.system('rm -rf '+visi+'.flagversions')
	os.system('rm -rf '+cal_table_prefix+'*.phself')
	os.system('rm -rf '+cal_table_prefix+'*.ampself')
	while cont0=='y':
		cont='y'
		if attemptnum==1:
			#make copy of data set for selfcal
			print 'Making copy of data set for selfcal...'
			os.system('rm -rf '+scname+'_selfcal_'+'presc'+'.ms')
			os.system('rm -rf '+scname+'_selfcal_'+'presc'+'.ms.flagversions')
			os.system('rm -rf '+scname+'_selfcal_'+'postsc'+'.ms')
			os.system('rm -rf '+scname+'_selfcal_'+'postsc'+'.ms.flagversions')
			split(vis=visi,outputvis=scname+'_selfcal_'+'presc'+'.ms',
		      datacolumn='data',field='',spw='')
			selfcalvis=scname+'_selfcal_'+'presc'+'.ms'
			#save flags
			flagmanager(vis=selfcalvis,mode='save',versionname='before_selfcal',merge='replace')
			#light clean to get source model
			print 'Interactive Clean to get source model...'
			os.system('rm -rf '+my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean0*')
			clean(vis=selfcalvis, imagename=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean0',field='',spw='',interactive=True,\
				cell=mycell, imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,mode='mfs',niter=0,nterms=1,outlierfile=outlierf,multiscale=multiscale,robust=robust)
			raw_input('Please press enter when ready to continue.')
		else:
			selfcalvis=scname+'_selfcal_'+'presc'+'.ms'
			raw_input('Please press enter when ready to continue.')
		#solve for phase gains
		while cont=='y':
			solu=raw_input('Please enter solution interval in sec. e.g., 10-->')
			if solu in ['inf','int']:
				soluu=solu
			else:
				soluu=solu+'s'
			if os.path.isdir(cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself'):
				print 'This selfcal table already exists.'
			else:
				print 'Solving for phase selfcal solutions at ', solu,' sec interval...'
				if combi=='y':
					gaincal(vis=selfcalvis, caltable=cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself',\
						field='',selectdata=False,solint=soluu,refant=ref_ant, minblperant=4,\
						minsnr=3, solnorm=True,append=False, gaintype='G',calmode='p',combine='spw')
				else:
					gaincal(vis=selfcalvis, caltable=cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself',\
						field='',selectdata=False,solint=soluu,refant=ref_ant, minblperant=4,\
						minsnr=3, solnorm=True,append=False, gaintype='G',calmode='p')
				#plot solutions
				print 'Plotting solutions. Check that they smoothly vary with time.'
				plotcal(caltable=cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself',xaxis='time',yaxis='phase',field='',subplot=331,\
					overplot=False,clearpanel='Auto',iteration='antennas',plotsymbol='o',plotcolor='blue',\
					showgui=True,figfile='')
				raw_input('Please press enter when ready to continue.')
			cont=raw_input('Do you wish to try another solution interval?y or n-->')
			if cont =='y':
				print 'Redoing phase solutions...'
				
		scchoice=raw_input('What solution interval do you want to apply? e.g., 10-->')
		sctable=cal_table_prefix+'_'+str(scchoice)+'sec'+str(attemptnum)+'.phself'
		print 'Applying phase selfcal solutions at ', scchoice,' sec interval...'
		if combi=='y':
			applycal(vis=selfcalvis, field='',spw='',selectdata=False, gaintable= [sctable],gainfield=[''],\
				interp=['nearest'],spwmap=[[0]], calwt=[False])
		else:
			applycal(vis=selfcalvis, field='',spw='',selectdata=False, gaintable= [sctable],gainfield=[''],\
				interp=['nearest'], calwt=[False])
		flagmanager(vis=selfcalvis,mode='save',versionname='after_phase'+str(attemptnum))
		print 'Interactive Cleaning selfcaled data...'
		os.system('rm -rf '+my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1*')
		clean(vis=selfcalvis, imagename=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1',field='',spw='',interactive=True,\
			cell=mycell, imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,mode='mfs',niter=0,nterms=1,outlierfile=outlierf,multiscale=multiscale,robust=robust)
		raw_input('Please press enter when ready to continue.')
		print 'Viewing selfcaled image...'
		scim=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1.image'
		imview(scim)
		raw_input('Please press enter when ready to continue.')
		cont0=raw_input('Do you want to continue phase selfcal?y or n-->')
		if cont0=='y':
			attemptnum=attemptnum+1
	#amp selfcal
	ampsc=raw_input('Do you want to try amp selfcal? (y or n)-->')
	amp_sc_int='y'
	if ampsc=='y':
		while amp_sc_int=='y':
			soluamp=raw_input('Please enter solution interval in sec. e.g., inf-->')
			if soluamp in ['inf','int']:
				soluampp=soluamp
			else:
				soluampp=soluamp+'s'
			if os.path.isdir(cal_table_prefix+'_'+soluamp+'sec'+'.ampself'):
				print 'This selfcal table already exists.'
			else:
				print 'Solving for amp selfcal solutions at ', soluamp,' sec interval...'
				if combi=='y':
					gaincal(vis=selfcalvis, caltable=cal_table_prefix+'_'+soluamp+'sec'+'.ampself',\
						field='',selectdata=False,solint=soluampp,refant=ref_ant, minblperant=4,minsnr=3,\
						solnorm=True,append=False, gaintype='G',calmode='ap',combine='spw')
				else:
					gaincal(vis=selfcalvis, caltable=cal_table_prefix+'_'+soluamp+'sec'+'.ampself',\
						field='',selectdata=False,solint=soluampp,refant=ref_ant, minblperant=4,minsnr=3,\
						solnorm=True,append=False, gaintype='G',calmode='ap')
				print 'Plotting solutions. Check that they smoothly vary with time.'
				plotcal(caltable=cal_table_prefix+'_'+soluamp+'sec'+'.ampself',xaxis='time',\
					yaxis='amp',field='',subplot=331,overplot=False,clearpanel='Auto',\
					iteration='antennas',plotsymbol='o',plotcolor='blue',showgui=True,figfile='')
				raw_input('Please press enter when ready to continue.')
			amp_sc_int=raw_input('Do you wish to try another solution interval?y or n-->')
			if amp_sc_int=='y':
				print 'Redoing amp solutions...'
		scchoiceamp=raw_input('What solution interval do you want to apply? e.g., inf-->')
		sctablea=cal_table_prefix+'_'+str(scchoiceamp)+'sec'+'.ampself'
		print 'Applying phase selfcal solutions at ', scchoiceamp,' sec interval...'
		if combi=='y':
			applycal(vis=selfcalvis, field='',spw='',selectdata=False, gaintable= [sctablea],\
				gainfield=[''],interp=['nearest'],spwmap=[[0]], calwt=[False])
		else:
			applycal(vis=selfcalvis, field='',spw='',selectdata=False, gaintable= [sctablea],\
				gainfield=[''],interp=['nearest'], calwt=[False])
		flagmanager(vis=selfcalvis,mode='save',versionname='after_apcal')
		print 'Interactive Cleaning selfcaled data...'
		clean(vis=selfcalvis, imagename=my_dir+target+'_'+date+'_'+bband+'_phaseampsc'+'_clean1',\
			field='',spw='',interactive=True,cell=mycell, imsize=myimsize,gain=0.1,weighting=weighting,\
			threshold=mythreshold,mode='mfs',niter=0,nterms=1,outlierfile=outlierf,\
			multiscale=multiscale,robust=robust)
		raw_input('Please press enter when ready to continue.')
		print 'Viewing selfcaled image...'
		scim=my_dir+target+'_'+date+'_'+bband+'_phaseampsc'+'_clean1.image'
		imview(scim)
		raw_input('Please press enter when ready to continue.')
		del_amp=raw_input('Do you want to revert to pre amp selfcal? y or n-->')
		if del_amp=='y':
			flagmanager(vis=selfcalvis, mode='restore',versionname='after_phase'+str(attemptnum))
			clearcal(selfcalvis)
			delmod(selfcalvis,field='',otf=True)
			if combi=='y':
				applycal(vis=selfcalvis, field='',spw='',selectdata=False, gaintable= [sctable],\
					gainfield=[''],interp=['nearest'],spwmap=[[0]], calwt=[False])
			else:
				applycal(vis=selfcalvis, field='',spw='',selectdata=False, gaintable= [sctable],\
					gainfield=[''],interp=['nearest'], calwt=[False])
	print 'Splitting out final selfcaled data set...'
	os.system('rm -rf '+scname+'_selfcal_'+'postsc'+'.ms')
	split(vis=selfcalvis,outputvis=scname+'_selfcal_'+'postsc'+'.ms',field='',spw='')
	selfcalvisfin=scname+'_selfcal_'+'postsc'+'.ms'
	return(selfcalvisfin,scim)

def imfit_point(pbimage='',my_dir='',stokes='',help='F'):
	if help=='T':
		print 'arguments for imfit function are: imagename,my_dir,stokes'
		return
	if stokes=='I':
		ind_st=0
	elif stokes=='Q':
		ind_st=1
	elif stokes=='U':
		ind_st=2
	elif stokes=='V':
		ind_st=3
	elif stokes=='':
		ind_st=0
	else:
		raise Exception('Please enter valid stokes param; I, Q, U, or V.')
	print 'Viewing image...'
	imview(pbimage)
	raw_input('Please press enter when ready to continue.')
	box=raw_input('Please enter target box. e.g., x1,y1,x2,y2-->')
	annulus_rad_inner=raw_input('Please enter inner annulus radius for RMS estimation (e.g., 20)-->')
	annulus_rad_outer=raw_input('Please enter outer annulus radius for RMS estimation (e.g., 30)-->')
	x_sizel=float(box.split(',')[0])
	x_sizeu=float(box.split(',')[2])
	y_sizel=float(box.split(',')[1])
	y_sizeu=float(box.split(',')[3])
	cen_annulus='['+str(((x_sizeu-x_sizel)/2.)+x_sizel)+'pix,'+str(((y_sizeu-y_sizel)/2.)+y_sizel)+'pix]'
	cen_radius='['+str(annulus_rad_inner)+'pix,'+str(annulus_rad_outer)+'pix]'
	print 'Fitting point source...'
	beamMajor = imhead(imagename=pbimage,mode='get',hdkey='beammajor')
	beamMajor = str(beamMajor['value'])+'arcsec'
	beamMinor = imhead(imagename=pbimage,mode='get',hdkey='beamminor')
	beamMinor = str(beamMinor['value'])+'arcsec'
	beamPA = imhead(imagename=pbimage,mode='get',hdkey='beampa')
	beamPA = str(beamPA['value'])+'deg'
	imstatOut = imstat(imagename=pbimage,box=box)
	peak = imstatOut['max']
	peak = str(peak[0])
	peakPosValue = imstatOut['maxpos']
	peakPosX = str(peakPosValue[0])
	peakPosY = str(peakPosValue[1])
	tempFile = open(my_dir+'tempfile_fit.txt','w')
	mystring = str(peak+', '+peakPosX+', '+peakPosY+', '+beamMajor+', '+beamMinor+', '+beamPA)
	tempFile.write(mystring)
	tempFile.close()
	result_box1=imstat(imagename=pbimage,region='annulus['+cen_annulus+','+cen_radius+']')
	if stokes=='I' or stokes=='':
		fit_res=imfit(imagename=pbimage, box=box, estimates=tempFile.name, append=False, overwrite = True)
		os.system('rm -rf '+my_dir+'tempfile_fit.txt')
		if fit_res['results']['nelements']>0:
			return(fit_res['results']['component0']['flux']['value'][ind_st],fit_res['results']['component0']['flux']['error'][ind_st],\
				fit_res['results']['component0']['flux']['unit'],fit_res['results']['component0']['spectrum']['frequency']['m0']['value'],result_box1['rms'][0])
		else:
			return(0,0,0,0,0)
	else:
		return(0,0,0,0,result_box1['rms'][0])

def writeDict(dicti, filename,logdate):
	ordd=OrderedDict(dicti)
	with open(filename, "a") as f:
		f.write(logdate+':'+'\n')
		for i in ordd.keys():            
			f.write(i + " : " + str(ordd[i]) + "\n")


def setjy_parse(spw):
	'''After setjy call, parses latest casa log to find values needed for polarization calibration
	INPUT: None
	OUTPUT: lower flux, upper flux, lower freq, upper freq'''
	spwv=spw.split('~')[0]
	newest = max(glob.iglob('casa*.log'), key=os.path.getctime)
	#newest=files='/Users/atetarenk/Desktop/test_sjparse.txt'
	file_setjy=open(newest,'r')
	readit=file_setjy.read()
	lines=readit.splitlines()
	line1=[line for line in lines if "Scaling spw(s) ["+spwv in line]
	#2015-08-14 16:35:17 INFO imager Scaling spw(s) [8, 9, 10, 11, 12, 13, 14, 15]'s model image by channel to  I = 5.7914, 5.51751, 5.2727 Jy @(6.93716e+09, 7.44916e+09, 7.95916e+09)Hz for visibility prediction (a few representative values are shown).
	ind1=line1[0].index('I =')
	ind2=line1[0].index('Hz')
	line2a=line1[0]
	line2=line2a[ind1:ind2+2].split('@')
	fluxes=line2[0].strip('I = ').strip(' Jy ').split(',')
	freqs=line2[1].strip('(').strip(')Hz').split(',')
	return(float(fluxes[0]),float(fluxes[2]),float(freqs[0])/1e6,float(freqs[2])/1e6)

def last_field_parse(listobs_file):
	'''Parses listibs for last field and first antenna
	INPUT: listobs file
	OUTPUT: strings of values'''
	
	file_listobs=open(listobs_file,'r')
	readit=file_listobs.read()
	lines=readit.splitlines()
	line1=[line for line in lines if "nRows =" in line]
	ind1=lines.index(line1[0])
	line2=lines[ind1-1]
	last_field=filter(None,line2.split(' '))[4]
	line1a=[line for line in lines if "Antennas:" in line]
	ind1a=lines.index(line1a[0])
	line2a=lines[ind1a+3]
	first_ant=filter(None,line2a.split(' '))[1]
	return(last_field,first_ant)

def source_dict_create(listobs_file):
	'''Parses a list obs file and creates a dictionary of source attributes 
	INPUT: listobs file
	OUTPUT: dictionary
	'''
	file_listobs=open(listobs_file,'r')
	readit=file_listobs.read()
	lines=readit.splitlines()
	fields_line=[line for line in lines if "Fields:" in line]
	scan_line=[line for line in lines if "ScanIntent" in line]
	num_fields=int(fields_line[0].split(':')[1])
	fields_section=lines[lines.index(fields_line[0])+2:lines.index(fields_line[0])+num_fields+2]
	scans_section=lines[lines.index(scan_line[0])+1:lines.index(fields_line[0])-1]

	newDict = {}
	newDict['Fields']={}
	for i in fields_section:
		entry=i.split(' ')
		entry2 = filter(None, entry)
		newDict['Fields'][entry2[0]]={}
		newDict['Fields'][entry2[0]]['Name']=entry2[2]
	for i in scans_section:
		entrys=i.split(' ')
		entry2s = filter(None, entrys)
		intents=entry2s[-1]
		integs=entry2s[-2].strip(']')
		scants=entry2s[0:3]
		i1=intents.strip('[').strip(']')
		i2=i1.split(',')
		if 'Intent' not in newDict['Fields'][entry2s[4]]:
			newDict['Fields'][entry2s[4]]['Intent']=i2
			if 'CALIBRATE_BANDPASS#UNSPECIFIED' in i2 or 'CALIBRATE_FLUX#UNSPECIFIED' in i2:
				val1=scants[0]
				if '/' in val1:
					val1=scants[0].split('/')[1]
				newDict['Fields'][entry2s[4]]['scantimes']=val1+'~'+scants[2]
			else:
				foo=2
		else:
			if 'CALIBRATE_BANDPASS#UNSPECIFIED' in i2 or 'CALIBRATE_FLUX#UNSPECIFIED' in i2:
				newDict['Fields'][entry2s[4]]['scantimes']=val1+'~'+scants[2]
			else:
				foo=2
	file_listobs.close()
	return(newDict,integs)



	



