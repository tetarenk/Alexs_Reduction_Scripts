import numpy as np
import os
import analysisUtils as au
import casac
from tasks import *
from taskinit import *
from collections import OrderedDict

def phselfcal(visi='',mycell='',mynterms='',myimsize='',mythreshold='',ref_ant='',my_dir='',target='',\
	date='',bband='',combi='',outlierf='',multiscale='',robust='',weighting='',help='F'):
	'''Note: (1) 'int' is integration time, 'inf' is scan time (start here)
		(2) Avoid Amp self cal if possible, only use 'inf' interval
		(3) If get lots of failed solutions, try:
		 (a) gaintype='T' to combine polarizations
		 (b) combine='spw' or 'scan' or both'''
	if help=='T':
		print 'arguments for selfcal function are: msname,cellsize,nterms,imsize,threashold,ref_ant,my_dir,target,date,bband,combi,outlierf,multiscale,robust,weighting'
		return
	cal_table_prefix=my_dir+target+'_'+date+'_'+bband
	cont0='y'
	attemptnum=1
	scname=visi.strip('.ms')
	os.system('rm -rf '+visi+'.flagversions')
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
				cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,mode='mfs',niter=0,nterms=mynterms,outlierfile=outlierf,multiscale=multiscale,robust=robust)
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
			cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,mode='mfs',niter=0,nterms=mynterms,outlierfile=outlierf,multiscale=multiscale,robust=robust)
		raw_input('Please press enter when ready to continue.')
		print 'Viewing selfcaled image...'
		if mynterms>1:
			scim=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1.image.tt0'
		else:
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
			field='',spw='',interactive=True,cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,\
			threshold=mythreshold,mode='mfs',niter=0,nterms=mynterms,outlierfile=outlierf,\
			multiscale=multiscale,robust=robust)
		raw_input('Please press enter when ready to continue.')
		print 'Viewing selfcaled image...'
		if mynterms>1:
			scim=my_dir+target+'_'+date+'_'+bband+'_phaseampsc'+'_clean1.image.tt0'
		else:
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
	if stokes=='I':
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



	



