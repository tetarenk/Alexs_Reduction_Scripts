import numpy as np
import os
import analysisUtils as au
import casac
from tasks import *
from taskinit import *
from collections import OrderedDict

def phselfcal(visi,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
	date,bband,combi,spw,outlierf):
	cal_table_prefix=my_dir+target+'_'+date+'_'+bband
	cont0='y'
	attemptnum=1
	scname=visi.strip('.ms')
	while cont0=='y':
		cont='y'
		#make copy of data set for selfcal
		print 'Making copy of data set for selfcal...'
		os.system('rm -rf '+scname+'_selfcal_'+str(attemptnum)+'.ms')
		split(vis=visi,outputvis=scname+'_selfcal_'+str(attemptnum)+'.ms',
	      datacolumn='data',field='',spw='')
		selfcalvis=scname+'_selfcal_'+str(attemptnum)+'.ms'
		#light clean to get source model
		print 'Interactive Clean to get source model...'
		os.system('rm -rf '+my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean0*')
		clean(vis=selfcalvis, imagename=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean0',field='',spw=spw,interactive=T,\
			cell=[mycell], imsize=myimsize,gain=0.1,weighting='natural',threshold=mythreshold,mode='mfs',niter=0,nterms=mynterms,outlierfile=outlierf)
		raw_input('Please press enter when ready to continue.')
		#solve for phase gains
		while cont=='y':
			solu=raw_input('Please enter solution interval in sec. e.g., 10-->')
			if os.path.isdir(cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself'):
				print 'This selfcal table already exists.'
			else:
				print 'Solving for phase selfcal solutions at ', solu,' sec interval...'
				if combi=='y':
					gaincal(vis=selfcalvis, caltable=cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself',\
						field='',selectdata=False,solint=solu+'s',refant=ref_ant, minblperant=4,\
						minsnr=3, solnorm=True,append=False, gaintype='G',calmode='p',combine='spw')
				else:
					gaincal(vis=selfcalvis, caltable=cal_table_prefix+'_'+solu+'sec'+str(attemptnum)+'.phself',\
						field='',selectdata=False,solint=solu+'s',refant=ref_ant, minblperant=4,\
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
		print 'Interactive Cleaning selfcaled data...'
		os.system('rm -rf '+my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1*')
		clean(vis=selfcalvis, imagename=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1',field='',spw=spw,interactive=T,\
			cell=[mycell], imsize=myimsize,gain=0.1,weighting='natural',threshold=mythreshold,mode='mfs',niter=0,nterms=mynterms,outlierfile=outlierf)
		raw_input('Please press enter when ready to continue.')
		print 'Viewing selfcaled image...'
		if mynterms>1:
			scim=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1.image.tt0'
		else:
			scim=my_dir+target+'_'+date+'_'+bband+'_phasesc'+str(attemptnum)+'_clean1.image'
		imview(scim)
		raw_input('Please press enter when ready to continue.')
		cont0=raw_input('Do you want to redo selfcal?y or n-->')
		if cont0=='y':
			attemptnum=attemptnum+1
	return(selfcalvis,scim)

def imfit_point(pbimage,my_dir):
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
	fit_res=imfit(imagename=pbimage, box=box, estimates=tempFile.name, append=False, overwrite = True)
	result_box1=imstat(imagename=pbimage,region='annulus['+cen_annulus+','+cen_radius+']')
	os.system('rm -rf '+my_dir+'tempfile_fit.txt')
	return(fit_res['results']['component0']['flux']['value'][0],fit_res['results']['component0']['flux']['error'][0],\
		fit_res['results']['component0']['flux']['unit'],fit_res['results']['component0']['spectrum']['frequency']['m0']['value'],result_box1['rms'][0])

def writeDict(dicti, filename,logdate):
	ordd=OrderedDict(dicti)
	with open(filename, "a") as f:
		f.write(logdate+':'+'\n')
		for i in ordd.keys():            
			f.write(i + " : " + str(ordd[i]) + "\n")
	


