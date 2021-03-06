#################################################
#CASA Continuum Imaging/Selfcal Script
#################################################
'''CASA script to be used for imaging and selfcal of VLA/SMA continuum data
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) continuum image  -- [target]_[obsDate]_[band]_[subband]_clean1.image(.tt0).pbcor
        (2) (optional) Stokes Cube -- [target]_[obsDate]_[band]_[subband]_polcube_IQUV.image(.tt0).pbcor
        (3) (optional) Individual Stokes Images -- [target]_[obsDate]_[band]_[subband]_polcube.[I,Q,U, or V]
        (4) (optional) Polarization PA and Fractional Polarization images -- [target]_[obsDate]_[band]_[subband]_polcube.[PA or FP]
        (5) File of flux densities from image/UV plane fitting -- [target]_[obsDate]_[subbband]_imfit(uvfit)_results.txt
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
       - All output images are also converted to fits format (just append .fits to end of images 1 above)
       - This script images and selfcals an already calibrated CASA MS and fits a point source in the image/uv plane.
Written by: Alex J. Tetarenko
Last Updated: Oct. 2018
Works in CASA-5.3 now!'''

print '##################################################'
print 'Welcome to Alexs CASA Continuum Imaging/Selfcal Script'
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
from utils_CASA import imfit_point, phselfcal
from ekoch_casa_tools import set_imagermode,has_field,set_cellsize,set_imagesize,find_expected_beams,getBaselinePercentile,get_mosaic_info
import uvmultifit as uvm
import analysisUtils as au

#define output directory
my_dir='/mnt/bigdata/tetarenk/VLA_1758_1740/'
if not os.path.isdir(my_dir):
	os.system('mkdir '+my_dir)
	os.system('chown ubuntu '+my_dir)
	os.system('chmod -R u+r '+my_dir) 
	os.system('chmod -R u+w '+my_dir)
	os.system('chmod -R u+x '+my_dir)
print 'You have set your output directory to ', my_dir
print 'All output images & intermediate data products are put in this directory.\n'

#param file location
param_dir_file='/mnt/bigdata/tetarenk/CASA_reduction_scripts/params_cont_img.txt'
print 'You have set your param file to ', param_dir_file
print 'Please make sure all parameters are correct, they will change for each data set!\n'

#################################################
#Defining Parameters Section
#################################################
print 'Reading in parameters...'
#read in param file
f=open(param_dir_file)
data_params=imp.load_source('data_params','',f)
f.close()
#data set params
ms_name=data_params.ms_name
obsDate=data_params.obsDate
target=data_params.target
spw=data_params.spw
band=data_params.band
subband=data_params.subband
ref_ant=data_params.ref_ant
do_pol=data_params.do_pol
#general image params
use_auto=data_params.use_auto
mythreshold=data_params.mythreshold
myimsize=data_params.myimsize
mycell=data_params.mycell
mynterms=data_params.mynterms
decon=data_params.decon
myniter=data_params.myniter
mystokes=data_params.mystokes
outlierf=data_params.outlierf
combi=data_params.combi
multiscale=data_params.multiscale
robust=data_params.robust
weighting=data_params.weighting
#mask options
mymask=data_params.mymask
#uv fitting
uv_fit=data_params.uv_fit
uv_initp=data_params.uv_initp
#################################################

#################################################
#Imaging Section
#################################################
if use_auto=='T':
	myimsize=set_imagesize(ms_name,0,'0')
	mycell=set_cellsize(ms_name,0)
print 'Imaging...'
if mymask=='':
	os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1*')
	print 'Using interactive mode so you can make a mask...'
	print 'Cleaning...'
	tclean(vis=ms_name, imagename=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1',field='',spw=spw,interactive=True,\
		cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,deconvolver=decon,gridder='standard',\
		specmode='mfs',niter=myniter,nterms=mynterms,stokes=mystokes,outlierfile=outlierf,scales=multiscale,robust=robust)
else:
	os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1*')
	print 'Cleaning...'
	tclean(vis=ms_name, imagename=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1',field='',mask=mymask,spw=spw,interactive=False,\
		cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,deconvolver=decon,gridder='standard',\
		specmode='mfs',niter=myniter,nterms=mynterms,stokes=mystokes,outlierfile=outlierf,scales=multiscale,robust=robust)
if mynterms>1 and decon=='mtmfs':
	imagen=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.image.tt0'
	imagenpb=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.pb.tt0'
else:
	imagen=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.image'
	imagenpb=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.pb'
print 'Correcting for PB...'
os.system('rm -rf '+imagen+'.pbcor')
os.system('rm -rf '+imagen+'.pbcor.fits')
immath(imagename=[imagen,imagenpb],\
	expr='IM0/IM1',outfile=imagen+'.pbcor')
print 'Making fits image...'
exportfits(imagename=imagen+'.pbcor',fitsimage=imagen+'.pbcor.fits')
imagen=imagen+'.pbcor'
dofit=raw_input('Do you want to fit a point source?y or n-->')
if dofit=='y':
	flux,err,unit,freq,err_real=imfit_point(imagen,my_dir,'I')
	print 'Flux density of ',flux,' +/- ',err, unit
	print 'Local RMS in image is: ',err_real,' Jy'
dopsc=raw_input('Do you want to do phase selfcal?y or n-->')
if dopsc=='y':
	selfcal,scim=phselfcal(ms_name,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
obsDate,subband,combi,spw,outlierf,multiscale,robust,weighting)
	flux_sc,err_sc,unit_sc,freq_sc,err_real_sc=imfit_point(scim,my_dir)

if do_pol=='y':
	os.system('rm -rf '+my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV*')
	print 'Imaging polarization cube...'
	if mymask=='':
		tclean(vis=ms_name, imagename=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV',
			field='',spw=spw,interactive=True,cell=[mycell], imsize=myimsize,gain=0.1,
			weighting=weighting,threshold=mythreshold,specmode='mfs',niter=myniter,nterms=mynterms,
			stokes='IQUV',scales=multiscale,robust=robust,outlierfile=outlierf,deconvolver='clarkstokes',gridder='standard')
	else:
		tclean(vis=ms_name, imagename=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV',
			field='',mask=mymask,spw=spw,interactive=False,cell=[mycell], imsize=myimsize,gain=0.1,
			weighting=weighting,threshold=mythreshold,specmode='mfs',niter=myniter,nterms=mynterms,
			stokes='IQUV',scales=multiscale,robust=robust,outlierfile=outlierf,deconvolver='clarkstokes',gridder='standard')
	if mynterms>1:
		imagenpol=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV.image.tt0'
		imagenpolpb=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV.pb.tt0'
	else:
		imagenpol=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV.image'
		imagenpolpb=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_polcube_IQUV.pb'
	print 'Correcting for PB...'
	os.system('rm -rf '+imagenpol+'.pbcor')
	os.system('rm -rf '+imagenpol+'.pbcor.fits')
	immath(imagename=[imagenpol,imagenpolpb],
		expr='IM0/IM1',outfile=imagenpol+'.pbcor')
	print 'Making fits image...'
	exportfits(imagename=imagenpol+'.pbcor',fitsimage=imagenpol+'.pbcor.fits')
	print 'Extracting Seperate Stokes Images...'
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.I.')
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.Q.')
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.U.')
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.V.')
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP')
	#os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube_Q.flux')
	#os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP.pbcor')
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.PA')
	os.system('rm -rf '+my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.FP')
	imsubimage(imagename=imagenpol,outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.I',stokes='I')
	imsubimage(imagename=imagenpol,outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.Q',stokes='Q')
	imsubimage(imagename=imagenpol,outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.U',stokes='U')
	imsubimage(imagename=imagenpol,outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.V',stokes='V')
	fluxQ,errQ,unitQ,freqQ,errQ_real=imfit_point(my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.Q',my_dir,'Q')
	fluxU,errU,unitU,freqU,errU_real=imfit_point(my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.U',my_dir,'U')
	fluxI,errI,unitI,freqI,errI_real=imfit_point(my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.I',my_dir,'I')
	immath(outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP',mode='poli',
		imagename=[my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.Q',my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.U'],
		sigma='0.0Jy/beam')
	#immath(imagename=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube_IQUV.flux',
		#outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube_Q.flux',expr='IM0',stokes='Q')
	#immath(imagename=[my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP',my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube_Q.flux'],\
		#expr='IM0/IM1',outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP.pbcor')
	immath(outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.PA',mode='pola',
		imagename=[my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.Q',my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.U'],
		polithresh=str(errU_real)+'Jy/beam')
	immath(outfile=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.FP',mode='evalexpr',
		imagename=[my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.I',my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.Q',
		my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.U'],expr='sqrt((IM1^2+IM2^2)/IM0[IM0>'+str(3.*errI_real)+']^2)')
	print 'Making FITS files...'
	exportfits(imagename=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP',fitsimage=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.LP.fits')
	exportfits(imagename=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.PA',fitsimage=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.PA.fits')
	exportfits(imagename=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.FP',fitsimage=my_dir+target+'_'+obsdate+'_'+band+'_'+subband+'_polcube.FP.fits')
#writing imfit result to file
if dopsc=='n'and dofit=='y':
	print 'Writing imfit results to file...'
	resul_file=open(my_dir+target+'_'+obsDate+'_'+subband+'_imfit_results.txt','w')
	resul_file.write('#band freq flux err unit err_real\n')
	resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freq,flux,err,unit,err_real))
	resul_file.close()
elif dopsc=='y'and dofit=='y':
	print 'Writing imfit results to file...'
	resul_file=open(my_dir+target+'_'+obsDate+'_'+subband+'_imfit_results.txt','w')
	resul_file.write('#band freq flux err unit err_real\n')
	resul_file.write('Pre-selfcal:\n')
	resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freq,flux,err,unit,err_real))
	resul_file.write('Post selfcal:\n')
	resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freq_sc,flux_sc,err_sc,unit_sc,err_real_sc))
	resul_file.close()
###########################################

###########################################
#UVfitting
###########################################
if uv_fit=='T':
	phcen_rad=vishead(vis=ms_name,mode='get',hdkey='ptcs')
	phcen=au.rad2radec(phcen_rad[0]['r1'][0][0][0],phcen_rad[0]['r1'][1][0][0],hmsdms=True).replace(',','')
	if uv_initp != '':
		file_uv0=open(uv_initp)
		uv_num=sum(1 for line in file_uv0)-1
		file_uv0.close()
		initp_array=np.loadtxt(uv_initp)
		init_uv=[]
		if uv_num==1:
			init_uv.extend([initp_array[0],initp_array[1],initp_array[2]])
		else:
			for jj in range(0,uv_num):
				init_uv.extend([initp_array[jj,0],initp_array[jj,1],initp_array[jj,2]])
	else:
		init_uv=[]
		uv_num=1
	comp_uv=[]
	var_uv=[]
	for jj in range(0,uv_num):
		comp_uv.append('delta')
		var_uv.extend(['p['+str(3*jj)+'],p['+str(3*jj+1)+'],p['+str(3*jj+2)+']'])
	#
	print 'Performing UV fitting...'	
	stokes_param=mystokes
	comb=ms_name.strip('.ms')
	mstransform(vis=ms_name, outputvis=comb+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
	uvvis=comb+'_mstrans.ms'
	if uv_initp=='':
		fitfulluv=uvm.uvmultifit(vis=uvvis, spw=spw, column = "data",\
			uniform=False, model=comp_uv,stokes = stokes_param, var=var_uv, phase_center =phcen,\
			outfile = my_dir+'modelfit.dat', OneFitPerChannel=False ,\
			cov_return=False,finetune=False, method="levenberg")
	else:
		fitfulluv=uvm.uvmultifit(vis=uvvis, spw=spw, column = "data",\
			uniform=False, model=comp_uv,stokes = stokes_param, var=var_uv, phase_center =phcen,\
			outfile = my_dir+'modelfit.dat', OneFitPerChannel=False ,p_ini=init_uv,\
			cov_return=False,finetune=False, method="levenberg")
	src_uv_init=fitfulluv.result['Parameters'][2]
	src_uv_err=fitfulluv.result['Uncertainties'][2]
	print 'Flux density of ',src_uv_init,' +/- ',src_uv_err, 'Jy'
	print 'Writing uvfit results to file...'
	resuluv_file=open(my_dir+target+'_'+obsDate+'_'+subband+'_uvfit_results.txt','w')
	resuluv_file.write('{0} {1} {2} {3}\n'.format(band,subband,src_uv_init,src_uv_err))
	resuluv_file.close()
###########################################

print 'Cleaning up...'
os.system('rm -rf casa*.log')
os.system('rm -rf *.last')
print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'
#################################################
