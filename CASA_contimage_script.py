#################################################
#CASA Continuum Imaging/Selfcal Script
#################################################
'''CASA script to be used for imaging and selfcal of VLA/SMA continuum data
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) continuum image  -- [target]_[obsDate]_[band]_[subband]_clean1.image(.tt0).pbcor
        (2) File of flux densities from image/UV plane fitting -- [target]_[obsDate]_[subbband]_imfit(uvfit)_results.txt
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
       - All output images are also converted to fits format (just append .fits to end of images 1 above)
       - This script images and selfcals an already calibrated CASA MS and fits a point source in the image/uv plane.
Written by: Alex J. Tetarenko
Last Updated: Jan 5 2016'''

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
param_dir_file='/home/ubuntu/CASA_reduction_scripts/params_cont_img.txt'
print 'You have set your param file to ', param_dir_file
print 'Please make sure all parameters are correct, they will change for each data set!\n'

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
obsDate=data_params.obsDate
target=data_params.target
spw=data_params.spw
band=data_params.band
subband=data_params.subband
ref_ant=data_params.ref_ant
#general image params
use_auto=data_params.use_auto
mythreshold=data_params.mythreshold
myimsize=data_params.myimsize
mycell=data_params.mycell
mynterms=data_params.mynterms
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
	clean(vis=ms_name, imagename=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1',field='',spw=spw,interactive=T,\
		cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
		mode='mfs',niter=0,nterms=mynterms,stokes=mystokes,outlierfile=outlierf,multiscale=multiscale,robust=robust)
else:
	os.system('rm -rf '+my_dir+target+'_'+date+'_'+band_low+'_clean1*')
	print 'Cleaning...'
	clean(vis=ms_name, imagename=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1',field='',mask=mymask,spw=spw,interactive=F,\
		cell=[mycell], imsize=myimsize,gain=0.1,weighting=weighting,threshold=mythreshold,\
		mode='mfs',niter=myniter,nterms=mynterms,stokes=mystokes,outlierfile=outlierf,multiscale=multiscale,robust=robust)
if mynterms>1:
	imagen=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.image.tt0'
else:
	imagen=my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.image'
print 'Correcting for PB...'
os.system('rm -rf '+imagen+'.pbcor')
os.system('rm -rf '+imagen+'.pbcor.fits')
immath(imagename=[imagen,my_dir+target+'_'+obsDate+'_'+band+'_'+subband+'_clean1.flux'],\
	expr='IM0/IM1',outfile=imagen+'.pbcor')
print 'Making fits image...'
exportfits(imagename=imagen+'.pbcor',fitsimage=imagen+'.pbcor.fits')
imagen=imagen+'.pbcor'
dofit=raw_input('Do you want to fit point source?y or n-->')
if dofit=='y':
	flux,err,unit,freq,err_real=imfit_point(imagen,my_dir)
	print 'Flux density of ',flux,' +/- ',err, unit
	print 'Local RMS in image is: ',err_real,' Jy'
dopsc=raw_input('Do you want to do phase selfcal?y or n-->')
if dopsc=='y':
	selfcal,scim=phselfcal(ms_name,mycell,mynterms,myimsize,mythreshold,ref_ant,my_dir,target,\
obsDate,subband,combi,spw,outlierf,multiscale,robust,weighting)
	flux_sc,err_sc,unit_sc,freq_sc,err_real_sc=imfit_point(scim,my_dir)

#writing imfit result to file
if dopsc=='n' and dofit=='y':
	print 'Writing imfit results to file...'
	resul_file=open(my_dir+target+'_'+obsDate+'_'+subband+'_imfit_results.txt','w')
	resul_file.write('{0} {1} {2} {3} {4} {5}\n'.format(band,freq,flux,err,unit,err_real))
	resul_file.close()
elif dopsc=='y' and dofit=='y':
	print 'Writing imfit results to file...'
	resul_file=open(my_dir+target+'_'+obsDate+'_'+subband+'_imfit_results.txt','w')
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
	print 'Performing UV fitting...'
	comp_uv='delta'
	stokes_param=mystokes
	comb=ms_name.strip('.ms')
	mstransform(vis=ms_name, outputvis=comb+'_mstrans.ms', combinespws=True, spw='',datacolumn='data')
	uvvis=comb+'_mstrans.ms'
	fitfulluv=uvm.uvmultifit(vis=uvvis, spw=spw, column = "data", \
		uniform=False, model=[comp_uv],stokes = stokes_param, \
		var=['p[0],p[1],p[2]'],OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
	src_uv_init=fitfulluv.result['Parameters'][2]
	src_uv_err=fitfulluv.result['Uncertainties'][2]
	print 'Flux density of ',src_uv_init,' +/- ',src_uv_err, 'Jy'
	print 'Writing uvfit results to file...'
	resuluv_file=open(my_dir+target+'_'+obsDate+'_'+subband+'_uvfit_results.txt','w')
	resuluv_file.write('{0} {1} {2} {3}\n'.format(band,subband,src_uv_init,src_uv_err))
	resuluv_file.close()
###########################################

print 'Cleaning up...'
os.system('rm -rf *.log')
os.system('rm -rf *.last')
print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'
#################################################
