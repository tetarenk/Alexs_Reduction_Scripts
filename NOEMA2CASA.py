##########################################
#NOEMA/PdBI data in CASA Script
##########################################
'''CASA script to be used for importing, viewing visibilities, and imaging calibrated NOEMA/PdBI data in CASA.
INPUT: Parameters defined in section below, uvfits files of calibratedNOEMA/PdBI data sets
OUTPUT: (1) MS for each NOEMA/PdBI data set
		(2) Continuum images for each data set
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
	   - All output images are also converted to fits format (just append .fits to end of images above)
	   - This script is intended to be used with calibrated data that has been converted to uvfits format
	     following instructions here; http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf
Written by: Alex J. Tetarenko
Last Updated: November 13 2017'''


#needed packages
import os
import glob
from ekoch_casa_tools import set_imagermode,has_field,set_cellsize,set_imagesize,find_expected_beams,getBaselinePercentile,get_mosaic_info

############################
#Defining variables section
############################
my_dir='/Users/atetarenk/Desktop/CASA_reduction_scripts/NOEMA_uvfits_files/'
target='V404Cyg'
obsDate='jun26'
band='1mm'
mycell='0.5arcsec'
mythresh='1mJy'
myimsize=64
use_auto='F'
############################

#get list of uvfits files for importing
print 'Getting list of uvfits files...'
uvfits_files=[]
for name in glob.glob(my_dir+'*'):
    uvfits_files.append(name)

#import,plot visibilities, and image data set in casa
print 'Importing uvfits files...'
for i in range(0,len(uvfits_files)):
	print 'Importing ', uvfits_files[i],' ...'
	msname=uvfits_files[i].strip('.fits')
	importuvfits(fitsfile=uvfits_files[i],vis=my_dir+'noema_'+target+'_'+msname+'.ms')
	print 'Listobs...'
	listobs(my_dir+'noema_'+target+'_'+msname+'.ms',listfile=my_dir+'noema_'+target+'_'+msname+'_listobs.txt')
	os.system('pluma 'my_dir+'noema_'+target+'_'+msname+'_listobs.txt &')
	raw_input('Please press enter when ready to continue.')
	print 'Plotting visibility data...'
	plotms(vis=my_dir+'noema_'+target+'_'+msname+'.ms',xaxis="time",yaxis="amp",ydatacolumn='data',iteraxis="antenna")
	raw_input('Please press enter when ready to continue.')
	print 'Interactive cleaning...'
	if use_auto=='T':
		myimsize=set_imagesize(my_dir+'noema_'+target+'_'+msname+'.ms','0','0')
		mycell=set_cellsize(my_dir+'noema_'+target+'_'+msname+'.ms','0')
	clean(vis=my_dir+'noema_'+target+'_'+msname+'.ms', imagename=my_dir+'noema_'+target+'_'+obsDate+'_'+band,\
		field='', interactive=True, cell=[mycell], imsize=myimsize,gain=0.1,weighting='natural',\
		threshold=mythresh,mode='mfs',myniter=0)
	imview(my_dir+'noema_'+target+'_'+obsDate+'_'+band+'.image')
	raw_input('Please press enter when ready to continue.')
