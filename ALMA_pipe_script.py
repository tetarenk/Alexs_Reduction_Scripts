#################################################
#ALMA CASA Pipeline Script
#################################################
'''CASA script to use pipeline to do a quick preliminary reduction of ALMA
continuum data. Use when you only have the raw asdm data directly after a trigger.
INPUT: Parameters below.
OUTPUT: (1) Calibrated Split MS for full band (and each spw) -- [target]_[obsDate]_[band]_cal.ms
        (2) Continuum images for full band (and each epw) -- [target]_[obsDate]_[band]_clean1.image.(.tt0)
        (3) Fits point source and prints flux/rms to terminal
NOTES: - All output images & intermediate data products are put in current working directory.
Written by: Alex J. Tetarenko
Last Updated: Aug 2019
Works in CASA-5.4 now!

USAGE: start CASA pipeline with casapy-5.4 --pipeline, enter parameters below, run script
'''

import glob
import analysisUtils as au
import pipeline.recipereducer as recipereducer
from utils_CASA import imfit_point

##################
#Input Parameters
##################
asdm_name='/PATH_TO_ASDM/'
target='Maxi1535'
obsDate='Sep2017'
band='B3'
mythreshold='1mJy'
myimsize=[4096]
mycell=['0.015arcsec']
mynterms=1
decon='hogbom'#hogbom if nterms=1,mtmfs if nterms>1,multiscale if using scales
myniter=5000
mystokes='I'
multiscale=[]
weighting='natural'
robust=2.0
##################

#start of script
#(1) Import asdm to create MS and generate Weblog
print "Creating MS and generating weblog..."
h_init()
hifa_importdata(vis=[asdm_name],pipelinemode='automatic')
h_save()

#(2) open weblog
print "Opening Weblog..."
newest = max(glob.iglob('pipeline-*'), key=os.path.getctime)
os.system('google-chrome '+'/html/index.html')
raw_input('Please press enter to continue when you are done.')

#(3)Check that you have “flux.csv”, “*flagonline.txt”, 
#“*flagtemplate.txt”, and “*flagcmds.txt” in the working directory.
print 'Check you have all ingredients for pipeline...'
print 'MS,flagversions,flux.csv.'
os.system('ls')
os.system('sublime flux.csv &')
raw_input('Please press enter to continue when you are done.')

#(4) Run pipeline-- only calibration
print 'Running pipeline...'
ms_name=max(glob.iglob('*.ms'), key=os.path.getctime)
recipereducer.reduce(vis=[ms_name],procedure='procedure_hifacal.xml')
print 'Opening log...'
os.system('google-chrome pipeline-procedure_hifacal/html/index.html')
raw_input('Please press enter to continue when you are done.')

#(5) Check calibrated data
print 'Check calibrated data...'
print 'Listobs...'
listobs(ms_name, listfile='listfile_afterpipecal.txt')
raw_input('Please press enter to continue when you are done.')
spws=raw_input('Please enter science spws, e.g., 17,19,21,24--> ')
fields=raw_input('Please enter science fields, e.g., 1,2--> ')

print 'Plotting spectra...'
plotms(vis=ms_name,yaxis='amp',xaxis='channel',spw=spws,field=fields, \
	avgtime='1e8',avgscan=True,coloraxis='spw',overwrite=True,datacolumn='corrected')
raw_input('Please press enter to continue when you are done.')

print 'Splitting to make target MS...'
split(vis=ms_name,outputvis=target+'_'+obsDate+'_'+band+'_cal.ms',spw=spws,\
	datacolumn='corrected',field=fields)
ms_name=target+'_'+obsDate+'_'+band+'_cal.ms'

print 'Cleaning...'
tclean(vis=ms_name,imagename=target+'_'+obsDate+'_'+band+'_clean1',specmode='mfs',\
	imsize=myimsize,cell=mycell,spw='',gain=0.1,deconvolver=decon,gridder='standard',niter=myniter,\
	nterms=mynterms,stokes=mystokes,scales=multiscale,robust=robust,\
	weighting=myweighting,interactive=True,threshold=mythreshold,\
	mask='')
if mynterms >1:
	im=target+'_'+obsDate+'_'+band+'_clean1.image.tt0'
else:
	im=target+'_'+obsDate+'_'+band+'_clean1.image'
print 'Fitting point source to image...'
fluxbb,errbb,unitbb,freqbb,errbb_real=imfit_point(im, os.getcwd()+'/',mystokes)
print 'Flux density  of ',fluxbb,' +/- ',errbb, unitbb
print 'Local RMS is: ',errbb_real,' Jy'

print 'Cleaning up...'
os.system('rm -rf casa*.log')
os.system('rm -rf *.last')
os.system('rm -rf *.png')
print '**************************'
print 'The script is finished.'
print '**************************'