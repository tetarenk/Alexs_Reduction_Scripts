############################################################
#SMA uvfits data to CASA MS
############################################################
'''CASA script to be used to convert individual spw uvfits files of SMA data to CASA MS.
INPUT: Parameters defined in section below, uvfits files of corrected SMA data sets.
OUTPUT: (1) MS for each SMA sideband -- [target]+_lsb.ms,[target]+_usb.ms
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
	   - This script is intended to be used on output from miriad2fits.sh (step1/2). 
	   - This is step 2/2 in the long MIRIAD based process to convert SMA data to a CASA MS.
Written by: Alex J. Tetarenko
Last Updated: Jan 4 2017'''

#needed packages
import os

############################
#Defining variables section
############################
my_dir='/Users/atetarenk/Desktop/CASA_reduction_scripts/'
#from miriad scripts
#dt=150519
#rx=0
lsb_fits_name=''
usb_fits_name=''
beg=1
end=50
#target source name
target='V404Cyg'
############################

#make directory for ms files
os.system('rm -rf '+my_dir+'spw_ms_files')
os.system('mkdir '+my_dir+'spw_ms_files')

#import each sw fits file into casa for both lower and upper sideband
myfilesl = []
myfilesu = []
for i in range(beg,end+1):
    fitsfilel = my_dir+'spw_fits_files/'+lsb_fits_name+'.spw_'+i+'.fits'
    fitsfileu = my_dir+'spw_fits_files/'+usb_fits_name+'.spw_'+i+'.fits'
    msfilel = my_dir+'spw_ms_files'+target+'_lsb_spw_'+i+'.ms'
    msfileu = my_dir+'spw_ms_files'+target+'_usb_spw_'+i+'.ms'
    importuvfits(fitsfile=fitsfilel,
                 vis=msfilel)
    importuvfits(fitsfile=fitsfileu,
                 vis=msfileu)
    myfilesl.append(msfilel)
    myfilesu.append(msfileu)
#concatenate all individual spw MSs for each sideband
concat(vis=myfilesl,concatvis=my_dir+target+'_lsb.ms',timesort=True)
concat(vis=myfilesu,concatvis=my_dir+target+'_usb.ms',timesort=True)
#initialize scratch columns
clearcal(vis=my_dir+target+'_lsb.ms')
clearcal(vis=my_dir+target+'_usb.ms')
#do a listobs to check everything is correct
listobs(vis=my_dir+target+'_lsb.ms',verbose=False,listfile=my_dir+target+'_lsb_listfile.txt')
listobs(vis=my_dir+target+'_usb.ms',verbose=False,listfile=my_dir+target+'_usb_listfile.txt')
os.system('pluma '+my_dir+target+'_lsb_listfile.txt &')
os.system('pluma '+my_dir+target+'_usb_listfile.txt &')
