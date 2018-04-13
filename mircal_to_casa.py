'''Convert calibrated raw SMA data (MIR format) to CASA for imaging
INPUT: Path to uvfits files, and lists of source names and recievers input on line 12
OUTPUT: Concatenated CASA MS

REMINDER: (1) Need to download MIRFITStoCASA.py code, which is a version
          of the importfits task that propagates MIR weights correctly
          (2) Before this script is run, calibrated data must be written to
          UVFITS files in MIR, following instructions here:
          https://www.cfa.harvard.edu/rtdc/SMAdata/process/casa/convertcasa/
'''

import sys
#path to MIRFITStoCASA.py code
sys.path.append('/mnt/bigdata/tetarenk/CASA_reduction_scripts')
import MIRFITStoCASA


###############################
#User input
###############################
my_dir='/mnt/bigdata/tetarenk/SMA_maxi1820/uvfits_files/'
sources=['MAXIJ1820+070']
recievers=['230','240']
ms_name=my_dir+'MAXI1820_calibrated_SMA.ms'
###############################

allNames=[]
for sou in sources:
	for rx in recievers:
		for sb in ['L','U']:
			for i in ['1','2','3','4']:
				name=my_dir+sou+'_'+sb+'_S'+i+'_RX'+rx
				print('------converting '+name+' ....')
				MIRFITStoCASA.MIRFITStoCASA(UVFITSname=name+'.UVFITS', MSname=name+'.ms')
				allNames.append(name+'.ms')
print 'Concatenating all MSs...'
concat(vis=allNames,concatvis =ms_name,timesort=True)
print 'Running listobs as a check...'
listobs(ms_name,listfile=my_dir+'listobs.txt')
os.system('pluma '+my_dir+'listobs.txt &')
