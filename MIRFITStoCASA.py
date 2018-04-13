def MIRFITStoCASA(UVFITSname='None', MSname='None', verbose=True):
	#UVFITSname='GJ322TRACK1_230SUB_U.FITS'
	#MSname='GJ322TRACK1_230SUB_U.ms'
	#verbose=True
	
	############
	## Import needed packages
	############
	import pyfits as pf
	import os
	from casac import casac
	import importuvfits
	import numpy as np
	import listobs
	#import casac.casac as casacc
	if os.path.exists(MSname):
		os.system('rm -rf '+MSname)
	tb = casac.table()
	ms = casac.ms()
	
	############
	## Import UVFITS into CASA MS using CASA's importuvfits
	############
	importuvfits.importuvfits(UVFITSname, vis=MSname)
	
	############
	#Read wtscale from UVFITS header input
	############
	head=pf.getheader(UVFITSname)
	if float(pf.__version__[:3])<3.2:
		# THIS WORKS WITH CASA 4.7.2 and earlier
		nhistlines=len(head.get_history())
		hist=[]
		for i in np.arange(nhistlines):
			hist.append(str(head.get_history()[i]))
			histmerged='\n'.join(hist)
	else:
		# THIS WORKS WITH CASA 5.0.0 and later
		histmerged=str(head['HISTORY'])
	begstrind=int((histmerged.find('AIPS WTSCAL ='))+len('AIPS WTSCAL =')+1)
	endstrind=int(len(histmerged))
	endlineind=histmerged.find('\n', begstrind, endstrind)	
	wtscale=float(histmerged[begstrind:endlineind])
	
	############
	## Correcting weights using basically 'scaleweights' function in CASA's analysisUtils package
	############
	# Find number of data description IDs
	#try:
	tb.open(MSname,nomodify=False)
	#except:
	#	print "ERROR: failed to open ms tool on file "+MSname
	#        tb.close()
	#        return(3)
	recw = tb.getcol('WEIGHT')
	recw_sp = tb.getcol('WEIGHT_SPECTRUM')
	if verbose:
		print 'Multiplying weights in the MS by a factor '+str(wtscale)
	recw_sp*=wtscale
	# recw*=wtscale
	# Although importuvfits sets WEIGHT=SUM(WEIGHT_SPECTRUM) I don't believe that should be the case, e.g. it does not seem to be the case in ALMA data.
	# I found this also causes problems when channel averaging to form a 'continuum' in CASA. The weights in the continuum should be = SUM(WEIGHT_SPECTRUM) but they end up being equal to SUM(WEIGHT) in plotms, and that's incorrect.
	# So setting WEIGHT=WEIGHT_SPECTRUM as in ALMA data but *NOT* as done by importuvfits. 
	recw=recw_sp[:,0,:]
	tb.putcol('WEIGHT', recw)
	tb.putcol('WEIGHT_SPECTRUM', recw_sp)
	tb.putcol('SIGMA', np.sqrt(1.0/recw))
	tb.close()
	
	return 1
