#######################################
# JCMT SCUBA-2 Timing Script
#######################################
'''Fits a 2D Gaussian to every plane of a set of JCMT fits timing cubes to get a high time resuloution light curve.
This script uses a simple least squares algorithm implmented by astropy models to do the fitting quickly.
INPUT: my_dir: Output directory
       shortmaps: List of shortmap FITS cubes from STARLINK for each scan
       w: Beam Parameters from STARLINK; [BMAJ (arcsec), BMIN (arcsec), BPA (deg)]
       cal_im: (optional) Calibrator FITS image for directly fitting beam size
       ranges/ranges_cal: pixel ranges to search for source/calibrator in images
       mjd_file: data file with mosaiced cubes MJD times
OUTPUT: Light curve data file (my_dir/JCMT_timing_fit_mJy.txt)
	Light curve plot (my_dir/target_lightcurve.png)
NOTES: - The uncertainties output are simply from the least squares covariance matrix.
       - The beam is fixed in target fitting, so only flux and position vary (pointing can change from plane
         to plane in shortmaps cube so dont fix position!).
       - To get MJD file: (1) Open a mosaiced shortmaps cube in gaia, (2) click on central source pixel,
         (3) a window will pop up showing lightcurve, (4) you can extract the light curve for this pixel
         and save to an txt file. MAKE SURE TO USE SDF not FITS version to do this!

Written by: Alex J. Tetarenko
Last Updated: Apr 22, 2019'''

#needed packages
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from astropy.coordinates import Angle
from astropy import wcs
from astropy import units as u
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator


print 'Reading in parameters...'
####################################
#User input
####################################
#input/output directory
my_dir=''
#list of shortmaps fits cubes from starlink reduction for all scans
shortmaps=[my_dir+'shortmap_cube_cal.fits']
#calibrator scan fits file for beam fitting
cal_im=my_dir+'cal.fits'
#size of gaussian beam in arcsec,pa in deg from starlink cal reduction
w=[15.54,14.40,177.31]
#boxes to search for source x1,x2,y1,y2 in all scans
ranges=[[125,137,118,130]]
ranges_cal=[132,149,116,132]
#file with MJDs from mosaiced shortmaps cube
mjd_file=my_dir+'mjd_full.txt'
#lightcurve binsize
interv='5sec'
#########################################



#########################################
print 'Starting fitting process...'''

print 'Reading in fits headers...'
pixsize=Angle(abs(getheader(cal_im)['CDELT1'])*u.degree).arcsec
header_cal=getheader(cal_im)
wmapcal=wcs.WCS(header_cal)

cal_fit=raw_input('Do you want to fit beam size from calibrator? y or n?--> ')
if cal_fit=='y':
	#fit calibrator to get beam size
	data_psf=getdata(cal_im,0)[0,ranges_cal[2]:ranges_cal[3],ranges_cal[0]:ranges_cal[1]]
	print "Fitting calibator scan..."
	param_guess=[1.,data_psf.shape[0]/2.,data_psf.shape[0]/2.,w[0]/(2.*pixsize),w[1]/(2.*pixsize),w[2]]#amplitude,X,Y,width_x,width_y,rota
	psf_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':False,'y_stddev':False,'theta':True}
	M=models.Gaussian2D(param_guess[0],param_guess[1],param_guess[2],param_guess[3],\
		param_guess[4],param_guess[5],fixed=psf_fixed_params)
	x=np.arange(0,len(data_psf[0,:]))
	y=np.arange(0,len(data_psf[:,0]))
	X, Y = np.meshgrid(x, y)
	Z=data_psf
	lmf=fitting.LevMarLSQFitter()
	res=lmf(M,X,Y,Z)
	#print results
	print '\n'
	print 'Calibrator scan fit results:'
	print 'Cal Flux = ','%.2f'%(res.parameters[0]*1e3),' +/- ','%.2f'%(1e3*np.sqrt(np.diag(lmf.fit_info['param_cov']))[0]), 'mJy'
	print 'Beam Maj (arcsec) = ','%.2f'%(res.parameters[3]*2.*pixsize)
	print 'Beam Min (arcsec) = ','%.2f'%(res.parameters[4]*2.*pixsize)
	print 'Beam PA (deg) = ','%.2f'%(res.parameters[4])
	print 'Starlink Beam Maj, Min, PA = ',w[0],w[1],w[2],'\n'
	outfile=open(my_dir+'cal_fits_timing_res.txt','w')
	outfile.write('Calibrator scan fit results:\n')
	outfile.write('\n')
	outfile.write('Cal Flux = '+'%.2f'%(res.parameters[0]*1e3)+' +/- '+'%.2f'%(1e3*np.sqrt(np.diag(lmf.fit_info['param_cov']))[0])+'mJy\n')
	outfile.write('Beam Maj (arcsec) = '+'%.2f'%(res.parameters[3]*2.*pixsize)+'\n')
	outfile.write('Beam Min (arcsec) = '+'%.2f'%(res.parameters[4]*2.*pixsize)+'\n')
	outfile.write('Beam PA (deg) = '+'%.2f'%(res.parameters[4])+'\n')
	outfile.write('Input Starlink Beam Maj, Min, PA = '+str(w[0])+','+str(w[1])+','+str(w[2])+'\n')
	outfile.write('\n')
	outfile.close()
	#show results
	print 'Plotting result...'
	fig=plt.figure()
	plt.rc('xtick',color='w')
	plt.rc('ytick',color='w')
	mpl.rcParams['xtick.major.width']=1.5
	mpl.rcParams['xtick.major.size']=7
	mpl.rcParams['xtick.minor.size']=5
	mpl.rcParams['xtick.direction']='in'
	mpl.rcParams['ytick.direction']='in'
	mpl.rcParams['ytick.right']='on'
	mpl.rcParams['xtick.top']='on'
	mpl.rcParams['axes.edgecolor']='w'
	ax=plt.subplot(111,projection=wmapcal.celestial)
	ax.imshow(data_psf,origin='lower')
	ax.contour(res(X,Y),colors='w',levels=res.parameters[0]*np.array([0.2,0.5,0.8]))
	ax.coords['ra'].set_axislabel('Right Ascension')
	ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
	ax.coords['ra'].set_major_formatter('hh:mm:ss')
	ax.set_title('Calibrator Fit')
	plt.savefig(my_dir+'calibrator_fit.png')
	plt.show()
	raw_input('Please press enter when ready to continue.')

#fit target by fixing beam size in each plane of shortmaps cube
print 'Beginning to fit target...'
fix=raw_input('Do you want to fix the beam to the starlink values or the calibrator fit (if selected above)? star or cal?--> ')


flux=[]
error=[]
for i in range(0,len(shortmaps)):
	print "Fitting Target in scan ",i+1,'...'
	fitsim=shortmaps[i]
	if fix=='cal':
		param_guess=[1.,(abs(ranges[i][1]-ranges[i][0])/2.),(abs(ranges[i][3]-ranges[i][2])/2.),res.parameters[3],res.parameters[4],res.parameters[5]]
	elif fix=='star':
		param_guess=[1.,(abs(ranges[i][1]-ranges[i][0])/2.),(abs(ranges[i][3]-ranges[i][2])/2.),w[0]/(2.*pixsize),w[1]/(2.*pixsize),w[2]]
	tar_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':True,'y_stddev':True,'theta':True}
	data_target=getdata(fitsim,0)
	shape_data=data_target.shape
	for k in range(0,shape_data[0]):
		data_target=getdata(fitsim,0)[k,ranges[i][2]:ranges[i][3],ranges[i][0]:ranges[i][1]]
		x=np.arange(0,len(data_target[0,:]))
		y=np.arange(0,len(data_target[:,0]))
		X, Y = np.meshgrid(x, y)
		Z=data_target
		M=models.Gaussian2D(param_guess[0],param_guess[1],param_guess[2],param_guess[3],\
			param_guess[4],param_guess[5],fixed=tar_fixed_params)
		lmf=fitting.LevMarLSQFitter()
		res=lmf(M,X,Y,Z)
		flux.append(res.parameters[0]*1e3)
		error.append(1e3*np.sqrt(np.diag(lmf.fit_info['param_cov']))[0])

#read in MJDs
mjd_list=np.loadtxt(mjd_file)
#mjds now off by 12 hours coming out of Gaia?
mjd=mjd_list[:,0]+0.5
if len(mjd)!=len(flux):
	raise Exception('MJDs dont match fluxes, please double check your input files')

#write results to file
print 'Writing results to file...'
file_res=open(my_dir+'JCMT_timing_fit_mJy_'+interv+'.txt','w')
for j in range(0,len(mjd)):
	if flux[j] != np.nan:
		file_res.write(' {0} {1} {2}\n'.format(mjd[j],flux[j],error[j]))
file_res.close()

#show results
plt.rcdefaults()
print 'Plotting result...'
fig=plt.figure()
font={'family':'serif','weight':'bold','size' : '14'}
rc('font',**font)
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
ax=plt.subplot(111)
ax.errorbar(mjd,flux,yerr=error,marker='o',color='m',ls='')
#ax.set_title('Target Light Curve')
ax.set_ylabel('Flux Density (mJy)')
ax.set_xlabel('Time (HH:MM)')
locator=mdates.MinuteLocator(interval=30)
locator2=mdates.MinuteLocator(interval=5)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_minor_locator(locator2)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.setp(ax.get_xticklabels(),rotation=45,horizontalalignment='right')
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis='x',which='major', labelsize=9,length=5,width=1.5,top='on',bottom='on',pad=7)
ax.tick_params(axis='x',which='minor', labelsize=9,length=3.5,width=1.,top='on',bottom='on',pad=7)
ax.tick_params(axis='y',which='major', labelsize=10,length=5,width=1.5,right='on',left='on')
ax.tick_params(axis='y',which='minor', labelsize=9,length=3.5,width=1.,right='on',left='on',pad=7)
plt.savefig(my_dir+'target_lightcurve_'+interv+'.png',bbox_inches='tight')
plt.show()


print '\n'
print '*********************************************************'
print 'Script is finished. Please inspect output data products.'
print '*********************************************************'
