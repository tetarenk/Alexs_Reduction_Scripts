#######################################
# JCMT SCUBA-2 Flux Fitting Code
#######################################
'''Fits a 2D Gaussian to JCMT fits images using MCMC or a simple least squares algorithm.
INPUT: my_dir: Output directory
       fitsim: Target FITS image
       w: Beam Parameters from STARLINK; [BMAJ (arcsec), BMIN (arcsec), BPA (deg)]
       cal_im: (optional) Calibrator FITS image for directly fitting beam size
       ranges/ranges_cal: pixel ranges to search for source/calibrator in FITS images
OUTPUT: Results file (my_dir/fit_results.txt)
	Resulting images/models (my_dir/calibrator_fit.png,my_dir/target_fit.png)
NOTES: It is recommended to use image RMS as an uncertainty measurment on flux.
The beam is fixed in target fitting, so only flux and position vary.

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
from astropy.time import Time
import emcee
from astropy.utils.console import ProgressBar
import scipy.stats as ss

def confidenceInterval(y,sig):
	median=np.median(y)
	pct15=np.percentile(y,15)
	pct85=np.percentile(y,85)
	list1=np.array([median,median-pct15,pct85-median])
	return list1
def lp(p,data,error,fixp,guess,pixsize):
	#amp in mJy, xx and yy in pixels, bmaj and bmin in arcsec, bpa in deg
	amp,xx,yy,bmaj,bmin,bpa=p[0],p[1],p[2],p[3],p[4],p[5]
	mod0=models.Gaussian2D(amp,xx,yy,bmaj/(2.*pixsize),bmin/(2.*pixsize),bpa)
	xval=np.arange(0,len(data[0,:]))
	yval=np.arange(0,len(data[:,0]))
	Xval, Yval = np.meshgrid(xval, yval)
	mod1=mod0(Xval,Yval)
	re=-0.5*np.nansum(np.log(2*np.pi*error**2))-np.nansum((mod1-data)**2/(2*error**2))
	prior=prior_func(p,fixp,guess)
	return(re+prior)
def prior_func(pval,fixp,guess):
	pv=[]
	for i in range(0,len(fixp)):
		if fixp[i]==True:
			pv.append(guess[i])
		elif fixp[i]==False:
			pv.append(pval[i])
		else:
			raise ValueError('The fixed param array values can only be True or False')
	p=np.array(pv)
	amp,xx,yy,bmaj,bmin,bpa=p[0],p[1],p[2],p[3],p[4],p[5]
	prior=0.0
	prior += ss.norm.logpdf(amp,loc=guess[0],scale=10.)+ss.uniform.logpdf(amp,loc=0,scale=10e3)
	prior += ss.norm.logpdf(xx,loc=guess[1],scale=5.)+ss.uniform.logpdf(xx,loc=0,scale=2.*guess[1])
	prior += ss.norm.logpdf(yy,loc=guess[2],scale=5.)+ss.uniform.logpdf(xx,loc=0,scale=2.*guess[2])
	prior += ss.norm.logpdf(bmaj,loc=guess[3],scale=5.)+ss.uniform.logpdf(bmaj,loc=0,scale=30)
	prior += ss.norm.logpdf(bmin,loc=guess[4],scale=5.)+ss.uniform.logpdf(bmin,loc=0,scale=30)
	prior += ss.norm.logpdf(bpa,loc=guess[5],scale=10.)
	if np.isnan(prior):
		return(-np.inf)
	return(prior)
def mcmc_fit(data,error,guess,fixp,nBurn,nSample,flag,pixsize):
	ndim=6
	nwalkers=ndim*10
	p0=np.zeros((nwalkers,ndim))
	for i in np.arange(ndim):
		if fixp[i]==True:
			p0[:,i]=guess[i]
		elif fixp[i]==False:
			p0[:,i]=(((np.random.randn(nwalkers))*0.01)+guess[i])
	sampler = emcee.EnsembleSampler(nwalkers,ndim,lp,args=[data,error,fixp,guess,pixsize])
	with ProgressBar(nBurn) as bar:
		for i,result in enumerate(sampler.sample(p0,iterations=nBurn)):
			bar.update()
	pos,prob,state=result[0],result[1],result[2]
	sampler.reset()
	with ProgressBar(nSample) as bar:
		for i,result in enumerate(sampler.sample(p0,iterations=nSample)):
			bar.update()
	amp=confidenceInterval(sampler.flatchain[:,0],1)
	xx=confidenceInterval(sampler.flatchain[:,1],1)
	yy=confidenceInterval(sampler.flatchain[:,2],1)
	bmaj=confidenceInterval(sampler.flatchain[:,3],1)
	bmin=confidenceInterval(sampler.flatchain[:,4],1)
	bpa=confidenceInterval(sampler.flatchain[:,5],1)
	if flag=='y':
		plt.rcdefaults()
		fig=plt.figure(figsize=(10,3))
		for i in range(0, ndim):
			plt.subplot(1,6,i+1)
			patches=plt.hist(sampler.flatchain[:,i],bins=100)
		fig.subplots_adjust(hspace=0.5)
		plt.show()
		raw_input('Press enter to continue')
		fig=plt.figure(figsize=(10,3))
		for i in range(0, ndim):
			plt.subplot(1,6,i+1)
			plt.plot(sampler.chain[:,:,i].T)
		fig.subplots_adjust(hspace=0.5)
		plt.show()
	return(amp,xx,yy,bmaj,bmin,bpa)
	


print 'Reading in parameters...'
####################################
#User input
####################################
#input/output directory
my_dir=''
#fits images
fitsim=my_dir+'target.fits'
cal_im=my_dir+'cal.fits'
#flux guesses for cal and target in mJy
flux_guess=[5000,700]
#size of gaussian beam in arcsec,pa in deg
w=[15.54,14.40,177.31]
#boxes to search for source x1,x2,y1,y2
ranges=[45,55,44,56]
ranges_cal=[132,149,116,132]
#fit algorithm -- mcmc or lsq
fit_alg='mcmc'
outf=my_dir+'fit_results.txt'
#append to already existing outfile ('a') or create newfile ('w')
file_choice='w'
#########################################



#########################################
print 'Starting fitting process...'''

print 'Reading in fits headers...'
pixsize=Angle(abs(getheader(cal_im)['CDELT1'])*u.degree).arcsec
header_cal=getheader(cal_im)
header_target=getheader(fitsim)
wmapcal=wcs.WCS(header_cal)
wmaptar=wcs.WCS(header_target)
outfile=open(outf,file_choice)
outfile.write('JCMT SCUBA-2 Fit Results-->\n')
outfile.write('\n')

cal_fit=raw_input('Do you want to fit beam size from calibrator? y or n?--> ')
if cal_fit=='y':
	#fit calibrator to get beam size
	data_psf=getdata(cal_im,0)[0,ranges_cal[2]:ranges_cal[3],ranges_cal[0]:ranges_cal[1]]
	var_psf=getdata(cal_im,1)[0,ranges_cal[2]:ranges_cal[3],ranges_cal[0]:ranges_cal[1]]
	print "Fitting calibrator scan..."
	param_guess=[flux_guess[0]/1e3,data_psf.shape[0]/2.,data_psf.shape[0]/2.,w[0],w[1],w[2]]#amplitude,X,Y,width_x,width_y,rota
	psf_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':False,'y_stddev':False,'theta':True}
	x=np.arange(0,len(data_psf[0,:]))
	y=np.arange(0,len(data_psf[:,0]))
	X, Y = np.meshgrid(x, y)
	Z=data_psf
	if fit_alg=='lsq':
		M=models.Gaussian2D(param_guess[0],param_guess[1],param_guess[2],param_guess[3]/(2.*pixsize),\
			param_guess[4]/(2.*pixsize),param_guess[5],fixed=psf_fixed_params)
		lmf=fitting.LevMarLSQFitter()
		res=lmf(M,X,Y,Z)
		amp_cal=[res.parameters[0]*1e3,1e3*np.sqrt(np.diag(lmf.fit_info['param_cov']))[0]]
		xx_cal=res.parameters[1]
		yy_cal=res.parameters[2]
		bmaj=res.parameters[3]*2.*pixsize
		bmin=res.parameters[4]*2.*pixsize
		bpa=res.parameters[5]
	elif fit_alg=='mcmc':
		amp_cal0,xx_cal0,yy_cal0,bmaj0,bmin0,bpa0=mcmc_fit(data_psf*1e3,var_psf*1e3,[param_guess[0]*1e3,param_guess[1],param_guess[2],param_guess[3],\
			param_guess[4],param_guess[5]],[False,False,False,False,False,False],500,1500,'n',pixsize)
		amp_cal=[amp_cal0[0],np.sqrt(amp_cal0[1]**2+amp_cal0[2]**2)]
		xx_cal=xx_cal0[0]
		yy_cal=yy_cal0[0]
		bmaj=bmaj0[0]
		bmin=bmin0[0]
		bpa=bpa0[0]
	#print results
	print '\n'
	print 'Calibrator scan fit results:'
	print 'Cal Flux = ','%.2f'%(amp_cal[0]),' +/- ','%.2f'%(amp_cal[1]), 'mJy'
	print 'Beam Maj (arcsec) = ','%.2f'%(bmaj)
	print 'Beam Min (arcsec) = ','%.2f'%(bmin)
	print 'Beam PA (deg) = ','%.2f'%(bpa)
	print 'Starlink Beam Maj, Min, PA = ',w[0],w[1],w[2],'\n'
	outfile.write('Calibrator scan fit results:\n')
	outfile.write('\n')
	outfile.write('Cal Flux = '+'%.2f'%(amp_cal[0])+' +/- '+'%.2f'%(amp_cal[1])+'mJy\n')
	outfile.write('Beam Maj (arcsec) = '+'%.2f'%(bmaj)+'\n')
	outfile.write('Beam Min (arcsec) = '+'%.2f'%(bmin)+'\n')
	outfile.write('Beam PA (deg) = '+'%.2f'%(bpa)+'\n')
	outfile.write('Input Starlink Beam Maj, Min, PA = '+str(w[0])+','+str(w[1])+','+str(w[2])+'\n')
	outfile.write('\n')

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
	#ax.contour(res(X,Y),colors='w',levels=res.parameters[0]*np.array([0.2,0.5,0.8]))
	res2=models.Gaussian2D(amp_cal[0],xx_cal,yy_cal,bmaj/(2.*pixsize),bmin/(2.*pixsize),bpa)
	ax.contour(res2(X,Y),colors='w',levels=amp_cal[0]*np.array([0.2,0.5,0.8]))
	ax.coords['ra'].set_axislabel('Right Ascension')
	ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
	ax.coords['ra'].set_major_formatter('hh:mm:ss')
	ax.set_title('Calibrator Fit')
	plt.savefig(my_dir+'calibrator_fit.png')
	plt.show()
	raw_input('Please press enter when ready to continue.')

#fit target by fixing beam size
data_target=getdata(fitsim,0)[0,ranges[2]:ranges[3],ranges[0]:ranges[1]]
vartarget=getdata(fitsim,1)[0,ranges[2]:ranges[3],ranges[0]:ranges[1]]
print "Fitting Target..."
fix=raw_input('Do you want to fix the beam to the starlink values or the calibrator fit (if selected above)? star or cal?--> ')
if fix=='cal':
	param_guess=[flux_guess[1]/1e3,data_target.shape[0]/2.,data_target.shape[0]/2.,bmaj,bmin,bpa]
elif fix=='star':
	param_guess=[flux_guess[1]/1e3,data_target.shape[0]/2.,data_target.shape[0]/2.,w[0],w[1],w[2]]
tar_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':True,'y_stddev':True,'theta':True}
x=np.arange(0,len(data_target[0,:]))
y=np.arange(0,len(data_target[:,0]))
X, Y = np.meshgrid(x, y)
Z=data_target
if fit_alg=='lsq':
	M=models.Gaussian2D(param_guess[0],param_guess[1],param_guess[2],param_guess[3]/(2.*pixsize),\
		param_guess[4]/(2.*pixsize),param_guess[5],fixed=tar_fixed_params)
	lmf=fitting.LevMarLSQFitter()
	res=lmf(M,X,Y,Z)
	amp_tar=[res.parameters[0]*1e3,1e3*np.sqrt(np.diag(lmf.fit_info['param_cov']))[0]]
	xxt=res.parameters[1]
	yyt=res.parameters[2]
elif fit_alg=='mcmc':
	ampt,xxt0,yyt0,bmajt,bmint,bpat=mcmc_fit(data_target*1e3,vartarget*1e3,[param_guess[0]*1e3,param_guess[1],param_guess[2],w[0],\
		w[1],w[2]],[False,False,False,True,True,True],500,1500,'n',pixsize)
	amp_tar=[ampt[0],np.sqrt(ampt[1]**2+ampt[2]**2)]
	xxt=xxt0[0]
	yyt=yyt0[0]

#print results
print '\n'
print 'Target scan fit results (fixed to '+fix+'):'
print 'Target Flux = ','%.4f'%(amp_tar[0]),' +/- ','%.4f' %(amp_tar[1]), 'mJy','\n'
outfile.write('Target fit results:\n')
outfile.write('\n')
outfile.write('Target Flux = '+'%.4f'%(amp_tar[0])+' +/- '+'%.4f'%(amp_tar[1])+'mJy')
outfile.close()



#show results
print 'Plotting result...'
fig=plt.figure()
ax=plt.subplot(111,projection=wmaptar.celestial)
ax.imshow(data_target,origin='lower')
res2=models.Gaussian2D(amp_tar[0],xxt,yyt,param_guess[3]/(2.*pixsize),param_guess[4]/(2.*pixsize),param_guess[5])
ax.contour(res2(X,Y),colors='w',levels=amp_tar[0]*np.array([0.2,0.5,0.8]))
ax.coords['ra'].set_axislabel('Right Ascension')
ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax.coords['ra'].set_major_formatter('hh:mm:ss')
ax.set_title('Target Fit (beam size/angle fixed)')
plt.savefig(my_dir+'target_fit.png')
plt.show()

print 'MJD midpoint:'
hdu=fits.open(fitsim)[0]
mjderr=((Time(hdu.header['DATE-END'],format='isot',scale='utc').mjd-Time(hdu.header['DATE-OBS'],format='isot',scale='utc').mjd)/2.)
print Time(hdu.header['DATE-OBS'],format='isot',scale='utc').mjd+mjderr, "+/-", mjderr



print '\n'
print '*********************************************************'
print 'Script is finished. Please inspect output data products.'
print '*********************************************************'
