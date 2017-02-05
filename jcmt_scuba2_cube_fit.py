import pyfits
import numpy as np
import pylab as pl
import math as ma
import matplotlib.pyplot as plt
import gaussfitter as gf
from astropy.table import hstack
from photutils import CircularAperture,aperture_photometry,CircularAnnulus
from mpl_toolkits.mplot3d.axes3d import Axes3D

def read_data_cube(fits_file,flag):
	#basic info on fits file
	info_data=pyfits.info(fits_file)
	#header info
	header_primary0 = pyfits.getheader(fits_file)
	#get data from cube
	data_cube = pyfits.getdata(fits_file,0)
	if flag =='cube':
		var_cube=pyfits.getdata(fits_file,2)
		header_primary1 = pyfits.getheader(fits_file,1)
		header_primary2 = pyfits.getheader(fits_file,2)
		header_primary3 = pyfits.getheader(fits_file,3)
	elif flag == 'full':
		var_cube=pyfits.getdata(fits_file,1)
		header_primary1 = pyfits.getheader(fits_file,1)
	else:
		var_cube=0.
	#check you get a numpy nd array
	ty=type(data_cube)
	#check what dimensions are (time bin,y,x)
	dim=data_cube.shape
	return(data_cube,ty,dim,var_cube)

def aper_photo(data,X,Y,pix_rad,var):
	positions=(X,Y)
	apertures = CircularAperture(positions, r=pix_rad)
	annulus_apertures = CircularAnnulus(positions, r_in=4.0, r_out=5.0)
	rawflux_table0 = aperture_photometry(data, apertures,error=var,method='exact')
	bkgflux_table0 = aperture_photometry(data, annulus_apertures,error=var,method='exact')
	rawflux_table = aperture_photometry(data, apertures,method='exact')
	bkgflux_table = aperture_photometry(data, annulus_apertures,method='exact')
	phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
	bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
	bkg_sum = bkg_mean * apertures.area()
	final_sum = phot_table['aperture_sum_raw'] - bkg_sum
	phot_table['residual_aperture_sum'] = final_sum
	#print(phot_table['residual_aperture_sum'][0])
	#print rawflux_table['aperture_sum'][0]
	return(phot_table['residual_aperture_sum'][0],\
		np.sqrt((rawflux_table0['aperture_sum_err'][0])**2+(bkgflux_table0['aperture_sum_err'][0])**2))

def gauss_fit(data,var,param_guess,fi,fp):
	height,amplitude,X,Y,width_x,width_y,rota=param_guess[0],param_guess[1],param_guess[2],param_guess[3],\
	param_guess[4],param_guess[5],param_guess[6]
	if fp =='f':
		psf_fit,im_fit=gf.gaussfit(data,error=var,params=(height,amplitude,X,Y,width_x,width_y,rota),\
			returnfitimage=True,returnmp=True,fixed=fi,return_all=1,rotate=1,vheight=1)
	else:
		psf_fit,im_fit=gf.gaussfit(data,error=var,returnfitimage=True,returnmp=True,fixed=fi,return_all=1,\
			rotate=1,vheight=1)
	return(psf_fit.params,im_fit,psf_fit.perror,psf_fit.covar)
def image_display(img_data,img_dem,name):
	if img_dem==None:
		im0=plt.imshow(img_data,cmap='jet',origin="lower")
		cbar=plt.colorbar(im0, orientation='vertical')
		plt.grid(True,color='w', lw=3)
		plt.title(name)
		plt.show()
	else:
		img_lowx,img_highx,img_lowy,img_highy=img_dem[0],img_dem[1],img_dem[2],img_dem[3]
		#show results 2d image
		im0=plt.imshow(img_data,cmap='jet',origin="lower")
		cbar=plt.colorbar(im0, orientation='vertical')
		plt.grid(True,color='w', lw=3)
		plt.title(name)
		plt.show()
		#show results surface plot
		figs=plt.figure()
		axs=figs.add_subplot(1,1,1,projection='3d')
		ps=axs.plot_surface(np.arange(img_lowx,img_highx,1),np.arange(img_lowy,img_highy,1),\
			img_data[img_lowx:img_highx,img_lowy:img_highy],\
			rstride=1,cstride=1,cmap=plt.cm.get_cmap('coolwarm',50),linewidth=0,antialiased=False)
		axs.view_init(30,260)
		axs.tick_params(which='major',axis="y", labelsize=16, length=6, width=2, pad=5)
		axs.tick_params(which='major',axis="x", labelsize=16, length=6, width=2, pad=5)
		axs.tick_params(which='major',axis="z", labelsize=16, length=6, width=2, pad=5)
		cb=plt.colorbar(ps,shrink=1.0)
		cb.ax.tick_params(labelsize=12)
		plt.xlabel('X',fontsize=16)
		plt.ylabel('Y',fontsize=16)
		plt.show()

print 'Reading in parameters...'
####################################
#User input
####################################
my_dir='/Users/atetarenk/Desktop/V404_mm_radio_timing_res/jcmt_timing/'
band=850
date='jun22'
do_timing='y'
#########
scan_list8=[36,37,41,42,46,47,52]#850
flux_ratio8=[2.359/630.371,2.359/630.371,2.359/630.371,2.359/630.371,2.347/645.264,2.347/645.264,2.347/645.264]#850
ranges8=[[123,133,129,139],[122,132,128,138],[110,125,118,135],[120,130,125,135],[110,135,130,155],\
[93,123,116,146],[107,123,121,137]]#850
w850=1.63
pixrad8=1.2
psf_dim8=[0,48,0,48]
psf_fit_guess8=[0.0,0.8,24,24,1.0,1.0,0.0]#height,amplitude,X,Y,width_x,width_y,rota
psf_string8=my_dir+'out_crl2688_32_cal_psf.fits'#or54
full_fit_guess8=[0.0,1.0,1.0,1.0,1.63,1.63,0.0]#height,amplitude,X,Y,width_x,width_y,rota
#########
scan_list4=[36,37,41,42,46]#450
flux_ratio4=[7.873/641.625,7.873/641.625,7.873/641.625,7.873/641.625,9.701/638.043]#450
ranges4=[[252,262,267,275],[250,260,266,276],[240,250,260,270],[240,250,245,255],[243,253,278,288]]#450
w450=2.17
pixrad4=1.4
psf_dim4=[0,50,0,50]
psf_fit_guess4=[0.0,0.8,25,25,1.0,1.0,0.0]#height,amplitude,X,Y,width_x,width_y,rota
psf_string4=my_dir+'out_crl2688_32_cal4_psf.fits'#or54
full_fit_guess4=[0.0,1.0,1.0,1.0,2.17,2.17,0.0]#height,amplitude,X,Y,width_x,width_y,rota
#########
psf_fixed_params=np.array([True, False, False, False, False, False, True], dtype=bool)
full_fixed_params=np.array([True, False, False, False, False,False, True], dtype=bool)
cube_fixed_params=np.array([True, False, True,True,True,True,True], dtype=bool)
####################################
if band==850:
	scan_list=scan_list8
elif band==450:
	scan_list=scan_list4

chi2=[]
image_fit=[]
amp=[]
err=[]
param=[]
param_err=[]
amp_p=[]
amp_pe=[]

print 'There are ', len(scan_list), ' target scans present.'

#fit psf from full data set
if band==850:
	data_psf,type_data_psf,shape_data_psf,var_psf=read_data_cube(psf_string8,'p')
	print "Fitting calibator scan..."
	psf_params,img_psf,psf_err,psf_cov=gauss_fit(data_psf[psf_dim8[0]:psf_dim8[1],psf_dim8[2]:psf_dim8[3]],None,\
		psf_fit_guess8,psf_fixed_params,'f')
	#print results
	print 'Calibrator scan fit results:'
	print 'h,amp,x,y,wx,wy,rot'
	print 'err(h,amp,x,y,wx,wy,rot)'
	print psf_params
	print psf_err
	#show results
	image_display(img_psf,psf_dim8,'fit_psf')
	raw_input('Please press enter when ready to continue.')
elif band==450:
	data_psf,type_data_psf,shape_data_psf,var_psf=read_data_cube(psf_string4,'p')
	print "Fitting calibator scan..."
	psf_params,img_psf,psf_err,psf_cov=gauss_fit(data_psf[psf_dim4[0]:psf_dim4[1],psf_dim4[2]:psf_dim4[3]],None,\
		psf_fit_guess4,psf_fixed_params,'f')
	#print results
	print 'Calibrator scan fit results:'
	print 'h,amp,x,y,wx,wy,rot'
	print 'err(h,amp,x,y,wx,wy,rot)'
	print psf_params
	print psf_err
	#show results
	image_display(img_psf,psf_dim4,'fit_psf')
	raw_input('Please press enter when ready to continue.')

for j in range(0,len(scan_list)):
	#fits file names
	if band==850:
		shortmap_string=my_dir+'v404cyg_'+str(scan_list[j])+'_shortmaps_cube_cal.fits'
		full_string=my_dir+'v404cyg_'+str(scan_list[j])+'_fullmap_cal.fits'
		#read in data cube, full map psf and full map
		print 'Reading in fits files...'
		data,type_data,shape_data,var=read_data_cube(shortmap_string,'cube')
		data_full,type_data_full,shape_data_full,var_full=read_data_cube(full_string,'full')
		#fit full data set
		print 'Fitting full target scan mosaic...'
		full_params,img_full,full_err,full_cov=gauss_fit(data_full[0,ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]],\
			var_full[0,ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]],full_fit_guess8,\
			full_fixed_params,'t')
	elif band==450:
		shortmap_string=my_dir+'v404cyg_'+str(scan_list[j])+'_shortmaps_cube_cal4f.fits'
		full_string=my_dir+'v404cyg_'+str(scan_list[j])+'_fullmap_cal4.fits'
		#read in data cube, full map psf and full map
		print 'Reading in fits files...'
		data,type_data,shape_data,var=read_data_cube(shortmap_string,'cube')
		data_full,type_data_full,shape_data_full,var_full=read_data_cube(full_string,'full')
		#fit full data set
		print 'Fitting full target scan mosaic...'
		full_params,img_full,full_err,full_cov=gauss_fit(data_full[0,ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]],\
			var_full[0,ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]],full_fit_guess4,\
			full_fixed_params,'t')
	#print results
	print 'Full scan ',scan_list[j], 'fit results:'
	print 'h,amp,x,y,wx,wy,rot'
	print 'err(h,amp,x,y,wx,wy,rot)'
	print full_params
	print full_err
	#show results
	image_display(img_full,None,'fit_full_scan'+str(scan_list[j]))
	raw_input('Please press enter when ready to continue.')
	#aperature photometry version
	#ap_phot,er_phot=aper_photo(data_full[0,ranges[j][0]:ranges[j][1],ranges[j][2]:ranges[j][3]],\
		#full_params[2],full_params[3],2,np.sqrt(var_full[0,ranges[j][0]:ranges[j][1],ranges[j][2]:ranges[j][3]]))

	#fit all planes of cube
	if do_timing=='y':
		print 'Fitting cube for scan ',scan_list[j],' ...'
		for kk in range(0,shape_data[0]):
			if band==850:
				cube_params,img_cube,cube_err,cube_cov=gauss_fit(data[kk,ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]],\
					np.sqrt(var[kk,ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]]),\
					[0.0,1.0,full_params[2],full_params[3],w850,w850,0.0],\
					cube_fixed_params,'f')
				chi2.append(((data[kk,ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]]-img_cube)**2)/(var[kk,\
					ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]]).sum())
			elif band==450:
				cube_params,img_cube,cube_err,cube_cov=gauss_fit(data[kk,ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]],\
					np.sqrt(var[kk,ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]]),\
					[0.0,1.0,full_params[2],full_params[3],w450,w450,0.0],\
					cube_fixed_params,'f')
				chi2.append(((data[kk,ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]]-img_cube)**2)/(var[kk,\
					ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]]).sum())
			#aperture photometry version
			if band==850:
				ap_params,ap_err=aper_photo((flux_ratio8[j])*data[kk,ranges8[j][0]:ranges8[j][1],ranges8[j][2]:ranges8[j][3]],\
					full_params[2],full_params[3],pixrad8,np.sqrt(((flux_ratio8[j])**2)*var[kk,ranges8[j][0]:ranges8[j][1],\
						ranges8[j][2]:ranges8[j][3]]))
				amp_p.append(ap_params*np.pi*((4.*pixrad8)**2))
				amp_pe.append(ap_err*np.pi*((4.*pixrad8)**2))
			elif band==450:
				ap_params,ap_err=aper_photo((flux_ratio4[j])*data[kk,ranges4[j][0]:ranges4[j][1],ranges4[j][2]:ranges4[j][3]],\
					full_params[2],full_params[3],pixrad4,np.sqrt(((flux_ratio4[j])**2)*var[kk,ranges4[j][0]:ranges4[j][1],\
						ranges4[j][2]:ranges4[j][3]]))
				amp_p.append(ap_params*np.pi*((2.*pixrad4)**2))
				amp_pe.append(ap_err*np.pi*((2.*pixrad4)**2))
			image_fit.append(img_cube)
			amp.append(cube_params[1])
			err.append(cube_err[1])
			param.append(cube_params)
			param_err.append(cube_err)

if do_timing=='y':
	print 'Reading in MJD values...'
	if band==850:
		d=np.loadtxt(my_dir+'JCMT_850_jun22.txt')
		mjd=d[:,0]
		print 'len test', len(mjd),len(amp)
	elif band==450:
		d=np.loadtxt(my_dir+'JCMT_450_jun22.txt')
		mjd=d[:,0]
		print 'len test', len(mjd),len(amp)

	print "Plotting results..."
	#plot gaussan fit and aperature photometry results
	fig=plt.figure()
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	plt.errorbar(mjd,np.array(amp),np.array(err),fmt='mo')
	plt.errorbar(mjd,np.array(amp_p),np.array(amp_pe),fmt='co')
	plt.show()
	raw_input('Please press enter when ready to continue.')

	print 'Writing results to file...'
	#write results to file
	file_res_g=open(my_dir+'JCMT_'+str(band)+'_'+date+'_gfit_redo17.txt','w')
	for jj in range(0,len(mjd)):
		if amp[jj] != np.nan:
			file_res_g.write(' {0} {1} {2}\n'.format(mjd[jj],amp[jj]*1000.,err[jj]*1000.))
	file_res_g.close()
	file_res_a=open(my_dir+'JCMT_'+str(band)+'_'+date+'_phot_redo17.txt','w')
	for ll in range(0,len(mjd)):
		if amp_p[ll] != np.nan:
			file_res_a.write(' {0} {1} {2}\n'.format(mjd[ll],amp_p[ll]*1000.,amp_pe[ll]*1000.))
	file_res_a.close()

	print 'Showing final lightcurve written to file...'
	#show final lightcurve that was written to file
	data=np.loadtxt(my_dir+'JCMT_'+str(band)+'_'+date+'_gfit_redo17.txt')
	datap=np.loadtxt(my_dir+'JCMT_'+str(band)+'_'+date+'_phot_redo17.txt')

	plt.figure()
	plt.errorbar(data[:,0],data[:,1],data[:,2],linestyle='',marker='o',color='m',label='Gauss')
	plt.errorbar(datap[:,0],datap[:,1],datap[:,2],linestyle='',marker='o',color='b',label='Phot')
	plt.ylim(-1000,6000)
	plt.legend(loc='upper left',numpoints=1)
	plt.show()

print 'Script is finished.'