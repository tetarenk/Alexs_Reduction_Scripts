#!/bin/bash

#JCMT SCUBA-2 reduction script for bright compact sources
#INPUT: parameters defined in section below, raw scuba2 data, map cropping param file (mypars)
#OUTPUT: (1) Maps of calibrted 850um and 450um scan(s) -- [target]_[date]_[scan]_[band]_fullmap_cal_crop(_mf).sdf/.fits
#        (2) Maps of 850um and 450um calibrator scan -- [cal]_[date]_[cal_scan]_[band]_fullmap_cal_crop(_mf).sdf
#        (3) Noise Maps -- [target]_[date]_[scan]_[band]_fullmap_cal_crop(_mf)_noi.sdf
#        (4) Timing Cubes (if requested) -- [target]_[date]_[scan]_[band]_shortmap_cube_cal.sdf/.fits
#        (5) Combined scans map (if scans >1) -- [target]_[date]_[scan]_[band]_fullmap_all.sdf/.fits
#        (6) Combined scans timing cube (if scans >1) -- [target]_[date]_[band]_shortmap_cube_cal_all.sdf/.fits
#        (7) Log file of FCFs, beam sizes, and maps RMS for individual scans -- output_results.logg
#Written by: Alex J Tetarenko
#Last Updated: April 26, 2019

############################
#Defining variables section
############################
data_dir=/path/to/raw
my_dir=/path/to/results
#file_lst=/export/data2/atetarenko/JCMT_maxi1820/cadcUrlList5.txt
date=20190422
#scans need to be 2 digits, e.g., scan 7 is 07
scan8=(64 65 69 70 75 76 79 80)
scan4=()
numscans=${#scan8[@]}
cal_scan4=72
cal_scan8=72
cal='crl2688'
target='cygx3'
do_timing="y"
do_mb='n'
shortm8=200
shortm4=400
mypars=/path/to/mypars.lis
############################


#download data
#echo 'Downloading data...'
#wget --http-user=[cadc_user] --http-password=[cadc_password] --content-disposition -i $file_lst

echo '###############################################'
echo 'Welcome to Alexs JCMT SCUBA-2 Reduction Script'
echo '###############################################'
echo ''

#start and intitalize starlink modules
echo 'Initializing Starlink software...'
source $STARLINK_DIR/etc/profile
eval 'kappa'
eval 'smurf'
eval 'convert'
export ORAC_DATA_OUT=$my_dir

#copy mapping recipe-only do once, add shortmap parameter to timing config files.
FILE1=$my_dir/dimmconfig'_'bright'_'compact.lis
FILE2=$my_dir/dimmconfig'_'bright'_'compact_shortmap8.lis
FILE3=$my_dir/dimmconfig'_'bright'_'compact_shortmap4.lis
if [ -f $FILE1 ]
then
	echo 'Point source mapping recipe already exists.'
else
	echo 'Copying point source mapping recipe to current directory for use...'
	cp -r $STARLINK_DIR/share/smurf/dimmconfig_bright_compact.lis $my_dir/dimmconfig'_'bright'_'compact.lis
fi
if [ $do_timing = "y" ];
then
	if [ -f $FILE2 ];
	then
		echo 'timing mapping recipe at 850um already exists.'
	else
		echo 'Copying 850um timing mapping recipe to current directory for use...'
		cp -r $STARLINK_DIR/share/smurf/dimmconfig_bright_compact.lis $my_dir/dimmconfig'_'bright'_'compact_shortmap8.lis
		echo "   shortmap=$shortm8" >> $my_dir/dimmconfig'_'bright'_'compact_shortmap8.lis
	fi
	if [ ${#scan4[@]} -gt 0 ]
	then
		if [ -f $FILE3 ];
		then
			echo 'timing mapping recipe at 450um already exists.'
		else
			echo 'Copying 450um timing mapping recipe to current directory for use...'
			cp -r $STARLINK_DIR/share/smurf/dimmconfig_bright_compact.lis $my_dir/dimmconfig'_'bright'_'compact_shortmap4.lis
			echo "   shortmap=$shortm4" >> $my_dir/dimmconfig'_'bright'_'compact_shortmap4.lis
		fi
	fi
fi
read -r -p 'Please press enter when ready to continue >>> '
Nscans=$(expr $numscans - 1)
for i in $(seq 0 $Nscans)
do
	scan8=${scan8[$i]}
	if [ ${#scan4[@]} -eq 0 ]
	then
		echo '****************************************************'
		echo 'Reducing 850um scan '$scan8
		echo '****************************************************'
	else
		scan4=${scan4[$i]}
		echo '****************************************************'
		echo 'Reducing 850um scan '$scan8' and 450um scan '$scan4
		echo '****************************************************'
	fi
	read -r -p 'Please press enter when ready to continue >>> '
	#create file first (one triangle is not exist yet, two is append)
	echo 'Create list of files for target scan...'
	rm -rf $my_dir/$target'_'$date'_'850'_'$scan8.lst
	ls $data_dir/s8*$date'_'000$scan8*.sdf > $my_dir/$target'_'$date'_'850'_'$scan8.lst
	#read -r -p 'Please press enter when ready to continue >>> '
	#make full scan maps
	rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap.sdf'
	echo 'Making 850um map for full target scan...'
	makemap in=^$my_dir/$target'_'$date'_'850'_'$scan8.lst out=$my_dir/$target'_'$date'_'$scan8'_850_fullmap' config=^$my_dir/dimmconfig'_'bright'_'compact.lis
	#read -r -p 'Please press enter when ready to continue >>> '
	if [ ${#scan4[@]} -gt 0 ]
	then
		rm -rf $my_dir/$target'_'$date'_'450'_'$scan4.lst
		ls $data_dir/s4*$date'_'000$scan4*.sdf > $my_dir/$target'_'$date'_'450'_'$scan4.lst
		rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap.sdf'
		echo 'Making 450um map for full target scan...'
		makemap in=^$my_dir/$target'_'$date'_'450'_'$scan4.lst out=$my_dir/$target'_'$date'_'$scan4'_450_fullmap' config=^$my_dir/dimmconfig'_'bright'_'compact.lis
		#read -r -p 'Please press enter when ready to continue >>> '
	fi


	#calibrate to Jy
	#reduce cal scan
	echo 'Reducing calibrator scans...'
	echo 'Making list of files for calibrator scan...'
	FILE_CAL8=$my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap.sdf'
	
	if [ -f $FILE_CAL8 ]
	then
		echo '850um calibrator scan already reduced.'
	else
		ls $data_dir/s8*$date'_'000$cal_scan8*.sdf > $my_dir/$cal'_'$date'_'850'_'$cal_scan8.lst 
		echo 'Making 850um map for cal scan...'
		makemap in=^$my_dir/$cal'_'$date'_'850'_'$cal_scan8.lst out=$my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap' config=^$my_dir/dimmconfig'_'bright'_'compact.lis
		#read -r -p 'Please press enter when ready to continue >>> '
		# run picard recipe,sf prints to screen and log file (-log sf)
		echo 'Running SCUBA2_CHECK_CAL recipe for 850um cal scan...'
		picard -log sf SCUBA2_CHECK_CAL $my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap.sdf'
		read -r -p 'Please enter fcf for 850um >>> ' fcf8
		read -r -p 'Please enter BMAJ (arcsec) for 850um >>> ' bmaj8
		read -r -p 'Please enter BMIN (arcsec) for 850um >>> ' bmin8
		read -r -p 'Please enter BPA (deg) for 850um >>> ' bpa8
		echo 'FCF850_'$scan8'='$fcf8 > $my_dir/output_results.logg
		echo 'BMAJ850_'$scan8'='$bmaj8 >> $my_dir/output_results.logg
		echo 'BMIN850_'$scan8'='$bmin8 >> $my_dir/output_results.logg
		echo 'BPA850_'$scan8'='$bpa8 >> $my_dir/output_results.logg
		cmult in=$my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap.sdf' out=$my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap_cal.sdf' scalar=$fcf8
		ndf2fits $my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap_cal.sdf' $my_dir/$cal'_'$date'_'$cal_scan8'_850_fullmap_cal.fits'
		read -r -p 'Please press enter when ready to continue >>> '
		## Storing new FCF for 850: 549.081 +/-   9.185 Jy/beam/pW (cf 537: 2.2% higher)
	fi
	if [ ${#scan4[@]} -gt 0 ]
	then
		FILE_CAL4=$my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap.sdf'
		if [ -f $FILE_CAL4 ]
		then
			echo '450um calibrator scan already reduced.'
		else
			ls $data_dir/s4*$date'_'000$cal_scan4*.sdf > $my_dir/$cal'_'$date'_'450'_'$cal_scan4.lst
			echo 'Making 450um map for cal scan...'
			makemap in=^$my_dir/$cal'_'$date'_'450'_'$cal_scan4.lst out=$my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap' config=^$my_dir/dimmconfig'_'bright'_'compact.lis
			#read -r -p 'Please press enter when ready to continue >>> '
			# run picard recipe,sf prints to screen and log file (-log sf)
			echo 'Running SCUBA2_CHECK_CAL recipe for 450um cal scan...'
			picard -log sf SCUBA2_CHECK_CAL $my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap.sdf'
			## Storing new FCF for 850: 549.081 +/-   9.185 Jy/beam/pW (cf 537: 2.2% higher)
			read -r -p 'Please enter fcf for 450um >>> ' fcf4
			read -r -p 'Please enter BMAJ (arcsec) for 450um >>> ' bmaj4
			read -r -p 'Please enter BMIN (arcsec) for 450um >>> ' bmin4
			read -r -p 'Please enter BPA (deg) for 450um >>> ' bpa4
			echo 'FCF450_'$scan4'='$fcf4 >> $my_dir/output_results.logg
			echo 'BMAJ450_'$scan4'='$bmaj4 >> $my_dir/output_results.logg
			echo 'BMIN450_'$scan4'='$bmin4 >> $my_dir/output_results.logg
			echo 'BPA450_'$scan4'='$bpa4 >> $my_dir/output_results.logg
			cmult in=$my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap.sdf' out=$my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap_cal.sdf' scalar=$fcf4
			ndf2fits $my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap_cal.sdf' $my_dir/$cal'_'$date'_'$cal_scan4'_450_fullmap_cal.fits'
			read -r -p 'Please press enter when ready to continue >>> '
		fi
	fi
	echo 'Running cmult to do absolute flux calibration on target maps...'
	if [ ${#scan4[@]} -gt 0 ]
	then
		rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal.sdf'
		cmult in=$my_dir/$target'_'$date'_'$scan4'_450_fullmap.sdf' out=$my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal.sdf' scalar=$fcf4
	fi
	rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal.sdf'
	cmult in=$my_dir/$target'_'$date'_'$scan8'_850_fullmap.sdf' out=$my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal.sdf' scalar=$fcf8
		
	#crop maps
	echo 'Cropping target maps...'
	if [ ${#scan4[@]} -gt 0 ]
	then
		rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf'
		picard -log sf -recpars $mypars CROP_SCUBA2_IMAGES $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal.sdf'
		#mv $target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf' $my_dir
		#mv $target'_'$date'_'$scan4'_450_fullmap_cal_crop_psf.sdf' $my_dir
	fi
	rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf'
	picard -log sf -recpars $mypars CROP_SCUBA2_IMAGES $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal.sdf'
	#mv $target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf' $my_dir
	#mv $target'_'$date'_'$scan8'_850_fullmap_cal_crop_psf.sdf' $my_dir
	

	#matched beam filter
	if [ $do_mb = 'y' ]
	then
		echo 'Applying matched beam filter to better find point sources in target maps...'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf'
		picard -log sf SCUBA2_MATCHED_FILTER $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf'
		mv $target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf' $my_dir
		if [ ${#scan4[@]} -gt 0 ]
		then
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf'
			picard -log sf SCUBA2_MATCHED_FILTER $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf'
			mv $target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf' $my_dir
		fi
	fi

	#make noise map, get noise estimate of map
	if [ $do_mb = 'y' ]
	then
		echo 'Making noise map and esimating noise in target map at 850um...'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf_snr.sdf'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf_noi.sdf'
		
		
		makesnr in=$my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf' out=$my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf_snr.sdf' minvar=1.0e-20
		stats comp=err $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf'
		read -r -p 'Please enter RMS (mean) for 850um >>> ' rms8
		echo 'RMS850_'$scan4'='$rms8 >> $my_dir/output_results.logg
		read -r -p 'Please press enter when ready to continue >>> '
		ndfcopy $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf' comp=err $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf_noi.sdf'
		if [ ${#scan4[@]} -gt 0 ]
		then
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf_snr.sdf'
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf_noi.sdf'
			echo 'Making noise map and esimating noise in target map at 450um...'
			makesnr in=$my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf' out=$my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf_snr.sdf' minvar=1.0e-20
			stats comp=err $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf'
			read -r -p 'Please enter RMS (mean) for 450um >>> ' rms4
			echo 'RMS450_'$scan4'='$rms4 >> $my_dir/output_results.logg
			read -r -p 'Please press enter when ready to continue >>> '
			ndfcopy $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf' comp=err $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf_noi.sdf'
		fi
	else
		echo 'Making noise map and esimating noise in target map at 850um...'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_snr.sdf'
		rm -rf $my_dir/$targe'_'$datet'_'$scan8'_850_fullmap_cal_crop_noi.sdf'
		rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_snr.sdf'
		rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_noi.sdf'
		makesnr in=$my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf' out=$my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_snr.sdf' minvar=1.0e-20
		stats comp=err $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf'
		read -r -p 'Please enter RMS (mean) for 850um >>> ' rms8
		echo 'RMS850_'$scan8'='$rms8 >> $my_dir/output_results.logg
		read -r -p 'Please press enter when ready to continue >>> '
		ndfcopy $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf' comp=err $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_noi.sdf'
		if [ ${#scan4[@]} -gt 0 ]
		then
			echo 'Making noise map and esimating noise in target map at 450um...'
			makesnr in=$my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf' out=$my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_snr.sdf' minvar=1.0e-20
			stats comp=err $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf'
			read -r -p 'Please enter RMS (mean) for 450um >>> ' rms4
			echo 'RMS450_'$scan4'='$rms4 >> $my_dir/output_results.logg
			read -r -p 'Please press enter when ready to continue >>> '
			ndfcopy $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf' comp=err $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_noi.sdf'
		fi
	fi

	#view maps
	if [ $do_mb = 'y' ]
	then
		echo 'Viewing 850um target map...'
		gaia $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf'
		read -r -p 'Please press enter when ready to continue >>> '
		echo 'Viewing 850um noise map...'
		gaia $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf_noi.sdf'
		read -r -p 'Please press enter when ready to continue >>> '
		if [ ${#scan4[@]} -gt 0 ]
		then
			echo 'Viewing 450um target map...'
			gaia $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf'
			read -r -p 'Please press enter when ready to continue >>> '
			echo 'Viewing 450um noise map...'
			gaia $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf_noi.sdf'
			read -r -p 'Please press enter when ready to continue >>> '
		fi
	else
		echo 'Viewing 850um target map...'
		gaia $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf'
		read -r -p 'Please press enter when ready to continue >>> '
		echo 'Viewing 850um noise map...'
		gaia $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_noi.sdf'
		read -r -p 'Please press enter when ready to continue >>> '
		if [ ${#scan4[@]} -gt 0 ]
		then
			echo 'Viewing 450um target map...'
			gaia $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf'
			read -r -p 'Please press enter when ready to continue >>> '
			echo 'Viewing 450um noise map...'
			gaia $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_noi.sdf'
			read -r -p 'Please press enter when ready to continue >>> '
		fi
	fi

	#convert to fits
	echo 'Converting target maps to fits...'
	if [ $do_mb = 'y' ]
	then
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.fits'
		ndf2fits $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.sdf' $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop_mf.fits'
		if [ ${#scan4[@]} -gt 0 ]
		then
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.fits'
			ndf2fits $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.sdf' $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop_mf.fits'
		fi
	else
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.fits'
		ndf2fits $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.sdf' $my_dir/$target'_'$date'_'$scan8'_850_fullmap_cal_crop.fits'
		if [ ${#scan4[@]} -gt 0 ]
		then
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.fits'
			ndf2fits $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.sdf' $my_dir/$target'_'$date'_'$scan4'_450_fullmap_cal_crop.fits'
		fi
	fi

	#high time res
	if [ $do_timing = 'y' ]
	then
		echo 'Starting timing analysis on scans...'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_shortmap.sdf'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube.sdf'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube_cal.sdf'
		rm -rf $my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube_cal.fits'
		echo 'Making 850 shortmaps cube...'
		makemap in=^$my_dir/$target'_'$date'_'850'_'$scan8.lst out=$my_dir/$target'_'$date'_'$scan8'_850_shortmap' config=^$my_dir/dimmconfig'_'bright'_'compact_shortmap8.lis
		#read -r -p 'Please press enter when ready to continue >>> '
		echo 'Stacking all 850um shortmaps...'
		stackframes $my_dir/$target'_'$date'_'$scan8'_850_shortmap.more.smurf.shortmaps' sort=true sortby=MJD-AVG $my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube.sdf'
		#read -r -p 'Please press enter when ready to continue >>> '
		FILE_CALTIM8=$my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap.sdf'
		if [ -f $FILE_CALTIM8 ]
		then
			echo '850 calibrator shortmaps already created.'
		else
			echo 'Making 850um shortmaps for cal scan...'
			makemap in=^$my_dir/$cal'_'$date'_'850'_'$cal_scan8.lst out=$my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap' config=^$my_dir/dimmconfig'_'bright'_'compact_shortmap8.lis
			#read -r -p 'Please press enter when ready to continue >>> '
			echo 'Stacking all 850um cal shortmaps...'
			stackframes $my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap.more.smurf.shortmaps' sort=true sortby=MJD-AVG $my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap_cube.sdf'
			#echo 'Running SCUBA2_CHECK_CAL recipe for 850um cal shortmap cube...'
			#picard -log sf SCUBA2_CHECK_CAL $my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap_cube.sdf'
			fcf8short=$fcf8
			cmult in=$my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap_cube.sdf' out=$my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap_cube_cal.sdf' scalar=$fcf8short
			ndf2fits $my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap_cube_cal.sdf' $my_dir/$cal'_'$date'_'$cal_scan8'_850_shortmap_cube_cal.fits'
			#read -r -p 'Please enter fcf for 850um >>> ' fcf8short
			read -r -p 'Please press enter when ready to continue >>> '
		fi
		echo 'Running cmult to do absolute flux calibration on target shortmaps cube...'
		cmult in=$my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube.sdf' out=$my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube_cal.sdf' scalar=$fcf8short
		#read -r -p 'Please press enter when ready to continue >>> '
		if [ ${#scan4[@]} -gt 0 ]
		then
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_shortmap.sdf'
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube.sdf'
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube_cal.sdf'
			rm -rf $my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube_cal.fits'
			echo 'Making 450 shortmaps cube...'
			makemap in=^$my_dir/$target'_'$date'_'450'_'$scan4.lst out=$my_dir/$target'_'$date'_'$scan4'_450_shortmap' config=^$my_dir/dimmconfig'_'bright'_'compact_shortmap4.lis
			#read -r -p 'Please press enter when ready to continue >>> '
			echo 'Stacking all 450um shortmaps...'
			stackframes $my_dir/$target'_'$date'_'$scan4'_450_shortmap.more.smurf.shortmaps' sort=true sortby=MJD-AVG $my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube.sdf'
			#read -r -p 'Please press enter when ready to continue >>> '
			FILE_CALTIM4=$my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap.sdf'
			if [ -f $FILE_CALTIM4 ]
			then
				echo '450 calibrator shortmaps already created.'
			else
				echo 'Making 450um shortmaps for cal scan...'
				makemap in=^$my_dir/$cal'_'$date'_'450'_'$cal_scan4.lst out=$my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap' config=^$my_dir/dimmconfig'_'bright'_'compact_shortmap4.lis
				#read -r -p 'Please press enter when ready to continue >>> '
				echo 'Stacking all 450um cal shortmaps...'
				stackframes $my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap.more.smurf.shortmaps' sort=true sortby=MJD-AVG $my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap_cube.sdf'
				#echo 'Running SCUBA2_CHECK_CAL recipe for 450um cal shortmap cube...'
				#picard -log sf SCUBA2_CHECK_CAL $my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap_cube.sdf'
				#read -r -p 'Please enter fcf for 450um >>> ' fcf4short
				fcf4short=$fcf4
				cmult in=$my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap_cube.sdf' out=$my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap_cube_cal.sdf' scalar=$fcf4short
				ndf2fits $my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap_cube_cal.sdf' $my_dir/$cal'_'$date'_'$cal_scan4'_450_shortmap_cube_cal.fits'
				read -r -p 'Please press enter when ready to continue >>> '
			fi
			echo 'Running cmult to do absolute flux calibration on target shortmaps cube...'
			cmult in=$my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube.sdf' out=$my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube_cal.sdf' scalar=$fcf4short
		fi
		echo 'Converting cubes to fits...'
		ndf2fits $my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube_cal.sdf' $my_dir/$target'_'$date'_'$scan8'_850_shortmap_cube_cal.fits'
		if [ ${#scan4[@]} -gt 0 ]
		then
			ndf2fits $my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube_cal.sdf' $my_dir/$target'_'$date'_'$scan4'_450_shortmap_cube_cal.fits'
		fi
	fi
done

#to combine scans if needed
rm -rf $my_dir/all_scan8.lst
if [ $do_mb = 'y' ]
then
	ls $my_dir/$target'_'$date'_'*'_850_fullmap_cal_crop_mf.sdf' > $my_dir/all_scan8.lst	
else
	ls $my_dir/$target'_'$date'_'*'_850_fullmap_cal_crop.sdf' > $my_dir/all_scan8.lst	
fi
NUMOFLINESscan8=$(cat $my_dir/all_scan8.lst | wc -l )
if [ $NUMOFLINESscan8 -gt 1 ]
then
	echo 'Mosaicing all target 850 um scans...'
	rm -rf $my_dir/$target'_'$date'_all_850_fullmap_all.sdf'
	rm -rf $my_dir/$target'_'$date'_all_850_fullmap_all.fits'
	wcsmosaic in=^$my_dir/all_scan8.lst out=$my_dir/$target'_'$date'_all_850_fullmap_all.sdf'
	ndf2fits $my_dir/$target'_'$date'_all_850_fullmap_all.sdf' $my_dir/$target'_'$date'_all_850_fullmap_all.fits'
	makesnr in=$my_dir/$target'_'$date'_all_850_fullmap_all.sdf' out=$my_dir/$target'_'$date'_all_850_fullmap_all_snr.sdf' minvar=1.0e-20
	ndfcopy $my_dir/$target'_'$date'_all_850_fullmap_all.sdf' comp=err $my_dir/$target'_'$date'_all_850_fullmap_all_noi.sdf'
	stats comp=err $my_dir/$target'_'$date'_all_850_fullmap_all.sdf'
	read -r -p 'Please enter RMS (mean) for 850um >>> ' rmsfull8
	echo 'RMSFULL850='$rmsfull8 >> $my_dir/output_results.logg
else
	echo 'Only 1 850um scan present.'
fi

if [ ${#scan4[@]} -gt 0 ]
then
	rm -rf $my_dir/all_scan4.lst
	if [ $do_mb = 'y' ]
	then
		ls $my_dir/$target'_'$date'_'*'_450_fullmap_cal_crop_mf.sdf' > $my_dir/all_scan4.lst
	else
		ls $my_dir/$target'_'$date'_'*'_450_fullmap_cal_crop.sdf' > $my_dir/all_scan4.lst
	fi
	NUMOFLINESscan4=$(cat $my_dir/all_scan4.lst | wc -l )
	if [ $NUMOFLINESscan4 -gt 1 ]
	then
		echo 'Mosaicing all target 450 um scans...'
		rm -rf $my_dir/$target'_'$date'_all_450_fullmap_all.sdf'
		rm -rf $my_dir/$target'_'$date'_all_450_fullmap_all.fits'
		wcsmosaic in=^$my_dir/all_scan4.lst out=$my_dir/$target'_'$date'_all_450_fullmap_all.sdf'
		ndf2fits $my_dir/$target'_'$date'_all_450_fullmap_all.sdf' $my_dir/$target'_'$date'_all_450_fullmap_all.fits'
		makesnr in=$my_dir/$target'_'$date'_all_450_fullmap_all.sdf' out=$my_dir/$target'_'$date'_all_450_fullmap_all_snr.sdf' minvar=1.0e-20
		ndfcopy $my_dir/$target'_'$date'_all_450_fullmap_all.sdf' comp=err $my_dir/$target'_'$date'_all_450_fullmap_all_noi.sdf'
		stats comp=err $my_dir/$target'_'$date'_all_450_fullmap_all.sdf'
		read -r -p 'Please enter RMS (mean) for 450um >>> ' rmsfull4
		echo 'RMSFULL450='$rmsfull4 >> $my_dir/output_results.logg
	else
		echo 'Only 1 450um scan present.'
	fi
fi



#to combine shortmaps cube scans if needed

if [ $do_timing = 'y' ]
then
	rm -rf $my_dir/all8.lst
	ls $my_dir/$target'_'$date'_'*'_850_shortmap_cube_cal.sdf' > $my_dir/all8.lst
	NUMOFLINES8=$(cat $my_dir/all8.lst | wc -l )
	if [ $NUMOFLINES8 -gt 1 ]
	then
		echo 'Mosaicing all 850um shortmap cubes for all target scans...'
		rm -rf $my_dir/$target'_'$date'_850_shortmap_cube_cal_all.sdf'
		rm -rf $my_dir/$target'_'$date'_850_shortmap_cube_cal_all.fits'
		wcsmosaic in=^$my_dir/all8.lst out=$my_dir/$target'_'$date'_850_shortmap_cube_cal_all.sdf'
		ndf2fits $my_dir/$target'_'$date'_850_shortmap_cube_cal_all.sdf' $my_dir/$target'_'$date'_850_shortmap_cube_cal_all.fits'
	else
		echo 'Only 1 850um scan present.'	
	fi
fi
if [ ${#scan4[@]} -gt 0 ]
then
	if [ $do_timing = 'y' ]
	then
		rm -rf $my_dir/all4.lst
		ls $my_dir/$target'_'$date'_'*'_450_shortmap_cube_cal.sdf' > $my_dir/all4.lst
		NUMOFLINES4=$(cat $my_dir/all4.lst | wc -l )
		if [ $NUMOFLINES4 -gt 1 ]
		then
			echo 'Mosaicing all 450um shortmap cubes for all target scans...'
			rm -rf $my_dir/$target'_'$date'_450_shortmap_cube_cal_all.sdf'
			rm -rf $my_dir/$target'_'$date'_450_shortmap_cube_cal_all.fits'
			wcsmosaic in=^$my_dir/all4.lst out=$my_dir/$target'_'$date'_450_shortmap_cube_cal_all.sdf'
			ndf2fits $my_dir/$target'_'$date'_450_shortmap_cube_cal_all.sdf' $my_dir/$target'_'$date'_450_shortmap_cube_cal_all.fits'
		else
			echo 'Only 1 450um scan present.'	
		fi
	fi
fi


#remove extra stuff that starlink creates
echo 'Cleaning up...'
echo 'Please make sure to remove temp directory created by starlink, it is usally named a string of numbers...'
#rm -rf $my_dir/adam*
#rm -rf adam*

echo '********************************************************************'
echo 'The script is finished. Please inspect the resulting data products.'
echo '********************************************************************'
