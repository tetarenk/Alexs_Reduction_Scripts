#!/bin/bash

#MIRIAD script to read in pre-processed SMA data and do antenna (optional) and Tsys corrections.
#This is step 1/2 in converting SMA data to a CASA MS through MIRIAD and its output is used 
#by the fits2casa.py (step 2/2).
#INPUT: parameters defined in section below, raw SMA data, antennas file (optional)
#OUTPUT: (1) [lsb_name].tsys.spw_*.fits, [usb_name].tsys.spw_*.fits
#		 	 in spw_fits_files directory. These data sets are corrected for antenna 
#			 positions (optional) and Tsys.
#Written by: Alex J Tetarenko
#Last Updated: Jan 9 2017

############################
#Defining variables section
############################
#sma_type='swarm'
my_dir=/mnt/bigdata/tetarenk/test_scripts
script_dir=/home/ubuntu/CASA_reduction_scripts
end_spw=50
do_ant='y'
antennas=$my_dir/antennas_150617.txt
#vars from csh scripts
#dt=150519
#rx=0
#output data sets from pre-processing scripts.
lsb_name=150617_130601_BHXRB_lsb
usb_name=150617_130601_BHXRB_usb
############################


#start miriad; must be path from bashrc!!
source /home/ubuntu/miriad_home/miriad/miriad_cvs/miriad_start.sh

#run pre-processing script-not available for users outside cfa yet!
#if [ $sma_type = 'swarm' ]
#then
	#$script_dir/swarmload.csh
#elif [ $sma_type = 'asic' ]
#then
	#$script_dir/asicload.csh
#else
	#echo 'Please enter swarm or asic'
#fi
#untar data sets output from pre-processing script
#tar xvf $my_dir/$dt'_'$rx.lsb.tar
#tar xvf $my_dir/$dt'_'$rx.usb.tar
read -r -p 'Please press enter when ready to continue >>> '

#antenna corrections if needed
if [ $do_ant = 'y' ]
then
	uvedit vis=$my_dir/$lsb_name apfile=$antennas out=$my_dir/$lsb_name.ant options=sma
	read -r -p 'Please press enter when ready to continue >>> '
	uvedit vis=$my_dir/$usb_name apfile=$antennas out=$my_dir/$usb_name.ant options=sma
	read -r -p 'Please press enter when ready to continue >>> '
#tsys correction
	smafix vis=$my_dir/$lsb_name.ant out=$my_dir/$lsb_name.ant.tsys device=/xs xaxis=time yaxis=systemp nxy=2,4 options=tsyscorr
	smafix vis=$my_dir/$usb_name.ant out=$my_dir/$usb_name.ant.tsys device=/xs xaxis=time yaxis=systemp nxy=2,4 options=tsyscorr
	read -r -p 'Please press enter when ready to continue >>> '
#write out fits files for individual spws
	for i in $(seq 0 $end_spw) 
	do 
		fits in=$my_dir/$lsb_name.ant.tsys out=$my_dir/$lsb_name.ant.tsys.spw'_'$i.fits op=uvout select=win\($i\)
		fits in=$my_dir/$usb_name.ant.tsys out=$my_dir/$usb_name.ant.tsys.tsys.spw'_'$i.fits op=uvout select=win\($i\)
	done
	read -r -p 'Please press enter when ready to continue >>> '
#put in directory
	mkdir $my_dir/spw_fits_files
	read -r -p 'Please press enter when ready to continue >>> '
	mv $my_dir/$lsb_name.ant.tsys.spw'_'*.fits $my_dir/spw'_'fits'_'files
	mv $my_dir/$usb_name.ant.tsys.spw'_'*.fits $my_dir/spw'_'fits'_'files
else
	#tsys correction
	echo 'no'
	read -r -p 'Please press enter when ready to continue >>> '
	smafix vis=$my_dir/$lsb_name out=$$my_dir/$lsb_name.tsys device=/xs xaxis=time yaxis=systemp nxy=2,4 options=tsyscorr
	smafix vis=$my_dir/$usb_name out$my_dir/$usb_name.tsys device=/xs xaxis=time yaxis=systemp nxy=2,4 options=tsyscorr
#write out fits files for individual spws
	for i in $(seq 0 $end_spw) 
	do 
		fits in=$my_dir/$lsb_name.tsys out=$my_dir/$lsb_name.tsys.spw'_'$i.fits op=uvout select=win\($i\)
		fits in=$my_dir/$usb_name.tsys out=$my_dir/$usb_name.tsys.spw'_'$i.fits op=uvout select=win\($i\)
	done
#put in directory
	mkdir $my_dir/spw_fits_files
	mv $my_dir/$lsb_name.tsys.spw'_'*.fits $my_dir/spw'_'fits'_'files
	mv $my_dir/$usb_name.tsys.spw'_'*.fits $my_dir/spw'_'fits'_'files
fi
