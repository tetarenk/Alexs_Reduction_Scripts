#! /bin/csh -f
#
#jhz 2015-07-20
#jhz 2015-07-22
#jhz 2015-07-23
#jhz 2015-07-24
#jhz 2015-07-24 for swarm+asic data
echo "                                                                      "
echo "######################################################################"
echo "# swarmload.csh: script to convert SMA data produced from hybrid      "
echo "# correlator in new format into Miriad format, averaging the SWARM    "
echo "# data every two channels. The output data sets of both ASIC and      "
echo "# SWARM data can be further reduced with the SMA Miriad released      "
echo "# prior to 2014-11-15 developed for ASIC correlator with 4 GHz BW     "
echo "# check web site below for a stable version of SMA Miriad 1.5.1:      "
echo "# click here for download (RH 6/CentOS 6)                             " 
echo "######################################################################"
echo "                                                                      "
############################################################################
#
#begin of users' input
#
#set the name of the SMA raw data
#
set file = 150519_09:50:17
echo $file
#
#the observing date for the prefix of the output file names 
#
set dt    = 150519 
#
#set the data path, e.g. for RTDC in Cambridge, MA
#
set DATPATH = /2014/science/mir_data/ 
#
#end of users' input
#
#####################################
#                                   #
# The pre-processing routine starts #
#                                   #
#####################################
#
#set the receiver
#
set rxid  = 0       #0 for 230 GHz
                    #1 for 345 GHz
                    #2 for 400 GHz
                    #3 for 690 GHz
set rx    = rx${rxid}
#
#set the sideband
#
set sb    = lsb
set nn    = 0       #0 for lsb 
                    #1 for usb
#
#the handle of switch sideband
#
OTHERSIDEBAND:
set fname = "$dt""_""$rx".$sb
#
#start the pipeline to convert new SMA format to Miriad 
#
SMALOD:
echo "Loading "$sb "data: "$fname
\rm -fr $fname
smalod in=${DATPATH}/$file out=$dt rxif=$rxid \
   options=circular sideband=$nn
#
#separating ASIC data
#
SPLIT-ASIC:
echo "splitting ASIC data to: "$fname.asic
\rm -rf $fname.asic
swarmsplt vis=$fname  out=$fname.asic options=asic 
#
#separating the 1st SWARM chunk & reducing spectral resolution by a factor of 2
#
SPLIT-SWARM1:
echo "splitting data from SWARM chunk1 (s49) to: "$fname.sw1
\rm -rf $fname.sw1
swarmsplt vis=$fname  out=$fname.sw1 options=sw1 line=channel,8192,1,2
#
#separating the 2nd SWARM chunk & reducing spectral resolution by a factor of 2
#
SPLIT-SWARM2:
echo "splitting data from SWARM chunk2 (s50) to: "$fname.sw2
\rm -rf $fname.sw2
swarmsplt vis=$fname  out=$fname.sw2 options=sw2 line=channel,8192,1,2
echo " "
#
#switching sideband from lsb to usb
#
if ($sb == lsb) then
set sb = usb
set nn = 1
goto OTHERSIDEBAND
endif
#
#reporting the output files from this pre-process
#
echo " "
echo "Successfully converted SMA ASIC and SWARM into Miriad:"
echo "                         "
echo "1) $dt""_""$rx"".lsb.asic"
echo "2) $dt""_""$rx"".usb.asic"
echo "3) $dt""_""$rx"".lsb.sw1 "
echo "4) $dt""_""$rx"".usb.sw1 "
echo "3) $dt""_""$rx"".lsb.sw2 "
echo "4) $dt""_""$rx"".usb.sw2 "
echo "                         "
#
#tar the subsets of the output files directories
#
TAR:
echo "Tar the data files       "
\rm -rf $dt"_"$rx.lsb.tar $dt"_"$rx.usb.tar 
echo "                         "
tar cvf $dt"_"$rx.lsb.tar $dt"_"$rx.lsb.asic $dt"_"$rx.lsb.sw1 $dt"_"$rx.lsb.sw2
tar cvf $dt"_"$rx.usb.tar $dt"_"$rx.usb.asic $dt"_"$rx.usb.sw1 $dt"_"$rx.usb.sw2
TARREPORT:
echo " "
echo "Successfully tarred the preprocessed the subsets of data into:"
echo " "
echo "LSB data: $dt""_""$rx"".lsb.tar" 
echo "USB data: $dt""_""$rx"".usb.tar"
echo " "                      
#
#clean up the intermediate working files
#
echo "Cleaning the working files":
echo ""
echo "$dt""_""$rx.lsb.asic"
echo "$dt""_""$rx.lsb.sw1"
echo "$dt""_""$rx.lsb.sw2"
echo "$dt""_""$rx.usb.asic"
echo "$dt""_""$rx.usb.sw1"
echo "$dt""_""$rx.usb.sw2"
#
\rm -rf $dt"_"$rx.lsb.asic $dt"_"$rx.usb.asic $dt"_"$rx.lsb.sw1 $dt"_"$rx.usb.sw1 $dt"_"$rx.lsb.sw2 $dt"_"$rx.usb.sw2
echo " "
echo "End the pre-processing successfully!"
echo " "
##########################
#                        #
# End of the pre-process #
#                        #
##########################
exit