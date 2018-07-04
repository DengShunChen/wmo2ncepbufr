#!/bin/ksh
set -x
Ifile=/nwpr/gfs/a398/tmpgps/20160608/bfrPrf_C001.2016.160.00.04.G10_0001.0007_bufr

#-------------------------------------------------------------------------
export PATH=/users/xa09/x86_64/perl-5.24.0/bin:${PATH}
echo $PATH
dump_bufr.pl $Ifile > bfrPrf_C001.2016.170.00.16.G12_0001.0007_bufr.txt

exit
