#!/bin/bash
#
# compare two pi files, ignoring newlines and whitespaces
#
if [ $# -ne 2 ]
then
	echo $0: usage: $0: PI_file1 PI_file2
	exit 1
fi
pi_file1=$1
pi_file2=$2
pi_file1_temp=/tmp/`echo $pi_file1 | sed -e s'?/?-?g'`.$$.txt
pi_file2_temp=/tmp/`echo $pi_file2 | sed -e s'?/?-?g'`.$$.txt
tr -d ' \n' < $pi_file1 > $pi_file1_temp
tr -d ' \n' < $pi_file2 > $pi_file2_temp
echo cmp $pi_file1_temp $pi_file2_temp
cmp $pi_file1_temp $pi_file2_temp
/bin/rm -f $pi_file1_temp $pi_file2_temp
