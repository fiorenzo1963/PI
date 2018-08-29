#!/bin/bash
#
# compare two pi files, ignoring newlines and whitespaces
#
if [ $# -ne 2 ]
then
	echo $0: usage: $0: PI_file1 PI_file2
fi
pi_file1=$1
pi_file2=$2
tr -d ' \n' < $pi_file1 > /tmp/pi_file1.$$.str
tr -d ' \n' < $pi_file2 > /tmp/pi_file2.$$.str
cmp /tmp/pi_file1.$$.str /tmp/pi_file2.$$.str
/bin/rm -f /tmp/pi_file1.$$.str /tmp/pi_file2.$$.str
