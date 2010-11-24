#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

##############################################################
# 
#  This is a script that replaces resp.
#  chrom 23, 24, 25 with X,Y,M. 
#  
#
#  Written by Annelies Cassiers
#
##############################################################

if {[llength $argv] < 2} {
	puts "Format is: filename, output filname "
	exit 1
}

if {[catch "open [lindex $argv 0] r" fileid]} {
	puts "Could not open input file - $fileid"
	exit 1
}

if {[catch "open [lindex $argv 1] w" fileid_out]} {
	puts "Could not open outputfile - $fileid_out"
	exit 1
}


set in [gets $fileid]
set headers [split $in "\t"]

set column 0
foreach head $headers {
	switch -glob -nocase $head {
		chr?? {set chrom $column}
	}
	incr column
}
puts $fileid_out [join $headers \t]


# read real first line 
set in [gets $fileid]
while {![eof $fileid]} {
	set line [split $in "\t"]
	set char [lindex $line $chrom]
	if {$char == 23} {
		set char "X"
	}
	if {$char == 24} {
		set char "Y"
	}
	if {$char == 25} {
		set char "M"
	}
	set line [lreplace $line $chrom $chrom $char]
	puts $fileid_out [join $line \t]
	set in [gets $fileid]
}
	

close $fileid
close $fileid_out

return 0