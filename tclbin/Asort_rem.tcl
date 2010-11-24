#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that sorts a given database on chromEnd,
#  replaces resp. chrom X,Y,M with 23, 24, 25.
#  and removes possible random chromosome regions
#
#  Written by Annelies Cassiers
#
##############################################################




if {[llength $argv] < 2} {
	puts "Format is: filename, output filename "
	exit 1
}


if {[catch "open [lindex $argv 0] r" fileid]} {
	puts "Could not open input file - $fileid"
	exit 1
}
set in [gets $fileid]
set headers [split $in "\t"]
set temp [tempfile get]


if {[catch "open $temp a+" fileid_temp]} {
	puts "Could nog open temp file - $fileid_temp"
	exit 1
}

set column 0
foreach head $headers {
	switch -glob -nocase $head {
		chr?? {set chrom $column}
	}
	incr column
}

puts $fileid_temp [join $headers \t]


# read real first line 
set in [gets $fileid]

while {![eof $fileid]} {
	set line [split $in "\t"]
	if {![regexp "_" [lindex $line $chrom]]} {
		set char [regexp -inline -- {[0-9]+|X|Y|M} [lindex $line $chrom]]	
		if {$char == "X"} {
			set char 23
		}
		if {$char == "Y"} {
			set char 24
		}
		if {$char == "M"} {
			set char 25
		}
		set line [lreplace $line $chrom $chrom $char]
		puts $fileid_temp [join $line \t]
		set in [gets $fileid]
	} else {
		set in [gets $fileid]
	}

}
	
close $fileid
close $fileid_temp

if {[catch {exec cg select -s {chrom chromEnd} $temp [lindex $argv 1]} errmsg]} {
	puts "Something went wrong while sorting - $errmsg"
}

tempfile clean
return 0


