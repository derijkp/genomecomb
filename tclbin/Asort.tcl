#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that sorts a given database on chromEnd,
#  replaces resp. chrom X,Y,M with 24, 25, 23.
#  and removes possible random chromosome regions
#
#  Written by Annelies Cassiers
#
##############################################################




if {[llength $argv] < 2} {
	puts "Format is: directory (eg '/home/user'), build_version (eg 'hg18') "
	exit 1
}
set path [lindex $argv 0]
set build [lindex $argv 1]
set contents [glob -nocomplain -directory ${path}/Databases/s_${build} *.tsv]
if {![file isdirectory ${path}/Databases/${build}_clean]} {
	file mkdir "${path}/Databases/${build}_clean"
}
foreach file $contents {
	
	puts "Sorting $file......."
	if {[catch "open $file r" fileid]} {
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
				set char 24
			}
			if {$char == "Y"} {
				set char 25
			}
			if {$char == "M"} {
				set char 23
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

	if {[catch {exec cg select -s {chrom chromEnd} $temp ${path}/Databases/${build}_clean/[file tail $file]} errmsg]} {
		puts "Something went wrong while sorting - $errmsg"
	}

	tempfile clean

}

return 0


