#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that sorts the input file on chromEnd,
#  replaces resp. chrom X,Y,M with 24, 25, 23.
#  removes possible random chromosome regions
#  and only selects the first few necessery columns
#
#  Written by Annelies Cassiers
#
##############################################################

if {[llength $argv] < 2} {
	puts "Format is: filename, output filename, chromColumnNumber, BeginColumnNumber, EndColumnNumber "
	exit 1
}
#[file rootname [lindex $argv 0]]_prepared.tsv

set temp [tempfile get]
set temp2 [tempfile get]

if {[catch "exec cut -f[lindex $argv 2],[lindex $argv 3],[lindex $argv 4] [lindex $argv 0] > $temp" errmsg]} {
	puts "could not select first columns from inputfile - $errmsg"
	return 1
}

if {[catch "open $temp r" fileid]} {
	puts "Could not open input file - $fileid"
	exit 1
}

if {[catch "open $temp2 a+" fileid_temp]} {
	puts "Could nog open temp file - $fileid_temp"
	exit 1
}

set in [gets $fileid]
set headers {chromosome begin end}
puts $fileid_temp [join $headers \t]


# read real first line 
set in [gets $fileid]

while {![eof $fileid]} {
	set line [split $in "\t"]
	if {![regexp "_" [lindex $line 0]]} {
		set char [regexp -inline -- {[0-9]+|X|Y|M} [lindex $line 0]]	
		if {$char == "X"} {
			set char 24
		}
		if {$char == "Y"} {
			set char 25
		}
		if {$char == "M"} {
			set char 23
		}
		set line [lreplace $line 0 0 $char]
		puts $fileid_temp [join $line \t]
		set in [gets $fileid]
	} else {
		set in [gets $fileid]
	}

}
	
close $fileid
close $fileid_temp

if {[catch {exec cg select -s {chromosome end} $temp2 [lindex $argv 1]} errmsg]} {
	puts "Something went wrong while sorting - $errmsg"
}

tempfile clean
return 0
