#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that downloads a number of databases
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

if {![file isdirectory $path/Databases/$build]} {
	file mkdir $path/Databases/$build
}


#downloading with wget

set UCSCdb { 		{omimGene} \
		{phastConsElements44way} \
		{oreganno} \
		{rnaGene} \
		{tRNAs} \
		{tfbsConsSites} \
		{targetScanS} \
		{evofold} \
}

set i 0

while {$i < [llength $UCSCdb]} {

	if {![file isfile $path/Databases/$build/[lindex $UCSCdb $i].tsv]} {
		puts "----------------------------------------------------"
		puts "Downloading [lindex $UCSCdb $i].tsv....."
		if {[catch "exec wget --tries=45 --directory-prefix=$path/Databases/$build/ http://hgdownload.cse.ucsc.edu/goldenPath/$build/database/[lindex $UCSCdb $i].txt.gz " errmsg]} {
			puts "Downloading database succeeded"
		}
		if {[catch "exec wget --directory-prefix=$path/Databases/$build/ http://hgdownload.cse.ucsc.edu/goldenPath/$build/database/[lindex $UCSCdb $i].sql " errmsg]} {
			puts "Downloading header succeeded"
		}

		puts "Gunzipping [lindex $UCSCdb $i].tsv...."
		if {[catch "exec gunzip -f $path/Databases/$build/[lindex $UCSCdb $i].txt.gz " errmsg]} {
			puts "Gunzip failed - $errmsg"
		}

		puts "Giving [lindex $UCSCdb $i].tsv the right header...."
		if {[catch "open $path/Databases/$build/[lindex $UCSCdb $i].sql r" fileid]} {
			puts "Could not open sql file - $fileid"
			exit 1
		}
		set temp [tempfile get]
		if {[catch "open $temp w" file_outid]} {
			puts "Could not open sql file - $file_outid"
			exit 1
		}
	
		set line [gets $fileid]
		set header ""

		while {![eof $fileid]} {
			if {[lsearch $line "CREATE"] != -1} {
				set line [gets $fileid]
				set head [string trim [lindex $line 0] (`) ]
				while {$head != "KEY"} {	
					lappend header $head 
					set line [gets $fileid]
					set head [string trim [lindex $line 0] (`) ]
				}
			} 
			set line [gets $fileid]

		}	

		puts $file_outid [join $header \t]
		close $fileid
		close $file_outid

		if {[ catch "exec cat $temp $path/Databases/$build/[lindex $UCSCdb $i].txt > $path/Databases/$build/[lindex $UCSCdb $i].tsv" errmsg]} {
			puts "Something went wrong while adding header - $errmsg"
		}

		file delete -force $path/Databases/$build/[lindex $UCSCdb $i].txt $path/Databases/$build/[lindex $UCSCdb $i].sql

		puts "----------------------------------------------------"

	} else {
		puts "----------------------------------------------------"
		puts "The file '[lindex $UCSCdb $i].tsv' already exists..."
		puts "Skipping the download..."
		puts "----------------------------------------------------"
	}

	incr i
}



tempfile clean
exit 0