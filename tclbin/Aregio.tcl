#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral

########################################################
# 
#  This is a script that calls on the program
#  that does the region based annotation
#
#
#  Written by Annelies Cassiers
#
########################################################




if {[llength $argv] < 4} {
	puts "Format is: filename, outputfilename, completegenomics directory (eg '/home/user/dev/completegenomics'), database1 (eg /home/user/Databases/hg18_clean/tRNAs.tsv) ?db2 db3.... db? "
	exit 1
}
puts "Preparing input...."
set outpath [file dirname [lindex $argv 0]]/temp_[file tail [file rootname [lindex $argv 0]]]
set path [lindex $argv 2]

set tmp [tempfile get]
if {[catch "exec ./prepareInput.tcl [lindex $argv 0] $tmp" errmsg]} {
	puts "could not execute prepareInput.tcl - $errmsg"
	exit 1
}


if {[catch {exec cg select -h $tmp} headers]} {
	puts "error getting header from $tmp - $headers"
	exit 1	
}


set column 0
foreach head $headers {
	switch -glob $head {
		*chr* {set input_chrom_column $column}
		*CHR* {set input_chrom_column $column }
		*begin* {set input_begin_column $column}
		*end* {set input_end_column $column}
	
	}
	incr column
}


set i 3
set input_c ""

while { $i < $argc} {
	puts "Examening database [lindex $argv $i]....."
	if {[catch {exec cg select -h [lindex $argv $i]} headers]} {
		puts "error getting header from [lindex $argv $i] - $headers"
		exit 1	
	}
	set column 0
	#Not all files have the following columns, 
	#this way we can leave the colums out at the end without an error.
	set name_column 100			        		
	set score_column 100
	set allel_column 100

	foreach head $headers {
		switch -glob $head {
			*chrom {set chrom_column $column}
			chr*begin {set begin_column $column}
			chr*Start {set begin_column $column}
			chr*End {set end_column $column}
			name {set name_column $column}	
			score {set score_column $column}
			allel_freq {set allel_column $column}
		}
		incr column
	}
	lappend input_c [lindex $argv $i] $chrom_column $begin_column $end_column $name_column $score_column $allel_column 
	incr i
}


# Reg_analyze lopen over alle db's


puts "Running reg_analyze.c ....."
puts "This may take a while......"
puts "The log file will be written at [file dirname [lindex $argv 0]]/log_Aregio.txt.... "


if {![file isdirectory $outpath]} {
	file mkdir "$outpath"
	set k 100
}
exec $path/bin/reg_analyze $outpath/ $tmp $input_chrom_column $input_begin_column $input_end_column  {*}$input_c >& [file dirname [lindex $argv 0]]/log_Aregio.txt

puts "Done running reg_analyze.c!"
puts "Pasting the output files together to one annotation file....."
	

# Paste all files in directory /temp/
if {[catch {exec paste [lindex $argv 0]  [glob -nocomplain -directory $outpath *.tsv]	> [lindex $argv 1] } errmsg]} {
	puts "Something went wrong while creating output file - $errmsg"
}



# delete all dbfiles in directory /temp/ to prevent another pasting
# and delete directory if he didn't exist before this script
if {[info exists k]} {
	file delete -force $outpath/
} else {
	set i 3
	while { $i < $argc} {
		file delete -force $outpath/[file tail [lindex $argv $i]] 
		incr i
	}
}

tempfile clean
exit 0