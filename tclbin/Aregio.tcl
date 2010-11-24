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




if {[llength $argv] < 2} {
	puts "Format is: filename, outputfilename, directory (eg '/home/user'), completegenomics directory (eg '/home/user/dev/completegenomics') "
	exit 1
}
puts "Preparing input...."

set tmp [tempfile get]
if {[catch "exec ./prepareInput.tcl [lindex $argv 0] $tmp" errmsg]} {
	puts "could not execute prepareInput.tcl - $errmsg"
	exit 1
}

set dirpath [lindex $argv 2]
set path [lindex $argv 3]
if {[catch {exec cg select -h [lindex $argv 0]} headers]} {
	puts "error getting header from [lindex $argv 0] - $headers"
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


set db "${dirpath}/Databases/hg18_clean/phastConsElements44way.tsv \
	${dirpath}/Databases/hg18_clean/tfbsConsSites.tsv \
	${dirpath}/Databases/hg18_clean/omimGene.tsv \
	${dirpath}/Databases/hg18_clean/rnaGene.tsv \
	${dirpath}/Databases/hg18_clean/tRNAs.tsv \
	${dirpath}/Databases/hg18_clean/oreganno.tsv \
	${dirpath}/Databases/hg18_clean/targetScanS.tsv \
	${dirpath}/Databases/hg18_clean/evofold.tsv \
	${dirpath}/Databases/hg18_clean/CEU.tsv \
	${dirpath}/Databases/hg18_clean/CHB+JPT.tsv \
	${dirpath}/Databases/hg18_clean/YRI.tsv \
	${dirpath}/Databases/hg18_clean/hsa.tsv \
"


set input_c ""
foreach file $db {
	if {[catch {exec cg select -h $file} headers]} {
		puts "error getting header from $file - $headers"
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
lappend input_c $file $chrom_column $begin_column $end_column $name_column $score_column $allel_column 
}


# Reg_analyze lopen over alle db's


puts "Running reg_analyze.c ....."
puts "This may take a while......."
puts "Log_file will be written at $path/temp/log_Aregio.txt.... "

exec $path/bin/reg_analyze $tmp $input_chrom_column $input_begin_column $input_end_column  {*}$input_c >& $path/temp/log_Aregio.txt

puts "Done running reg_analyze.c!"
puts "Pasting the output files together to one annotation file....."


if {[catch {exec paste [lindex $argv 0] $path/temp/phastConsElements44way.tsv \
	$path/temp/tfbsConsSites.tsv \
	$path/temp/omimGene.tsv \
	$path/temp/rnaGene.tsv \
	$path/temp/tRNAs.tsv \
	$path/temp/oreganno.tsv \
	$path/temp/targetScanS.tsv \
	$path/temp/evofold.tsv \
	$path/temp/CEU.tsv \
	$path/temp/CHB+JPT.tsv \
	$path/temp/YRI.tsv \
	$path/temp/hsa.tsv \
	[lindex [split [lindex $argv 0] .] 0]_annovar.temp \
 	> [lindex $argv 1] } errmsg]} {
	puts "Something went wrong while creating output file - $errmsg"
}



file delete -force   $path/temp/phastConsElements44way.tsv \
	$path/temp/tfbsConsSites.tsv \
	$path/temp/omimGene.tsv \
	$path/temp/rnaGene.tsv \
	$path/temp/tRNAs.tsv \
	$path/temp/oreganno.tsv \
	$path/temp/targetScanS.tsv \
	$path/temp/evofold.tsv \
	$path/temp/CEU.tsv \
	$path/temp/CHB+JPT.tsv \
	$path/temp/YRI.tsv \
	$path/temp/hsa.tsv \
	[lindex [split [lindex $argv 0] .] 0]_annovar.temp

tempfile clean
exit 0