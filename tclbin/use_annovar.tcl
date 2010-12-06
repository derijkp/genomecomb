#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


#########################################################################
# 
#  This is a script that uses Annovar to do gene-based annotation
#  
#
#  Written by Annelies Cassiers
#
########################################################################


if {[llength $argv] < 2} {
	puts "Format is: filename, annovar directory (eg '/home/user/Downloads/annovar') ?-makeFile?"
	exit 1
}
# default setting is not making an output annotationfile
# -makeFile setting is making an output annotationfile

if {[lindex $argv 2] == "-makeFile"} {
	puts "Annotation file will be created."
} elseif {[lindex $argv 2] == "" } {
	puts "No annotation file will be created."
	puts "You can paste the output file in your command line to the input annotationfile."
} else {
	puts "Please use -makeLine if you wish to create an annotation file"
	puts "Otherwise use only 2 arguments"
	exit 1
}

set path [lindex $argv 1]

########################################
#
# Making a new input file from the old
#
########################################

puts "Making the right input file...."

if {[catch "open [lindex $argv 0] r" fileid]} {
	puts "Could not open input file - $fileid"
	exit 1
}

set temp [tempfile get]
if {[catch "open $temp w" fileid_out]} {
	puts "Could not open new inputfile - $fileid_out"
	exit 1
}

set in [gets $fileid]
set headers [split $in "\t"]
set column 0
foreach head $headers {
	switch -glob $head {
		chrom* {set chrom $column} 
		begin {set begin $column}
		end {set end $column}
		reference {set ref $column}
		type {set type $column}
	}
	incr column
}

set allel_place [lsearch -all $headers "alleleSeq*"]
set allel_count [llength $allel_place]

set in [gets $fileid]
while {![eof $fileid]} {
	set line [split $in "\t"]
	#replace empty columns with '-'
	foreach all $allel_place {
		if {[string length [lindex $line $all]] == 0} {
			set line [lreplace $line $all $all "-"]
		} elseif {[lindex $line $all] == "?"} {
			set line [lreplace $line $all $all "0"]
		}
	}

	if {[string length [lindex $line $ref]] == 0} {
		set line [lreplace $line $ref $ref "-"]
	} elseif {[lindex $line $ref] == "?"} {
		set line [lreplace $line $ref $ref "0"]
	}
	
	
	set i 0
	set allel {}
	while {$i < $allel_count} {
		lappend allel [lindex $line [lindex $allel_place $i]]
		incr i
	}

	set use_allels [lsort -unique $allel]
	set i 0
	while { $i < [llength $use_allels]}	{
		if {[lindex $line $ref] != [lindex $use_allels $i]} {
			#make new line
			set new_line "[lindex $line $chrom] [lindex $line $begin] [lindex $line $end] [lindex $line $ref] [lindex $use_allels $i] [lindex $line $type] " 
			puts $fileid_out [join $new_line \t]
		} else {
			#discard allele
		}
		incr i
	}	

	set in [gets $fileid]

}

close $fileid
close $fileid_out


###########################################
#
# Annovar laten lopen
#
###########################################

set db { {knownGene} {ensGene} }
set i 0

while {$i < [llength $db]} {

	# annovar lopen
	puts "Running annovar for [lindex $db $i] database ....."
	puts "this may take a while......"
	if {[catch "exec $path/annotate_variation.pl -geneanno -zerostart -dbtype [lindex $db $i] $temp $path/humandb/" errmsg]} {
		puts "annovar succeeded"
	}

	set typeA 0
	set nameA 1
	set chromA 2
	set beginA 3
	set endA 4
	set nucA 6
	set commentA 7

	
	set which_file { {variant_function} {exonic_variant_function} }
	set j 0
	while {$j < [llength $which_file]} {
		# data retrieval from annovar file
		puts "retrieving data from annovar file [lindex $db $i].[lindex $which_file $j]....."

		if {[catch "open $temp.[lindex $which_file $j] r" fileid]} {
			puts "Could not open .[lindex $which_file $i] inputfile - $fileid"
			exit 1
		}
		
		if {[catch "open [lindex $argv 0] r" fileid_old]} {
			puts "Could not open inputfile - $fileid_old"
			exit 1
		}

		if {[catch "open $temp.annot w" fileid_out]} {
			puts "Could not open outputfile - $fileid_out"
			exit 1
		}
	

		#skipping headers
		set in_old [gets $fileid_old]
	
		#reading real lines
		set in1 [gets $fileid]
		set line1 [split $in1 "\t"]
		set in_old [gets $fileid_old]
		set line_old [split $in_old "\t"]

		if {[lindex $which_file $j] == "variant_function" } {
			set header "[lindex $db $i]_type [lindex $db $i]_name"
			puts $fileid_out [join $header \t]
			set l 0
			set k 0
			while {![eof $fileid_old]} {
				if {[expr {$l%1000000}]==0} {	
					puts "reading line $l"
				}
								
				if {[lindex $line1 $chromA] == [lindex $line_old $chrom] && [lindex $line1 $beginA] == [lindex $line_old $begin] && [lindex $line1 $endA] == [lindex $line_old $end] && [lindex $line1 $commentA] == [lindex $line_old $type]} {
					set ntype [string map {" " "_"} [lindex $line1 $typeA]]
					set nname [string map {" " "_"} [lindex $line1 $nameA]]
					set newtype "$ntype"
					set newname "$nname"
					set in [gets $fileid]
					set line [split $in "\t"]
					incr k
					while {[lindex $line1 $beginA] == [lindex $line $beginA] && [lindex $line1 $endA] == [lindex $line $endA] && [lindex $line1 $commentA] == [lindex $line $commentA] }	{
						set in [gets $fileid]
						set line [split $in "\t"]
						incr k
					}
					set newline "$newtype $newname"
					puts $fileid_out [join $newline \t]
					set line1 $line
					set in_old [gets $fileid_old]
					set line_old [split $in_old "\t"]
					incr l
				} else {	
					#lijn  moet ergens van komen, dus deze houden we int geheugen, en lezen nieuwe lijn in van oude file
					# in output lege lijn plaatsen omdat er geen overeenkomst was met oude file
					set newline {" " " "}
					puts $fileid_out [join $newline \t]
					set in_old [gets $fileid_old]
					set line_old [split $in_old "\t"]
					incr l
				}
			}

		} elseif {[lindex $which_file $j] == "exonic_variant_function" } {
			set header "[lindex $db $i]_exonic_type [lindex $db $i]_exonic_effect"
			puts $fileid_out [join $header \t]

			set l 0
			set k 0
			while {![eof $fileid_old]} {
				if {[expr {$l%1000000}]==0} {	
					puts "reading line $l"
				}
								
				if {[lindex $line1 $chromA] == [lindex $line_old $chrom] && [lindex $line1 $beginA] == [lindex $line_old $begin] && [lindex $line1 $endA] == [lindex $line_old $end] && [lindex $line1 $commentA] == [lindex $line_old $type]} {
					set ntype [string map {" " "_"} [lindex $line1 $typeA]]
					set nname [string map {" " "_"} [lindex $line1 $nameA]]
					set newtype "[lindex $line1 $nucA]:$ntype"
					set newname "[lindex $line1 $nucA]:$nname"
					set in [gets $fileid]
					set line [split $in "\t"]
					incr k
					while {[lindex $line1 $beginA] == [lindex $line $beginA] && [lindex $line1 $endA] == [lindex $line $endA] && [lindex $line1 $commentA] == [lindex $line $commentA] }	{
						set ntype [string map {" " "_"} [lindex $line $typeA]]
						set nname [string map {" " "_"} [lindex $line1 $nameA]]
						lappend newtype "[lindex $line $nucA]:$ntype"
						lappend newname "[lindex $line $nucA]:$nname"
						set in [gets $fileid]
						set line [split $in "\t"]
						incr k
					}
					set n_type [join $newtype ";"]
					set n_name [join $newname ";"]
					set newline "$n_type $n_name"
					puts $fileid_out [join $newline \t]
					set line1 $line
					set in_old [gets $fileid_old]
					set line_old [split $in_old "\t"]
					incr l
				} else {	
					#lijn  moet ergens van komen, dus deze houden we int geheugen, en lezen nieuwe lijn in van oude file
					# in output lege lijn plaatsen omdat er geen overeenkomst was met oude file
					set newline {" " " "}
					puts $fileid_out [join $newline \t]
					set in_old [gets $fileid_old]
					set line_old [split $in_old "\t"]
					incr l
				}
			}

		}
		

		puts "Done reading line: $l"
		incr typeA
		incr nameA
		incr chromA 
		incr beginA
		incr endA 
		incr nucA 
		incr commentA
		
		close $fileid
		close $fileid_old
		close $fileid_out


		if {[catch "exec cp $temp.annot $temp.annot.[lindex $db $i].[lindex $which_file $j] " errmsg]} {
			puts "Copy annotated info to file failed - $errmsg"
		}
		incr j
	}
	

	incr i
}


######################################
#
# Running sift in annovar
#
######################################

# annovar lopen
puts "Running annovar for avsift database ....."
puts "this may take a while......"

if {[catch "exec $path/annotate_variation.pl -filter -sift_threshold 0 -zerostart -dbtype avsift $temp $path/humandb/" errmsg]} {
	puts "annovar succeeded"
}

######
#
# Creating header for sorting
#
######

if {[catch "open $temp.header w" fileout_header]} {
	puts "Could not open header file - $fileout_header"
	exit 1
}	
	

set headers "avsift score chrom begin end ref allel comment"
puts $fileout_header [join $headers \t]

close $fileout_header

if {[catch "exec cat $temp.header $temp.hg18_avsift_dropped > $temp.hg18_avsift_dropped_header" errmsg]} {
	puts "something went wrong while adding a header to input file - $errmsg"
	exit 1
}
if {[catch "exec cg select -s {chrom end} $temp.hg18_avsift_dropped_header $temp.hg18_avsift_dropped_sort" errmsg]} {
	puts "Something went wrong while sorting - $errmsg"
	exit 1
}


######
#
# extract info from annovar files
#
###### 

# data retrieval 
puts "retrieving data from annovar files....."
if {[catch "open $temp.hg18_avsift_dropped_sort r" fileid]} {
	puts "Could not open .hg18_avsift_dropped inputfile - $fileid"
	exit 1
}
	
if {[catch "open [lindex $argv 0] r" fileid_old]} {
	puts "Could not open old inputfile - $fileid_old"
	exit 1
}
	
if {[catch "open $temp.annot w" fileid_out]} {
	puts "Could not open outputfile - $fileid_out"
	exit 1
}

set header "avsift_score avsift_minimalScore"
puts $fileid_out [join $header \t]

set scoreA 1
set chromA 2
set beginA 3
set endA 4
set nucA 6
set commentA 7

#skipping headers
set in_old [gets $fileid_old]
set in1 [gets $fileid]

#reading real lines
set in1 [gets $fileid]
set line1 [split $in1 "\t"]
set in_old [gets $fileid_old]
set line_old [split $in_old "\t"]
				
set l 0
set k 0
while {![eof $fileid_old]} {
	if {[expr {$l%1000000}]==0} {	
		puts "reading line $l"
	}
							
	if {[lindex $line1 $chromA] == [lindex $line_old $chrom] && [lindex $line1 $beginA] == [lindex $line_old $begin] && [lindex $line1 $endA] == [lindex $line_old $end] && [lindex $line1 $commentA] == [lindex $line_old $type]} {
		set newscore "[lindex $line1 $nucA]:[lindex $line1 $scoreA]"
		set in [gets $fileid]
		set line [split $in "\t"]
		set min_score [lindex $line1 $scoreA]
		incr k
		while {[lindex $line1 $beginA] == [lindex $line $beginA] && [lindex $line1 $endA] == [lindex $line $endA] && [lindex $line1 $commentA] == [lindex $line $commentA] }	{
			lappend newscore "[lindex $line $nucA]:[lindex $line $scoreA]"
			if {$min_score > [lindex $line $scoreA]} {
				set min_score [lindex $line $scoreA]
			}
			set in [gets $fileid]
			set line [split $in "\t"]
			incr k
		}
		set n_score [join $newscore ";"]
		set newline "$n_score $min_score"
		puts $fileid_out [join $newline \t]
		set line1 $line
		set in_old [gets $fileid_old]
		set line_old [split $in_old "\t"]
		incr l
	} else {	
		set newline {" " " "}
		puts $fileid_out [join $newline \t]
		set in_old [gets $fileid_old]
		set line_old [split $in_old "\t"]
		incr l
	}
}



puts "Done reading line: $l"
		
close $fileid
close $fileid_old
close $fileid_out




######################################
#
# Pasting it all together
#
######################################

#depending on whether -makeFile is the 3th argument

if {[lindex $argv 2] == "-makeFile"} {
	puts "pasting files to the input annotation file -> [lindex [split [lindex $argv 0] .] 0]_annot_annovar.tsv"
	if {[catch "exec paste [lindex $argv 0] $temp.annot.knownGene.variant_function $temp.annot.knownGene.exonic_variant_function $temp.annot.ensGene.variant_function $temp.annot.ensGene.exonic_variant_function $temp.annot > [lindex [split [lindex $argv 0] .] 0]_annot_annovar.tsv" errmsg]} {		
		puts "Pasting annotated info to inputfile failed - $errmsg"
	}
} elseif {[lindex $argv 2] == "" } {
	puts "pasting files to a temporary output file -> [lindex [split [lindex $argv 0] .] 0]_annovar.temp"
	if {[catch "exec paste $temp.annot.knownGene.variant_function $temp.annot.knownGene.exonic_variant_function $temp.annot.ensGene.variant_function $temp.annot.ensGene.exonic_variant_function $temp.annot > [lindex [split [lindex $argv 0] .] 0]_annovar.temp" errmsg]} {		
		puts "Pasting annotated info to temporary outputfile failed - $errmsg"
	}
	puts "You can later paste the file to the annotation file. "	
}




tempfile clean
return 0