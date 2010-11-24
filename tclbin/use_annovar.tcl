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
	puts "Format is: filename, annovar directory (eg '/home/user/Downloads/annovar') "
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
		alleleSeq1 {set all1 $column}
		alleleSeq2	{set all2 $column}
	}
	incr column
}
set in [gets $fileid]
while {![eof $fileid]} {
	set line [split $in "\t"]
	#replace empty columns with '-'
	if {[string length [lindex $line $all1]] == 0} {
		set line [lreplace $line $all1 $all1 "-"]
	}
	if {[string length [lindex $line $all2]] == 0} {
		set line [lreplace $line $all2 $all2 "-"]
	}
	if {[string length [lindex $line $ref]] == 0} {
		set line [lreplace $line $ref $ref "-"]
	}
	if {[lindex $line $all1] != [lindex $line $all2]} {
		if {[lindex $line $ref] != [lindex $line $all2] && [lindex $line $ref] != [lindex $line $all1] } {
			# for easy use, we only incorperate allele1 in this case
			# for advanced use, you can incorperate both alleles, but at the end there is much more work and time putting them back together
			set new_line1 "[lindex $line $chrom] [lindex $line $begin] [lindex $line $end] [lindex $line $ref] [lindex $line $all1]"
			#set new_line2 "[lindex $line $chrom] [lindex $line $begin] [lindex $line $end] [lindex $line $ref] [lindex $line $all2]"		
			puts $fileid_out [join $new_line1 \t]
			#puts $fileid_out [join $new_line2 \t]
	
		} elseif {[lindex $line $ref] == [lindex $line $all2]} {
			set new_line1 "[lindex $line $chrom] [lindex $line $begin] [lindex $line $end] [lindex $line $ref] [lindex $line $all1]"
			puts $fileid_out [join $new_line1 \t]
		} else {
			set new_line2 "[lindex $line $chrom] [lindex $line $begin] [lindex $line $end] [lindex $line $ref] [lindex $line $all2]"		
			puts $fileid_out [join $new_line2 \t]
		}
	} else {
		set new_line "[lindex $line $chrom] [lindex $line $begin] [lindex $line $end] [lindex $line $ref] [lindex $line $all1]"		
		puts $fileid_out [join $new_line \t]

	}
	set in [gets $fileid]
}

close $fileid
close $fileid_out

#return 0

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

	# data retrieval from .human.variant_function
	puts "retrieving data from annovar files....."
	if {[catch "open $temp.variant_function r" fileid]} {
		puts "Could not open .variant_function inputfile - $fileid"
	}

	if {[catch "open $temp.exonic_variant_function r" fileid_ex]} {
		puts "Could not open .exonic_variant_function inputfile - $fileid_ex"
	}	

	if {[catch "open $temp r" fileid_old]} {
		puts "Could not open old inputfile - $fileid_old"
		exit 1
	}
	
	if {[catch "open $temp.annot w" fileid_out]} {
		puts "Could not open outputfile - $fileid_out"
		exit 1
	}
	
	set chrom 2
	set begin 3
	set end 4
	set type 0
	set name 1
	set j 1
	set header "[lindex $db $i]_type [lindex $db $i]_name [lindex $db $i]_exonic_type [lindex $db $i]_exonic_effect"
	puts $fileid_out [join $header \t]

	if {![file exists $temp.invalid_input]} {
		
		set in [gets $fileid]
		set ex [gets $fileid_ex]
		while {![eof $fileid]} {
			# Print every 1000000 lines
			if {[expr {$j%1000000}]==0} {	
				puts "reading line $j"
			}
			set line [split $in "\t"]
			set line_ex [split $ex "\t"]
			if {[string trimleft [lindex $line_ex 0] line] == $j} {
				set newName [string map {" " "_"} [lindex $line_ex 1]]
				set newline "[lindex $line $type] [lindex $line $name] $newName [lindex $line_ex 2]"
				puts $fileid_out [join $newline \t]
				set in [gets $fileid]
				set ex [gets $fileid_ex]
			} else {
				set newline "[lindex $line $type] [lindex $line $name] {} {} "
				puts $fileid_out [join $newline \t]
				set in [gets $fileid]
			}
			incr j
		}

	} else {
	
		set in1 [gets $fileid]
		set in2 [gets $fileid_old]
		set ex [gets $fileid_ex]
		while {![eof $fileid_old]} {
			if {[expr {$j%1000000}]==0} {	
				puts "reading line $j"
			}
			set line1 [split $in1 "\t"]
			set line2 [split $in2 "\t"]
			set line_ex [split $ex "\t"]
			if {[lindex $line1 $begin] == [lindex $line2 1]} {
				if {[string trimleft [lindex $line_ex 0] line] == $j} {
					set newName [string map {" " "_"} [lindex $line_ex 1]]
					set newline "[lindex $line1 $type] [lindex $line1 $name] $newName [lindex $line_ex 2]"
					puts $fileid_out [join $newline \t]
					set in1 [gets $fileid]
					set in2 [gets $fileid_old]
					set ex [gets $fileid_ex]
					
				} else {
					set newline "[lindex $line1 $type] [lindex $line1 $name] { } { }"
					puts $fileid_out [join $newline \t]
					set in1 [gets $fileid]
					set in2 [gets $fileid_old]
				}
			} else {
				set newline {" " " " " " " "}
				puts $fileid_out [join $newline \t]
				set in2 [gets $fileid_old]
			} 
			incr j	
		}
	}		
	close $fileid
	close $fileid_old
	close $fileid_ex
	close $fileid_out

	if {[catch "exec cp $temp.annot $temp.annot.[lindex $db $i] " errmsg]} {
		puts "Copy annotated info to file failed - $errmsg"
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

if {[catch "exec $path/annotate_variation.pl -filter -zerostart -dbtype avsift $temp $path/humandb/" errmsg]} {
	puts "annovar succeeded"
}
# temp file bekijken TEST
if {[catch "exec cp $temp /home/annelies/avsift_input " errmsg]} {
	puts "Copy annotated info to file failed - $errmsg"
}
if {[catch "exec cp $temp.hg18_avsift_dropped /home/annelies/avsift_dropped " errmsg]} {
	puts "Copy annotated info to file failed - $errmsg"
}

# data retrieval 
puts "retrieving data from annovar files....."
if {[catch "open $temp.hg18_avsift_dropped r" fileid]} {
	puts "Could not open .hg18_avsift_dropped inputfile - $fileid"
}
	
if {[catch "open $temp r" fileid_old]} {
	puts "Could not open old inputfile - $fileid_old"
	exit 1
}
	
if {[catch "open $temp.annot w" fileid_out]} {
	puts "Could not open outputfile - $fileid_out"
	exit 1
}

set header "avsift_score"
puts $fileid_out $header


set in1 [gets $fileid]
set in2 [gets $fileid_old]
set j 1

while {![eof $fileid_old]} {
	if {[expr {$j%1000000}]==0} {	
		puts "reading line $j"
	}
	set line1 [split $in1 "\t"]
	set line2 [split $in2 "\t"]
	
	if {[lindex $line1 3] == [lindex $line2 1]} {
		set newline "[lindex $line1 1]"
		puts $fileid_out [join $newline \t]
		set in1 [gets $fileid]
		set in2 [gets $fileid_old]
		
	} else {
		set newline {" "}
		puts $fileid_out [join $newline \t]
		set in2 [gets $fileid_old]
	} 
	incr j	

}
close $fileid
close $fileid_old
close $fileid_out

# temp file bekijken TEST

if {[catch "exec cp $temp.annot /home/annelies/avsift_annot " errmsg]} {
	puts "Copy annotated info to file failed - $errmsg"
}


######################################
#
# Pasting it all together
#
######################################

puts "pasting files to a temporary output file -> [lindex [split [lindex $argv 0] .] 0]_annovar.temp"
if {[catch "exec paste $temp.annot.knownGene $temp.annot.ensGene $temp.annot > [lindex [split [lindex $argv 0] .] 0]_annovar.temp" errmsg]} {		
	puts "Pasting annotated info to inputfile failed - $errmsg"
}

tempfile clean
return 0