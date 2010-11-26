#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


#########################################################################
# 
#  This is a script that sorts a given database on chromEnd
#  and removes possible overlap of regions.
#  Removal of overlap can be done by taking only the highest
#  scoring region (this is always done when score is available)
#  or taking all regions in 1 line (if score is not available) 
#
#  Written by Annelies Cassiers
#
########################################################################


if {[llength $argv] < 2} {
	puts "Format is: directory (eg '/home/user'), build_version (eg 'hg18') "
	exit 1
}
set path [lindex $argv 0]
set build [lindex $argv 1]
set contents [glob -nocomplain -directory ${path}/Databases/${build} *.tsv]
foreach file $contents {

	puts "Sorting $file ....."
	puts "This may take a few minutes....."

	set temp [tempfile get]

	if {[catch {exec ./Asort_rem.tcl $file $temp} errmsg]} {
		puts "Sorting failed - $errmsg"
		exit 1
	}

	if {[catch "open $temp r" fileid]} {
		puts "Could not open sorted file - $fileid"
		exit 1
	}

	if {[catch "open ${path}/Databases/${build}_clean/[file tail $file] w" fileid_out]} {
		puts "Could not open outputfile - $fileid_out"
		exit 1
	}

	seek $fileid 0 

	set in [gets $fileid]
	set headers [split $in "\t"]
	puts $fileid_out [join $headers \t]

	set column 0
	foreach head $headers {
		switch -glob $head {
			chr*Start {set begin $column}
			chr*End {set end $column}
			score {set score $column}
			name {set name $column}
			allel_freq	{set name $column}
		}
		incr column
	}

	set i 1

	if {[info exists score]} {
		
		set in1 [gets $fileid]
		set line1 [split $in1 "\t"]
		
		puts "Removing overlap with score......"
		puts "reading line $i "
		 
		while {![eof $fileid]} {
			incr i
			set in2 [gets $fileid]
			set line2 [split $in2 "\t"]
			set line1 [split $in1 "\t"]
		
			
			# Print every 100000 lines
			if {[expr {$i%100000}]==0} {	
				puts "reading line $i"
			}	

			if {[lindex $line1 $begin] < [lindex $line2 $begin]} {
				if {[lindex $line1 $end] == [lindex $line2 $end]} {
					set line1_output [lreplace $line1 $end $end [lindex $line2 $begin]]
					if {[lindex $line1 $score]>[lindex $line2 $score]} {
						set line2_output [join $line1 \t] 
					} else {
						set line2_output [join $line2 \t]
					} 
					set line2_output [lreplace $line2_output $begin $begin [lindex $line2 $begin]]
					set line2_output [lreplace $line2_output $end $end [lindex $line2 $end]]
					puts $fileid_out [join $line1_output \t]
					set in1 [join $line2_output \t]
				} elseif {[lindex $line1 $end] < [lindex $line2 $end] && [lindex $line2 $begin] > [lindex $line1 $end] } {
					puts $fileid_out [join $line1 \t]
					set in1 [join $line2 \t]
				} elseif {[lindex $line1 $end] < [lindex $line2 $end] && [lindex $line2 $begin] == [lindex $line1 $end]} {
					puts $fileid_out [join $line1 \t]
					set in1 [join $line2 \t]				
				} else {
					if {[lindex $line1 $score]>[lindex $line2 $score]} {
						set line2_output [lreplace $line2 $begin $begin [lindex $line1 $end]]
						puts $fileid_out [join $line1 \t]
						set in1 [join $line2_output \t]
					} else {
						set line1_output [lreplace $line1 $end $end [lindex $line2 $begin]]
						puts $fileid_out [join $line1_output \t]
						set in1 [join $line2 \t]
					} 		
				}
			} elseif {[lindex $line1 $begin] > [lindex $line2 $begin]} {	
				if {[lindex $line1 $end] == [lindex $line2 $end]} {
					set line1_output [lreplace $line2 $end $end [lindex $line1 $begin]]
					if {[lindex $line1 $score]>[lindex $line2 $score]} {
						set line2_output [join $line1 \t] 
					} else {
						set line2_output [join $line2 \t]
					} 			
					set line2_output [lreplace $line2_output $begin $begin [lindex $line1 $begin]]
					set line2_output [lreplace $line2_output $end $end [lindex $line2 $end]]
					puts $fileid_out [join $line1_output \t]
					set in1 [join $line2_output \t]
				} elseif {[lindex $line1 $begin] > [lindex $line2 $end]} {
					puts $fileid_out [join $line1 \t]
					set in1 [join $line2 \t]
				} else {
					if {[lindex $line1 $score]>[lindex $line2 $score]} {
						set line1_output [lreplace $line2 $end $end [lindex $line1 $begin]]
						set line2_output [join $line1 \t]
						set line3_output [lreplace $line2 $begin $begin [lindex $line1 $end]]
						puts $fileid_out [join $line1_output \t]
						puts $fileid_out [join $line2_output \t]
						set in1 [join $line3_output \t]
					} else {
						set in1 [join $line2 \t]
					} 
				}
			} elseif {[lindex $line1 $begin] == [lindex $line2 $begin]} {
				if {[lindex $line1 $end] == [lindex $line2 $end]} {
					if {[lindex $line1 $score]>[lindex $line2 $score]} {
						set line1_output [join $line1 \t] 
					} else {
						set line1_output [join $line2 \t]
					} 			
					set in1 [join $line1_output \t]
				} else {
					if {[lindex $line1 $score]>[lindex $line2 $score]} {
						set line1_output [join $line1 \t] 
					} else {
						set line1_output [join $line2 \t]
					} 				
					set line1_output [lreplace $line1_output $begin $begin [lindex $line1 $begin]]
					set line1_output [lreplace $line1_output $end $end [lindex $line1 $end]]
					set line2_output [lreplace $line2 $begin $begin [lindex $line1 $end]]
					puts $fileid_out [join $line1_output \t]
					set in1 [join $line2_output \t]
				}
			}
		}
	 
	} else {
		set in1 [gets $fileid]
		set line1 [split $in1 "\t"]
		puts "Removing overlap without score......"
		puts "reading line $i "

		while {![eof $fileid]} {
			incr i
			set in2 [gets $fileid]
			set line2 [split $in2 "\t"]
			set line1 [split $in1 "\t"]
			
			# Print every 100000 lines
			if {[expr {$i%100000}]==0} {	
				puts "reading line $i"
			}	
		
			if {[lindex $line1 $begin] < [lindex $line2 $begin]} {
				if {[lindex $line1 $end] == [lindex $line2 $end]} {
					set line1_output [lreplace $line1 $end $end [lindex $line2 $begin]]
					set line2_output [lreplace $line2 $name $name "[lindex $line1 $name],[lindex $line2 $name]" ]
					puts $fileid_out [join $line1_output \t]
					set in1 [join $line2_output \t]
				} elseif {[lindex $line1 $end] < [lindex $line2 $end] && [lindex $line2 $begin] > [lindex $line1 $end] } {
					puts $fileid_out [join $line1 \t]
					set in1 [join $line2 \t]
				} elseif {[lindex $line1 $end] < [lindex $line2 $end] && [lindex $line2 $begin] == [lindex $line1 $end]} {
					puts $fileid_out [join $line1 \t]
					set in1 [join $line2 \t]				
				} else {
					set line1_output [lreplace $line1 $end $end [lindex $line2 $begin]]
					set line2_output [lreplace $line2 $end $end [lindex $line1 $end]]
					set line2_output [lreplace $line2_output $name $name "[lindex $line1 $name],[lindex $line2 $name]" ]
					set line3_output [lreplace $line2 $begin $begin [lindex $line1 $end]]
					puts $fileid_out [join $line1_output \t]
					puts $fileid_out [join $line2_output \t]
					set in1 [join $line3_output \t]
				}
			} elseif {[lindex $line1 $begin] > [lindex $line2 $begin]} {
				if {[lindex $line1 $end] == [lindex $line2 $end]} {
					set line1_output [lreplace $line2 $end $end [lindex $line1 $begin]]
					set line2_output [lreplace $line2 $begin $begin [lindex $line1 $begin]]
					set line2_output [lreplace $line2_output $name $name "[lindex $line1 $name],[lindex $line2 $name]" ]
					puts $fileid_out [join $line1_output \t]
					set in1 [join $line2_output \t]
				} elseif {[lindex $line1 $begin] > [lindex $line2 $end]} {
					puts $fileid_out [join $line1 \t]
					set in1 [join $line2 \t]
				} else {
					set line1_output [lreplace $line2 $end $end [lindex $line1 $begin]]
					set line2_output [lreplace $line1 $name $name "[lindex $line1 $name],[lindex $line2 $name]" ]
					set line3_output [lreplace $line2 $begin $begin [lindex $line1 $end]]
					puts $fileid_out [join $line1_output \t]
					puts $fileid_out [join $line2_output \t]
					set in1 [join $line3_output \t]
				}
			} elseif {[lindex $line1 $begin] == [lindex $line2 $begin]} {
				if {[lindex $line1 $end] == [lindex $line2 $end]} {
					set line1_output [lreplace $line1 $name $name "[lindex $line1 $name],[lindex $line2 $name]" ]
					set in1 [join $line1_output \t]
				} else {
					set line1_output [lreplace $line1 $name $name "[lindex $line1 $name],[lindex $line2 $name]" ]
					set line2_output [lreplace $line2 $begin $begin [lindex $line1 $end]]
					puts $fileid_out [join $line1_output \t]
					set in1 [join $line2_output \t]

				}
			}
		}
	}



	close $fileid
	close $fileid_out
	puts "Done reading line $i "
	tempfile clean

}

exit 0
