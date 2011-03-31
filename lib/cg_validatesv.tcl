#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral
package require dict
package require cindex

#######################################################
# 
#  This is a script that designs the primer pairs 
#  for sequencing breakpoints of inversions.
#
#  Written by Annelies Cassiers
#
#######################################################

set maxnum 1000
set searchGenomeDB "/home/annelies/BigDisk/complgen/refseq/hg18/build36-ssa/"
set MAX_SIZE 4000

proc cg_validatesv_getSeqRef {chr patch breakpoint read} {
	set beginTarget [expr $breakpoint - 200]
	set endTarget [expr $breakpoint + $read - $patch + 250]
	if {[catch "exec FetchSequence.tcl $chr $beginTarget $endTarget" seqTarget]} {
		puts "fetching sequences failed - $seqTarget"
		exit 1
	}
	
	set beginL [expr $beginTarget - 400]
	set endL [expr $beginTarget - 1] ;#otherwise there will be overlap
	if {[catch "exec FetchSequence.tcl $chr $beginL $endL" seqL]} {
		puts "fetching sequences failed - $seqL"
		exit 1
	}

	set beginR [expr $endTarget + 1] ;#otherwise there will be overlap
	set endR [expr $endTarget + 400]
	if {[catch "exec FetchSequence.tcl $chr $beginR $endR" seqR]} {
		puts "fetching sequences failed - $seqR"
		exit 1
	}

	return [list $seqTarget $seqL $seqR]
}

proc cg_validatesv_getSeqInv {seq1 seq2} {
	set seqTarget1 [lindex $seq1 0]
	set seqL [lindex $seq1 1]
	set seqTarget2 [lindex $seq2 0]
	set seqR [string_reverse [string_change [lindex $seq2 1] [list A T G C C G T A]]]
	set seqTargetA [string range $seqTarget1 0 [expr [string length $seqTarget1]/2]]
	set subseq1 [string range $seqTarget2 0 [expr [string length $seqTarget2]/2]]
	set seqTargetB [string_reverse [string_change $subseq1 [list A T G C C G T A]]]
	set seqTarget $seqTargetA$seqTargetB

	return [list $seqTarget $seqL $seqR]
}


proc cg_validatesv_mask_seq {list_seq EVAL MIN} {
	set seqTarget [lindex $list_seq 0]
	set seqL [lindex $list_seq 1]
	set seqR [lindex $list_seq 2]
	
	# Search for homopolymer stretches
	# --------------------------------
	
	set len_stretch {} ; #so user knows how much nucleotides are located in stretches
	set hits [regexp -inline -indices -all {A{10,}|G{10,}|T{10,}|C{10,}} $seqL]
	foreach hit [lreverse $hits] {
		if {[expr [string length $seqL] - [lindex $hit end]] < $MIN} {
			lappend len_stretch [expr [lindex $hit end] - [lindex $hit 0] + 1]
		} else {
			set seqL [string range $seqL [lindex $hit end] end]
			break
		}
	}

	set hits [regexp -inline -indices -all {A{10,}|G{10,}|T{10,}|C{10,}} $seqR]
	foreach hit $hits {
		if {[lindex $hit 0] < $MIN} {
			lappend len_stretch [expr [lindex $hit end] - [lindex $hit 0] + 1] 
		} else {
			set seqR [string range $seqR 0 [lindex $hit 0]]
			break
		}
	}	
	set hits [regexp -inline -indices -all {A{10,}|G{10,}|T{10,}|C{10,}} $seqTarget]
	foreach hit $hits {
		lappend len_stretch [expr [lindex $hit end] - [lindex $hit 0] + 1]
	}



	# Mask similar sequences with FASTA
	# ---------------------------------
	set tempL [tempfile get]
	set tempR [tempfile get]
	if {[catch "open $tempL w" fileid_outL]} {
		puts "Could not open tempfile - $fileid_outL"
		exit 1
	}
	if {[catch "open $tempR w" fileid_outR]} {
		puts "Could not open tempfile - $fileid_outR"
		exit 1
	}
	puts $fileid_outL ">SeqL \n $seqL"
	close $fileid_outL
	puts $fileid_outR ">SeqR \n $seqR"
	close $fileid_outR
	if {[catch "exec FASTAwrap.tcl $tempL $tempR" err]} {
		puts "fetching sequences failed - $err"
		exit 1
	}

	#scrape FASTAoutput
	if {[catch "exec FASTAwrap.tcl $tempL $tempR -t 1 -e $EVAL" FASTAwrap]} {
		puts "something went wrong while FASTA'ing some files - $FASTAwrap"
		exit 1
	}
	set i 0
	#puts $FASTAwrap ; #TEST
	while {$i < [llength $FASTAwrap]} {
		set FASTAdict [lindex $FASTAwrap $i]
		if {[dict get $FASTAdict complement] == 0} {
			#normal string replace
			set seqL [string replace $seqL [dict get $FASTAdict hit_start] \
				[dict get $FASTAdict hit_end] [string repeat N [dict get $FASTAdict overlap_length]]]
			set seqR [string replace $seqR [dict get $FASTAdict query_start] \
				[dict get $FASTAdict query_end] [string repeat N [dict get $FASTAdict overlap_length]]]	
		} else {
			set seqL [string replace $seqL [dict get $FASTAdict hit_start] \
				[dict get $FASTAdict hit_end] [string repeat N [dict get $FASTAdict overlap_length]]]
			set seqR [string replace $seqR [dict get $FASTAdict query_end] \
				[dict get $FASTAdict query_start]  [string repeat N [dict get $FASTAdict overlap_length]]]
		}
		incr i
	}
	
	return [list $seqTarget $seqL $seqR $len_stretch]
}

proc cg_validatesv_runPrimer3 {list_seq ID ex_primer max_prim} {
	set target [lindex $list_seq 0]
	set seqL [lindex $list_seq 1]
	set seqR [lindex $list_seq 2]
	set seq "$seqL$target$seqR"
	set targetBegin [string length $seqL]
	set targetLength [string length $target]
	
	set temp_in [tempfile get]
	if {[catch "open $temp_in w" inid]} {
		puts "Could not open tempfile - $inid"
		exit 1
	}
	#when no existing primer is known, search for 2 primers
	#otherwise just for 1.
	if {$ex_primer == 0 } {
		puts $inid "SEQUENCE_ID=$ID"
		puts $inid "SEQUENCE_TEMPLATE=$seq"
		puts $inid "SEQUENCE_TARGET=$targetBegin,$targetLength"
		puts $inid "PRIMER_NUM_RETURN=$max_prim"
		puts $inid "PRIMER_MIN_SIZE=15"
		puts $inid "PRIMER_OPT_SIZE=18"
		puts $inid "PRIMER_MAX_SIZE=23"
		puts $inid "PRIMER_OPT_TM=60"
		puts $inid "PRIMER_MAX_NS_ACCEPTED=0"
		puts $inid "PRIMER_PRODUCT_SIZE_RANGE=300-1000"
		puts $inid "P3_FILE_FLAG=0" ; # if FLAG=1, extra output files will be created
		puts $inid "=" 
		close $inid
	} else { 
		puts $inid "SEQUENCE_ID=$ID"
		puts $inid "SEQUENCE_TEMPLATE=$seq"
		puts $inid "SEQUENCE_TARGET=$targetBegin,$targetLength"
		puts $inid "PRIMER_NUM_RETURN=$max_prim"
		puts $inid "PRIMER_MIN_SIZE=15"
		puts $inid "PRIMER_OPT_SIZE=18"
		puts $inid "SEQUENCE_PRIMER=$ex_primer"
		puts $inid "PRIMER_MAX_SIZE=23"
		puts $inid "PRIMER_OPT_TM=60"
		puts $inid "PRIMER_MAX_NS_ACCEPTED=0"		
		puts $inid "PRIMER_PRODUCT_SIZE_RANGE=300-1000"
		puts $inid "P3_FILE_FLAG=0"
		puts $inid "="
		close $inid
	}

	set temp_out [tempfile get]
	if {[catch "exec primer3_core -output=$temp_out $temp_in" err]} {
		puts "something went wrong while executing Primer3 - $err"
		exit 1
	}

	#return a dict of the output of Primer3	
	set outid [open $temp_out]
	set primerOut [read $outid]
	close $outid
	set primerDict [string map {"=" " "} $primerOut]
	#puts $primerDict ; #TEST
	return $primerDict
}

proc cg_validatesv_cindex_searchgenome {db pseq {add 0}} {
	global cindex_genome maxnum
	if {![info exists cindex_genome]} {
		set cindex_genome {}
		#puts "loading genome database" ; #TEST
		foreach file [lsort -dict [glob $db/*]] {
			set file [file root $file]
			set chr [lindex [split $file -] end]
			#puts "loading chr $chr" ; #TEST
			dict set cindex_genome $chr [cindex load $file]
		}
	}
	set results {}
	set numresults 0
	dict for {chr cindex} $cindex_genome {
		set temp [cindex locate $cindex $pseq $maxnum]
		incr numresults [llength $temp]
		if {$numresults > $maxnum} {error "found $numresults"}
		if {[llength $temp]} {
			if {$add} {set temp [lmath_calc $temp + $add]}
			dict set results $chr $temp
		}
	}
	return [list $numresults $results]
}

proc cg_validatesv_runEPCR {primer1 primer2 size inv} {
	global searchGenomeDB
	global MAX_SIZE
	if {$inv != 1} {set size 0}
	set seq2 [string range $primer2 [expr [string length $primer2] - 15] end]
	set seq1 [string range $primer1 [expr [string length $primer1] - 15] end]
	if {[catch "cg_validatesv_cindex_searchgenome $searchGenomeDB $seq2 " c_output2 ]} {error "found too many hits"}
	if {[catch "cg_validatesv_cindex_searchgenome $searchGenomeDB $seq1 " c_output1 ]} {error "found too many hits"}
	set hit 0	
	set chr_col1 0
	set chr_col2 0
	while {$chr_col1 < [llength [lindex $c_output1 1]] && $chr_col2 < [llength [lindex $c_output2 1]] } {
		set chr1 [lindex [lindex $c_output1 1] $chr_col1]
		set chr2 [lindex [lindex $c_output2 1] $chr_col2]
		if {$chr1 == $chr2} {
			foreach pos1 [lindex [lindex $c_output1 1] [expr $chr_col1 + 1]] {
				foreach pos2 [lindex [lindex $c_output2 1] [expr $chr_col2 + 1]] {
					if {[::tcl::mathfunc::abs [expr $pos1 - $pos2]] < [expr $size + $MAX_SIZE] } {incr hit}
					if {$hit > 1} {error "found more then 1 amplicon"}
				}
			}
			set chr_col1 [expr $chr_col1 +2]
			set chr_col2 [expr $chr_col2 +2]
		} elseif {$chr1 > $chr2} {
			set chr_col2 [expr $chr_col2 +2]
		} else {
			set chr_col1 [expr $chr_col1 +2]	
		}
	}
	return 0
}

proc cg_validatesv_getPrimerPairs {chr patchSize patchstart size breakpointL breakpointR READSIZE EVAL MIN} {	
	# get sequences around breakpoints
	set list_seq(1) [cg_validatesv_getSeqRef $chr $patchSize $breakpointL $READSIZE]
	set list_seq(3) [cg_validatesv_getSeqRef $chr $patchSize $breakpointR $READSIZE]
	#puts "target: [lindex $list_seq1 0]"  ; #TEST

	# mask repeats in flanking sequences 
	for {set x 1} {$x <= 3} {incr x} {
		if { $x == 2 } {continue} ; #Getting seq 2 in different manner
		set list [cg_validatesv_mask_seq $list_seq($x) $EVAL $MIN]
		set len_stretch${x} [lindex $list end]
		set list_seq${x}_mask [lreplace $list end end] ; #deleting stretch info from list
	}
	if {[llength $len_stretch1] > 0 || [llength $len_stretch3] > 0 } {
		puts "WARNING: There will be mononucleotide stretches in the amplicons"
		puts "Number nucleotides in such stretch surrounding \t left breakpoint: $len_stretch1"
		puts "\t \t \t \t \t \t right breakpoint: $len_stretch3"
		#puts $list_seq1_mask ; #TEST
	}
	set list_seq2_mask [cg_validatesv_getSeqInv $list_seq1_mask $list_seq3_mask]

	#try to get the best primer pairs for the inversion breakpoints
	set primerDict1 [cg_validatesv_runPrimer3 $list_seq1_mask "leftRefSeq_${chr}_${patchstart}" "0" "40"]
	#puts $primerDict1	;#TEST
	set i 0
	while {$i < [dict get $primerDict1 PRIMER_PAIR_NUM_RETURNED]} {
		set primer1L [dict get $primerDict1 PRIMER_LEFT_${i}_SEQUENCE]
		set primer1R [dict get $primerDict1 PRIMER_RIGHT_${i}_SEQUENCE]
		#run ucsc_epcr on the 2 primers. There has to be 1 amplicon amplified
		set inv 0
		if {[catch "cg_validatesv_runEPCR $primer1L $primer1R $size $inv" out ]} {incr i; continue}
		set primerDict2 [cg_validatesv_runPrimer3 $list_seq2_mask "leftInvSeq_${chr}_${patchstart}" $primer1L "20"]
		#puts "PD2: $primerDict2" ; #TEST
		set j 0
		while {$j < [dict get $primerDict2 PRIMER_RIGHT_NUM_RETURNED]} {
			set primer2R [dict get $primerDict2 PRIMER_RIGHT_${j}_SEQUENCE]
			#run cindex_searchgenome on the 2 primers. There has to be 1 amplicon amplified
			#these are the inverted primers so the size of the inversion has to be brought into account
			set inv 1
			if {[catch "cg_validatesv_runEPCR $primer1L $primer2R $size $inv" out ]} {incr j; continue}
			set primerDict3 [cg_validatesv_runPrimer3 $list_seq3_mask "leftInvSeq_${chr}_${patchstart}" $primer2R "20"]
			#puts "PR3: $primerDict3" ; #TEST
			set k 0
			while {$k < [dict get $primerDict3 PRIMER_RIGHT_NUM_RETURNED]} {	
				set primer3R [dict get $primerDict3 PRIMER_RIGHT_${k}_SEQUENCE]
				#run epcr on the 2 primers. There has to be 1 amplicon amplified
				set inv 0
				if {[catch "cg_validatesv_runEPCR $primer2R $primer3R $size $inv" out ]} {incr k; continue}
				#making output file of all the 4 primers and there info
				set primL "primL [dict get $primerDict1 PRIMER_LEFT_${i}_SEQUENCE] [string length [dict get $primerDict1 PRIMER_LEFT_${i}_SEQUENCE]]  	
					[dict get $primerDict1 PRIMER_LEFT_${i}_TM] [dict get $primerDict1 PRIMER_LEFT_${i}_GC_PERCENT] 
					[dict get $primerDict1 PRIMER_RIGHT_${i}_SEQUENCE] [string length [dict get $primerDict1 PRIMER_RIGHT_${i}_SEQUENCE]]
					[dict get $primerDict1 PRIMER_RIGHT_${i}_TM] [dict get $primerDict1 PRIMER_RIGHT_${i}_GC_PERCENT]
					[dict get $primerDict1 PRIMER_PAIR_${i}_PRODUCT_SIZE] " 
				set primR "primR [dict get $primerDict2 PRIMER_RIGHT_${j}_SEQUENCE] [string length [dict get $primerDict2 PRIMER_RIGHT_${j}_SEQUENCE]]  	
					[dict get $primerDict2 PRIMER_RIGHT_${j}_TM] [dict get $primerDict2 PRIMER_RIGHT_${j}_GC_PERCENT] 
					[dict get $primerDict3 PRIMER_RIGHT_${k}_SEQUENCE] [string length [dict get $primerDict3 PRIMER_RIGHT_${k}_SEQUENCE]]
					[dict get $primerDict3 PRIMER_RIGHT_${k}_TM] [dict get $primerDict3 PRIMER_RIGHT_${k}_GC_PERCENT]
					[dict get $primerDict2 PRIMER_PAIR_${j}_PRODUCT_SIZE] "
				return [list $primL $primR]
				incr k
			}
			incr j
		}
		incr i
	}
	return 1
}

proc cg_validatesv_help {} {
	set help [file_read $::appdir/lib/cg_validatesv.help]
	puts [string_change $help [list @BASE@ [get ::base {[info source]}]]]
}



proc cg_validatesv args {
	# set options
	# -----------
	set READSIZE 360
	set EVAL 5
	set MIN 100
		
	foreach {key value} $args {
		switch -- $key \
			"-h" - "--help" {cg_validatesv_help ; exit 0} \
			"-f" "set file $value" \
			"-r" "set READSIZE $value" \
			"-m" "set MIN $value" \
			"-o" "set file_out $value" \
			"-e" "set EVAL $value" ;
	}
	
	if {[llength $args] < 1 || [llength $args] > 5} {
		puts "Wrong number of arguments"
		cg_validatesv_help
		exit 1
	}	
	if {[catch "open $file r" fileid]} {
		puts "Could not open input file - $fileid"
		puts "Please note that the input file is a mandatory argument."
		cg_validatesv_help
		exit 1
	}
	if {![info exists file_out]} {
		set file_out [file rootname $file]_primers.tsv
	}

	# get headers
	# -----------
	set in [gets $fileid]
	set headers [split $in "\t"]

	set column 0
	foreach head $headers {
		switch -glob -nocase $head {
			patchstart {set PATCHSTART $column}
			patchend {set PATCHEND $column}
			height {set HEIGHT $column}
			chr {set CHR $column}
		}
		incr column
	}	

	# setting header in outputfile
	# ----------------------------
	if {[catch "open $file_out w" fileid_out]} {
		puts "Could not open output file - $fileid_out"
		exit 1
	}
	set line_out "chr patchStart patchEnd primer sequenceF sizePrimerF TmF GC-contentF sequenceR sizePrimerR TmR GC-contentR sizeAmplicon"
	puts $fileid_out [join $line_out \t]

	# Getting the primerpairs for each inversion
	# ------------------------------------------
	set in [gets $fileid]
	set inversion_count 1
	while {![eof $fileid]} {
		set line [split $in "\t"]
		set patchSize [expr [lindex $line $PATCHEND] - [lindex $line $PATCHSTART]]
		set breakpointL [lindex $line $PATCHEND]
		set breakpointR [expr [lindex $line $PATCHSTART] + [lindex $line $HEIGHT]]
		set size [expr $breakpointR - $breakpointL]
		set chr [lindex $line $CHR]
		set patchstart [lindex $line $PATCHSTART]
		puts "Analyzing inversion ${chr}_${patchstart} ..."
		puts "The predicted size for this inversion is $size"

		set primerPairs [cg_validatesv_getPrimerPairs $chr $patchSize $patchstart $size $breakpointL $breakpointR $READSIZE $EVAL $MIN]
		if {$primerPairs == 1} {
			puts "No primer pairs were found for this inversion"
			tempfile clean 
			set in [gets $fileid]
			incr inversion_count
				continue
		}
		foreach primer $primerPairs {
			set line_out "$chr $patchstart [lindex $line $PATCHEND] $primer"
				puts $fileid_out [join $line_out \t]
	
		}

		tempfile clean 
		set in [gets $fileid]
		incr inversion_count
	}
	close $fileid
	return 0


}
if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_select {*}$argv
}