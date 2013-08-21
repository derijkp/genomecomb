#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral
package require dict
package require cindex
catch {package require dbi}
package require dbi_sqlite3
package require BioTcl


#######################################################
# 
#  This is a script that designs the primer pairs 
#  for sequencing breakpoints of inversions.
#
#  Written by Annelies Cassiers
# Adapted by Peter De Rijk
#
#######################################################

set Temp 100
foreach rtype {
	partrepeat,low repeat,low 3repeat,low completerepeat,low
	partrepeat,multi repeat,multi 3repeat,multi completerepeat,multi
} {
	set rscore($rtype) $Temp
	incr Temp 100
}

proc makeprimers_primer3_core {seq rstart rend args} {
	set size 600
	set temperature 60
	set numreturn 5
	foreach {key value} $args {
		switch $key {
			-size {set size $value}
			-temperature {set temperature $value}
			-included_region {set included $value}
			-num_return {set numreturn $value}
			-task {
				if {![inlist {generic check_primers pick_primer_list pick_sequencing_primers pick_cloning_primers pick_discriminative_primers pick_pcr_primers pick_pcr_primers_and_hyb_probe pick_left_only pick_right_only pick_hyb_probe_only} $value]} {
					error "Unkown task $value"
				}
				set task $value
			}
			default {
				error "Unkown option $key"
			}
		}
	}
	global primer3 cachedir
	if {![info exists cachedir] || [file isdir $cachedir]} {
		set tempdir [tempdir]
	} else {
		set tempdir $cachedir
	}
puts  $tempdir/primer3input.txt
	set f [open $tempdir/primer3input.txt w]
	puts $f "SEQUENCE_ID=temp"
	puts $f "SEQUENCE_TEMPLATE=$seq"
	puts $f "PRIMER_MIN_SIZE=18"
	puts $f "PRIMER_OPT_SIZE=22"
	puts $f "PRIMER_MAX_SIZE=26"
	puts $f "PRIMER_MIN_TM=[expr {$temperature-1}]"
	puts $f "PRIMER_OPT_TM=$temperature"
	puts $f "PRIMER_MAX_TM=[expr {$temperature+1}]"
	puts $f "PRIMER_MAX_SELF_ANY=6"
	puts $f "PRIMER_MAX_END_STABILITY=9"
	puts $f "PRIMER_MAX_NS_ACCEPTED=0"
	puts $f "PRIMER_MAX_SELF_END=2"
	puts $f "PRIMER_MAX_POLY_X=3"
	puts $f "PRIMER_MIN_GC=20"
	puts $f "P3_FILE_FLAG=1"
	puts $f "PRIMER_LOWERCASE_MASKING=1"
	puts $f "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0"
	if {[info exists task]} {
		puts $f "PRIMER_TASK=$task"
	}
	if {[info exists included]} {
		puts $f "SEQUENCE_INCLUDED_REGION=$included"
	}
#	puts $f "PRIMER_TM_SANTALUCIA=1"
#	puts $f "PRIMER_SALT_CORRECTIONS=1"
#	puts $f "PRIMER_SALT_CONC="
#	puts $f "PRIMER_DIVALENT_CONC="
#	puts $f "PRIMER_DNTP_CONC="
#	puts $f "PRIMER_DNA_CONC=50.0"
	puts $f "PRIMER_NUM_RETURN=$numreturn"
	puts $f "PRIMER_PRODUCT_SIZE_RANGE=75-$size"
#	puts $f "PRIMER_THERMODYNAMIC_ALIGNMENT=1"
	puts $f "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$::externdir/primer3_config/"
	puts $f "SEQUENCE_TARGET=$rstart,[expr {$rend-$rstart}]"
	puts $f "="
	close $f
	set keep [pwd]
	cd $tempdir
	set result [exec primer3_core -strict_tags < $tempdir/primer3input.txt]
	set table {}
	foreach {key value} [split $result =\n] {
		switch -regexp $key {
			{^PRIMER_LEFT_?[0-9]*$} {set left [lindex [split $value ,] 0]}
			{^PRIMER_RIGHT_?[0-9]*$} {set right [lindex [split $value ,] 0]}
			{^PRIMER_LEFT_?[0-9]*_SEQUENCE$} {set leftseq $value}
			{^PRIMER_RIGHT_?[0-9]*_SEQUENCE$} {set rightseq $value}
			{^PRIMER_PRODUCT_SIZE_?[0-9]*$} {
				lappend table [list $left [expr {$right-[string length $rightseq]+1}] $leftseq $rightseq]
			}
		}
	}
	set left [lrange [split [string trim [file_read temp.for]] \n] 3 end]
	set right [lrange [split [string trim [file_read temp.rev]] \n] 3 end]
	set primerlist {}
	set result {}
	set min [expr {$rend - $size}]
	set max [expr {$rstart + $size}]
	foreach list [list $left $right] strand {+ -} {
		set temp {}
		foreach line $list {
			foreach {num seq start len} $line break
			if {$strand eq "-"} {
				set start [expr {$start-$len+1}]
			}
			lset line 2 $start
			set end [expr {$start + $len - 1}]
			if {$start < $min} continue
			if {$end > $max} continue
			set line [linsert $line 3 $end $strand]
			set line [linsert $line 6 0]
			lappend line ? ? clean {}
			lappend temp $line
		}
		lappend result $temp
	}
	cd $keep
	return $result
}

proc cg_validatesv_getSeqRef {dbdir chr patch breakpoint1 breakpoint2 read MIN} {
	set beginTarget [expr $breakpoint1 - 200]
	set endTarget [expr $breakpoint1 + $read - $patch + 250]

	set beginL [expr $beginTarget - 400]
	set endL [expr $beginTarget - 1] ;#otherwise there will be overlap

	set beginR [expr $endTarget + 1] ;#otherwise there will be overlap
	set endR [expr $endTarget + 400]

	#check that primers will fall into inversion
	if {$breakpoint1 < $breakpoint2} {
		#we're dealing with left breakpoint
		if {[expr $endTarget + $MIN] > $breakpoint2  } {
			set endTarget [expr $breakpoint2 - $MIN]
			set beginR [expr $endTarget + 1]
			set endR $breakpoint2
		} elseif {$endR > $breakpoint2 } {
			set endR $breakpoint2
		}

	} else {
		#dealing with right breakpoint
		if {[expr $beginTarget - $MIN] < $breakpoint2 } {
			set beginTarget [expr $breakpoint2 + $MIN]
			set beginL $breakpoint2
			set endL [expr $beginTarget - 1]
		} elseif {$beginL < $breakpoint2 } {
			set beginL $breakpoint2
		}
	}
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	puts "getrefseq $beginL $endL $beginTarget $endTarget $beginR $endR"
	set seqTarget [genome_get $fg $chr $beginTarget $endTarget]
	set seqL [genome_get $fg $chr $beginL $endL]
	set seqR [genome_get $fg $chr $beginR $endR]
	genome_close $fg
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

proc parseFASTA {INPUT} {
##################################################
# This procedure is made to parse FASTA output.
#
# written by Bart Aelterman 2010
##################################################
	regsub -all {<<<} $INPUT @ INPUT
	set INPUT [lindex [split $INPUT @] 0]
	regsub -all {>>>} $INPUT @ INPUT
	set INPUT_LIST [lrange [split $INPUT @] 1 end]
	set OUTLIST {}
	foreach {SUMMARY HITS_INFO} $INPUT_LIST {
		regexp {^(.*?),} $HITS_INFO MATCH QUERY_ID
		regsub -all {>>} $HITS_INFO @ HITS_INFO
		set ALIGNMENT_INFOS [lrange [split $HITS_INFO @] 1 end]
		foreach ALIGNMENT $ALIGNMENT_INFOS {
			regexp {^(.*?)\n} $ALIGNMENT MATCH HIT_ID 
			regsub -all {>} $ALIGNMENT @ ALIGNMENT
			set AL_LIST [split $ALIGNMENT @]
			foreach {AL_GENERAL QUERY_INFO HIT_INFO} $AL_LIST {break}
			if {![regexp {; fa_frame: (.*?)\n} $AL_GENERAL MATCH F_OR_R]} {
				puts "no F_OR_R variable found in line: $AL_GENERAL"
				return "no F_OR_R variable found in line: $AL_GENERAL"
			}
			if {$F_OR_R eq f} {
				set COMPLEMENT 0
			} else {
				set COMPLEMENT 1
			}
			#
			regexp {fa_initn: (.*?)\n} $AL_GENERAL MATCH INITN
			regexp {fa_init1: (.*?)\n} $AL_GENERAL MATCH INIT1
			regexp {fa_opt: (.*?)\n} $AL_GENERAL MATCH OPT
			regexp {fa_z-score: (.*?)\n} $AL_GENERAL MATCH ZSCORE
			regexp {fa_bits: (.*?)\n} $AL_GENERAL MATCH BITS
			regexp {fa_expect: (.*?)\n} $AL_GENERAL MATCH EVALUE
			regexp {[bs][sw]_ident: (.*?)\n} $AL_GENERAL MATCH IDENTITY
			if {![regexp {bs_sim: (.*?)\n} $AL_GENERAL MATCH SIMILAR]} {
				regexp {sw_gident: (.*?)\n} $AL_GENERAL MATCH SIMILAR
			}
			regexp {[bs][sw]_overlap: (.*?)\n} $AL_GENERAL MATCH OVERLAP
			#
			regexp {al_start: (.*?)\n} $QUERY_INFO MATCH QUERY_START
			regexp {al_stop: (.*?)\n} $QUERY_INFO MATCH QUERY_END
			regexp {al_display_start: .*?\n(.*?)$} $QUERY_INFO MATCH QUERYSEQ
			regsub -all {\n} $QUERYSEQ {} QUERYSEQ
			#
			regexp {al_start: (.*?)\n} $HIT_INFO MATCH HIT_START
			regexp {al_stop: (.*?)\n} $HIT_INFO MATCH HIT_END
			regexp {al_display_start: .*?\n(.*?)[0-9]*$} $HIT_INFO MATCH HITSEQ
			regsub -all {\n} $HITSEQ {} HITSEQ
			#
			set ALL_VALUES [list $INITN $INIT1 $OPT $ZSCORE $BITS $EVALUE $IDENTITY $SIMILAR $OVERLAP $QUERY_START $QUERY_END $QUERYSEQ $HIT_START $HIT_END $HITSEQ]
			set NEW_VALUES {}
			foreach VALUE $ALL_VALUES {
				set NEW_VALUE [string trim $VALUE]
				lappend NEW_VALUES $NEW_VALUE
			}
			foreach {INITN INIT1 OPT ZSCORE BITS EVALUE IDENTITY SIMILAR OVERLAP QUERY_START QUERY_END QUERYSEQ HIT_START HIT_END HITSEQ} $NEW_VALUES {break}
			#
			set OUTDICT [dict create query $QUERY_ID hit \
				$HIT_ID initn $INITN init1 $INIT1 opt $OPT \
				zscore $ZSCORE \ bits $BITS evalue $EVALUE \
				identity $IDENTITY similarity $SIMILAR \
				overlap_length $OVERLAP \
				hit_start $HIT_START hit_end $HIT_END query_start $QUERY_START \
				query_end $QUERY_END complement $COMPLEMENT hit_seq $HITSEQ \
				query_seq $QUERYSEQ]
			lappend OUTLIST $OUTDICT
		}
	}
	return $OUTLIST
}

package require dict

proc FASTAwrap {query database args} {
#################################################
# Program to automate FASTA runs
# This program was based on run_FASTA.tcl but
# this one is more flexible: it allows users to
# give options to alter alignment settings.
#
# written by Bart Aelterman 2010
#################################################
	# set variables and defaults
	set TMPPATH [tempdir]
	# args
	set ECUTOFF 5
	set BVAL 20
	set DVAL 20
	foreach {OPT VALUE} $args {
		switch -- $OPT {
			-b {set BVAL $VALUE}
			-d {set DVAL $VALUE}
			-e {set ECUTOFF $VALUE}
		}
	}
	# query and database filenames cannot be longer than 119 characters
	if {[string length $query] > 119} {
		puts "filename of query is too long. Copying query to ${TMPPATH}/query.fa"
		catch {file copy -force $query ${TMPPATH}/query.fa} result
		set query "${TMPPATH}/query.fa"
	}
	if {[string length $database] > 119} {
		puts "filename of database is too long. Copying database to ${TMPPATH}/database.fa"
		catch {file copy -force $database ${TMPPATH}/database.fa} result
		set database "${TMPPATH}/database.fa"
	}
	if {[catch {
		exec fasta34 -n -b $BVAL -d $DVAL -E $ECUTOFF -m 10 -H $query $database
	} result]} {
		error "error with fasta search occured:\n$result"
	}
	set OUTLIST [parseFASTA $result]
	return $OUTLIST
}

proc cg_validatesv_runFasta {list_seq_src EVAL} {
	global log ; #TEST
	set seqTarget [lindex $list_seq_src 0]
	set seqL [lindex $list_seq_src 1]
	set seqR [lindex $list_seq_src 2]
	puts $log "$seqL$seqTarget$seqR" ; #TEST
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
	puts $fileid_outL ">SeqL\n$seqL"
	close $fileid_outL
	puts $fileid_outR ">SeqR\n$seqR"
	close $fileid_outR
	#run FASTA - output = tclDict
	set FASTAwrap [FASTAwrap $tempL $tempR -e $EVAL]
	set i 0
	#puts $log $FASTAwrap ; #TEST
	#making distinction between masked and non-masked sequences
	set seqL_mask $seqL
	set seqR_mask $seqR
	while {$i < [llength $FASTAwrap]} {
		set FASTAdict [lindex $FASTAwrap $i]
		if {[dict get $FASTAdict complement] == 0} {
			#normal string replace
			set seqL_mask [string replace $seqL_mask [expr [dict get $FASTAdict query_start] - 1] \
				[expr [dict get $FASTAdict query_end] - 1] [string repeat N [expr [dict get $FASTAdict query_end] - [dict get $FASTAdict query_start] + 1 ]]]
			set seqR_mask [string replace $seqR_mask [expr [dict get $FASTAdict hit_start] - 1] \
				[dict get $FASTAdict hit_end] [string repeat N [expr [dict get $FASTAdict hit_end] - [dict get $FASTAdict hit_start] + 1 ]]]
		} else {
			set seqL_mask [string replace $seqL_mask [ expr [dict get $FASTAdict query_end] - 1] \
				[ expr [dict get $FASTAdict query_start] - 1] [string repeat N [expr [dict get $FASTAdict query_start] - [dict get $FASTAdict query_end] + 1 ]]]
			set seqR_mask [string replace $seqR_mask [ expr [dict get $FASTAdict hit_start] - 1] \
				[ expr [dict get $FASTAdict hit_end] - 1]  [string repeat N [expr [dict get $FASTAdict hit_end] - [dict get $FASTAdict hit_start] + 1 ]]]
		}
		incr i
	}
	return [list $seqTarget $seqL_mask $seqR_mask]
}

proc cg_validatesv_searchStretches {list_seq list_seq_mask MIN} {
	global prim_min
	global log ; #TEST
	set seqTarget [lindex $list_seq 0]
	set seqL [lindex $list_seq 1]
	set seqR [lindex $list_seq 2]

	set seqTarget_mask [lindex $list_seq_mask 0]
	set seqL_mask [lindex $list_seq_mask 1]
	set seqR_mask [lindex $list_seq_mask 2]

	# Search for homopolymer stretches
	# --------------------------------
	
	set len_stretch {} ; #so user knows how much nucleotides are located in stretches
	set hits [regexp -inline -indices -all {A{10,}|G{10,}|T{10,}|C{10,}} $seqL]
	#try to avoid the stretches in the sequences were the primers will be found
	puts $log "hits: $hits" ; #TEST
	if {![llength $hits]} {set found 1}

	foreach hit [lreverse $hits] {
		puts $log "aantal N: [regexp -all {N} [string range $seqL_mask [lindex $hit end] [string length $seqL]]]" ; #TEST
		puts $log "aantal normale stretches groter dan 50bp: [regexp -all {[AGTC]{[expr $MIN/2],}} [string range $seqL_mask [lindex $hit end] [string length $seqL]]]" ; #TEST
		if {[expr [string length $seqL] - [lindex $hit end] ] < $MIN && \
			[regexp -all {[AGTC]{$prim_min,}} [string range $seqL_mask [lindex $hit end] [string length $seqL]]] == 0 } {
			lappend len_stretch [expr [lindex $hit end] - [lindex $hit 0] + 1]
		} else {
			set seqL_mask [string range $seqL_mask [lindex $hit end] end]
			set found 1
			break
		}
	}

	set hits [regexp -inline -indices -all {A{10,}|G{10,}|T{10,}|C{10,}} $seqR]
	foreach hit $hits {
		if {[lindex $hit 0] < $MIN && \
			[regexp -all {[AGTC]{$prim_min,}} [string range $seqL_mask 0 [lindex $hit 0]]] == 0	} {
			lappend len_stretch [expr [lindex $hit end] - [lindex $hit 0] + 1] 
		} else {
			set seqR_mask [string range $seqR_mask 0 [lindex $hit 0]]
			set found 1
			break
		}
	}

	if {![info exists found]} {
		puts "WARNING: When there are no primers found, the cause is probably the similarity of sequences adjacent to the breakpoint"
	}
	unset -nocomplain found 
	
	set hits [regexp -inline -indices -all {A{10,}|G{10,}|T{10,}|C{10,}} $seqTarget]
	foreach hit $hits {
		lappend len_stretch [expr [lindex $hit end] - [lindex $hit 0] + 1]
	}

	puts $log "$seqL_mask$seqTarget$seqR_mask" ; #TEST

	
	return [list $seqTarget $seqL_mask $seqR_mask $len_stretch]

}


proc cg_validatesv_mask_seq {list_seq_src EVAL MIN fasta} {
	if {$fasta == 1} {
		set list_seq_mask [cg_validatesv_runFasta $list_seq_src $EVAL]
	} else {
		set list_seq_mask $list_seq_src
	}
	set list_seq_new [cg_validatesv_searchStretches $list_seq_src $list_seq_mask $MIN]
	return $list_seq_new
}

proc cg_validatesv_runPrimer3 {list_seq_mask ID ex_primer max_prim} {
	global max_amplicon externdir
	set target [lindex $list_seq_mask 0]
	set seqL [lindex $list_seq_mask 1]
	set seqR [lindex $list_seq_mask 2]
	set seq "$seqL$target$seqR"
	set targetBegin [string length $seqL]
	set targetLength [string length $target]
	#
	set temp_in [tempfile get]
	set temp_out [tempfile get]
	#
	if {[catch "open $temp_in w" inid]} {
		error "Could not open tempfile - $inid"
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
		puts $inid "PRIMER_PRODUCT_SIZE_RANGE=300-$max_amplicon"
		puts $inid "PRIMER_THERMODYNAMIC_ALIGNMENT=1"
		puts $inid "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$externdir/primer3_config/"
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
		puts $inid "PRIMER_PRODUCT_SIZE_RANGE=300-$max_amplicon"
		puts $inid "PRIMER_THERMODYNAMIC_ALIGNMENT=1"
		puts $inid "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$externdir/primer3_config/"
		puts $inid "P3_FILE_FLAG=0"
		puts $inid "="
		close $inid
	}
	if {[catch {exec primer3_core < $temp_in > $temp_out} err]} {
		error "something went wrong while executing Primer3 - $err"
	}
	#return a dict of the output of Primer3	
	set outid [open $temp_out]
	set primerOut [read $outid]
	close $outid
	set primerDict [string map {"=" " "} $primerOut]
	#puts $primerDict ; #TEST
	return $primerDict
}

proc cg_validatesv_cindex_searchgenome {DB pseq {add 0}} {
	global cindex_genome maxnum
	global log ; #TEST
	if {![info exists cindex_genome]} {
		set cindex_genome {}
		#puts "loading genome database" ; #TEST
		foreach file [ssort -natural [glob $DB/*]] {
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
		puts $log "numresult: $numresults" ;#TEST
		if {$numresults > $maxnum} {error "found $numresults"}
		if {[llength $temp]} {
			if {$add} {set temp [lmath_calc $temp + $add]}
			dict set results $chr $temp
		}
	}
	puts $log "end cindex" ;#TEST
	return [list $numresults $results]
}



proc cg_validatesv_searchAmplicon {seq1 seq2 size} {
	global searchGenomeDB
	global MAX_SIZE
	if {[catch {cg_validatesv_cindex_searchgenome $searchGenomeDB $seq2} c_output2 ]} {error "found too many hits"}
	if {[catch {cg_validatesv_cindex_searchgenome $searchGenomeDB $seq1} c_output1 ]} {error "found too many hits"}
	if {[lindex $c_output1 0] == 0} {return 0}
	if {[lindex $c_output2 0] == 0} {return 0}
	set hit 0	
	set chr_col1 0
	set chr_col2 0
	set cchr 0
	set cbegin1 0
	set cend1 0
	set cbegin2 0
	set cend2 0
	while {$chr_col1 < [llength [lindex $c_output1 1]] && $chr_col2 < [llength [lindex $c_output2 1]] } {
		set chr1 [lindex [lindex $c_output1 1] $chr_col1]
		set chr2 [lindex [lindex $c_output2 1] $chr_col2]
		if {$chr1 == $chr2} {
			foreach pos1 [lindex [lindex $c_output1 1] [expr $chr_col1 + 1]] {
				foreach pos2 [lindex [lindex $c_output2 1] [expr $chr_col2 + 1]] {
					if {[::tcl::mathfunc::abs [expr $pos1 - $pos2]] < [expr $size + $MAX_SIZE] } {
						incr hit 
						if {$hit > 1} {error "found more then 1 amplicon"}
						set cchr $chr1
						set cbegin1 $pos1
						set cend1 [expr $pos1 + 15]
						set cbegin2 $pos2
						set cend2 [expr $pos2 + 15]
					}
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
	
	return [list $hit $cchr $cbegin1 $cend1 $cbegin2 $cend2]


}

proc cg_validatesv_getRepScore {fts len start} {
	global rscore 
	set score 0
	list_foreach {fstart fend name} $fts {
		set repeatname [lindex $name 0 0]
		if {$repeatname eq "trf"} {
			set repeattype low
		} else {
			set repeattype multi
		}
		set cfstart [expr $fstart-$start]
		if {$cfstart < 0} {set cfstart 0}
		set cfend [expr $fend-$start]
		if {$cfend >= $len} {set cfend [expr {$len-1}]}
		set pct [expr {double($cfend+$cfstart)/$len}]
		if {($pct > 0.8) && ($repeatname eq "dust") || ([string range $repeatname 0 2] eq "Alu")} {
			set repeatloc completerepeat
			set score 100000
		} else {
			if {$pct > 0.9} {
				set repeatloc completerepeat
			} elseif {($cfend >= [expr {$len-1}]) && ($cfstart <= [expr {$len-6}])} {
				set repeatloc 3repeat
			} elseif {[expr {double($len - ($cfend - $cfstart))/$len}] > 0.6} {
				set repeatloc partrepeat
			} else {
				set repeatloc repeat
			}
			if {$rscore($repeatloc,$repeattype) > $score} {
				set score $rscore($repeatloc,$repeattype)
			}
		}
	}
	


	#puts "score: $score" ; #TEST
	return $score
}

proc log args {
	#puts "log: " ; #TEST
	#puts $args ; #TEST
	#puts  "" ; #TEST
}

proc cg_validatesv_repeatSearch {cchr cstart cend } {
	global archive
	global extraseq 
	global cachedir a
	unset -nocomplain a

	# get ensembl data for region
	# ---------------------------
	set emblfile $cachedir/${cchr}_${cstart}.embl
	if {![file exists $emblfile]} {
		set embl [ensembl_getregion $cchr [expr {$cstart-$extraseq}] [expr {$cend+$extraseq}] -archive $archive]
		file_write $emblfile $embl
	}
	catch {e destroy}
	catch {rename e {}}
	EmblFile new e
	e open $emblfile

	# store repeats from data in database
	# -----------------------------------
	catch {db destroy}
	dbi_sqlite3 db
	db open :memory:
	db exec {create table ft (id integer primary key, type text, chromosome text, start integer, end integer, name text)}
	set fts [e features 0]
	set id 1
	foreach ft $fts {
		foreach {type loc descr} $ft break
		set floc [lindex $loc 0]
		set start [dict get $floc start]
		if {[dict exists $floc complement]} {set complement 1} else {set complement 0}
		set end [dict get $floc end]
		if {$end < $start} continue
		if {[dict exists $floc acc]} continue
		set filter 0
		if {$type eq "repeat_region"} {
			db set [list ft $id] type repeat name [lindex [dict get [lindex $ft end] note] 0] chromosome $cchr start $start end $end
			incr id
		}
	}
	#puts [join [db exec {select * from ft}] \n]; #TEST
	
	# check whether primers falls into such repeat
	# --------------------------------------------
	set primstart $extraseq
	set primend [expr $extraseq + ($cend - $cstart)]
	set fts [db exec {
		select start,end,name from ft
		where start <= ? and end >= ? and type = 'repeat'
	} $primend $primstart]
	#puts $fts ; #TEST
	return $fts

}

proc cg_validatesv_runEPCR {primer1 primer2 prod_size size inv} {
	global PrimerPair 
	global log ; #TEST
	set primer1 [string toupper $primer1]
	set primer2 [string toupper $primer2]
 	# if it is searching for the inverted primers,
	# the reverse primer will be located at the other breakpoint
	# so you have to incorporate the size of the inversion
	if {$inv != 1} {
		set size 0
		set seq2 [seq_complement [string range [lindex $primer2 0] 0 14]]
	} else {
		set seq2 [string range [lindex $primer2 0] 0 14]
	}
	set seq1 [string range [lindex $primer1 0] [expr [string length [lindex $primer1 0]] - 15] end]
	set len1 [string length [lindex $primer1 0]]
	set len2 [string length [lindex $primer2 0]]
	if {[catch {cg_validatesv_searchAmplicon $seq1 $seq2 $size} hits]} {error $hits}
	if {[lindex $hits 0] < 1} {error "Found no amplicon"}

	#getting position of hit 
	set cchr [lindex $hits 1]
	set cbegin1 [lindex $hits 2]
	set cend1 [lindex $hits 3]
	set cbegin2 [lindex $hits 4]
	set cend2 [lindex $hits 5]

	#check the complements
	set seq1 [seq_complement $seq1]
	set seq2 [seq_complement $seq2]
	if {[catch "cg_validatesv_searchAmplicon $seq1 $seq2 $size" hits]} {error $hits}
	if {[lindex $hits 0] == 1} {error "Found amplicon with primer complements"}

	#check for repeats
	if {$inv != 1} {
		set prim_fts1 [cg_validatesv_repeatSearch $cchr $cbegin1 $cend1]
		set prim_fts2 [cg_validatesv_repeatSearch $cchr $cbegin2 $cend2] 
		#puts "fts1: $prim_fts1" ; #TEST
		if {[llength $prim_fts1] && [llength $prim_fts2]} {
			#if both primers fall into repeat, add both scores together 
			#and compare that score to the score's of the other primer pairs
			set rep_score1 [cg_validatesv_getRepScore $prim_fts1 $len1 $cbegin1]
			set rep_score2 [cg_validatesv_getRepScore $prim_fts2 $len2 $cbegin2]
			set rep_score [expr $rep_score2	+ $rep_score1]
			if {[info exists PrimerPair]} {
				if {[lindex $PrimerPair end] > $rep_score } {
					set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size $rep_score]
				}
			} else { 
				set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size $rep_score]
			}
			puts $log "PP: $PrimerPair" ; #TEST
			return "repeat"
		} elseif {[llength $prim_fts1]} {
			set rep_score [cg_validatesv_getRepScore $prim_fts1 $len1 $cbegin1]
			if {[info exists PrimerPair]} {
				if {[lindex $PrimerPair end] > $rep_score } {
					set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size $rep_score]
				}
			} else { 
				set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size $rep_score]
			}
			puts $log "PP: $PrimerPair" ; #TEST
			return "repeat"
		} elseif {[llength $prim_fts2]} {
			set rep_score [cg_validatesv_getRepScore $prim_fts2 $len2 $cbegin2]
			if {[info exists PrimerPair]} {
				if {[lindex $PrimerPair end] > $rep_score } {
					set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size $rep_score]
				}
			} else { 
				set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size $rep_score]
			}
			puts $log "PP: $PrimerPair" ; #TEST
			return "repeat"
		}
		set PrimerPair [list $primer1 $cbegin1 $primer2 $cbegin2 $prod_size]	
		puts $log "PP_noRepeat: $PrimerPair" ; #TEST
	}

	return 0
}

proc cg_validatesv_getSeq {dbdir chr breakpointL breakpointR} {
	#getting sequences for small inversions
	set beginTarget [expr $breakpointL - 200]
	set endTarget [expr $breakpointR + 200]
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	set seqTarget [genome_get $fg $chr $beginTarget $endTarget]
	set beginL [expr $beginTarget - 100]
	set endL [expr $beginTarget - 1] ;#otherwise there will be overlap
	set seqL [genome_get $fg $chr $beginL $endL]
	set beginR [expr $endTarget + 1] ;#otherwise there will be overlap
	set endR [expr $endTarget + 200]
	set seqR [genome_get $fg $chr $beginR $endR]
	genome_close $fg
	return [list $seqTarget $seqL $seqR]
}

proc cg_validatesv_getOnePair {dbdir chr patchstart breakpointL breakpointR size EVAL MIN fasta} {
	global log ; #TEST
	global PrimerPair
	unset -nocomplain PrimerPair

	puts $log "sequencing over the whole inversion...." ; #TEST
	#get sequences over breakpoints
	set list_seq_src [cg_validatesv_getSeq $dbdir $chr $breakpointL $breakpointR]
	
	# mask repeats in flanking sequences 
	set list [cg_validatesv_mask_seq $list_seq_src $EVAL $MIN $fasta]
	set len_stretch [lindex $list end]
	set list_seq_mask [lreplace $list end end] ; #deleting stretch info from list

	if {[llength $len_stretch] > 0 } {
		puts "WARNING: There will be mononucleotide stretches in the amplicons"
		puts "Number nucleotides in such stretch: $len_stretch"
		#puts $list_seq_mask ; #TEST
	}

	#try to get the best primer pairs for the inversion breakpoints
	set primerDict [cg_validatesv_runPrimer3 $list_seq_mask "leftRefSeq_${chr}_${patchstart}" 0 40]
	puts $log "PD: $primerDict"	;#TEST
	set i 0
	while {$i < [dict get $primerDict PRIMER_PAIR_NUM_RETURNED]} {
		set primerF [list [dict get $primerDict PRIMER_LEFT_${i}_SEQUENCE]	\
			[dict get $primerDict PRIMER_LEFT_${i}_TM]  [dict get $primerDict PRIMER_LEFT_${i}_GC_PERCENT] ]
		set primerR [list [dict get $primerDict PRIMER_RIGHT_${i}_SEQUENCE] \
			[dict get $primerDict PRIMER_RIGHT_${i}_TM]  [dict get $primerDict PRIMER_RIGHT_${i}_GC_PERCENT] ]
		set prod_size [dict get $primerDict PRIMER_PAIR_${i}_PRODUCT_SIZE]
		#run ucsc_epcr on the 2 primers. There has to be 1 amplicon amplified
		set inv 0
		if {[catch {cg_validatesv_runEPCR $primerF $primerR $prod_size $size $inv} out]} {
			#puts "1: $out"; #TEST
			incr i; continue
		}
		#picking the primer with no repeat
		if {$out != "repeat"} {set PP 1 ;  break}
		incr i
	}
		
	if {[info exists PrimerPair]} {
		#picking the best primer
		set primerF [lindex $PrimerPair 0]
		set beginF [lindex $PrimerPair 1]
		set primerR [lindex $PrimerPair 2]
		set beginR [lindex $PrimerPair 3]
		set prod_size [lindex $PrimerPair 4]
		set label "inv_${chr}_[string range $patchstart 0 [expr [string length $patchstart] - 7]]"	
		#making output file of all the 2 primers and there info
		set primF "${label}_F [lindex $primerF 0] {} {} {} {pg} {1} {Hs} $chr {} {} {} $beginF [lindex $primerF 1]	{}"
		set primR	"${label}_R [lindex $primerR 0] {} {} {} {pg} {1} {Hs} $chr {} {} {} $beginR [lindex $primerR 1] {} $prod_size "
		if {[llength $PrimerPair] > 5} {set repeat 1 }
		if {[info exists repeat] } {
			puts "WARNING: Primerpair is located in a repeat"
		}
		return [list $primF $primR]
	}
	
	return 1


}

proc cg_validatesv_getPrimerPairs {dbdir chr patchSize patchstart size breakpointL breakpointR READSIZE EVAL MIN fasta} {	
	global log ; #TEST
	global PrimerPair
	unset -nocomplain PrimerPair
	# get sequences around breakpoints
	set list_seq(1) [cg_validatesv_getSeqRef $dbdir $chr $patchSize $breakpointL $breakpointR $READSIZE $MIN]
	set list_seq(3) [cg_validatesv_getSeqRef $dbdir $chr $patchSize $breakpointR $breakpointL $READSIZE $MIN]
	#puts "target: [lindex $list_seq1 0]"  ; #TEST

	#puts $log "pre: [lindex $list_seq(3) 0][lindex $list_seq(3) 1][lindex $list_seq(3) 2]" ; #TEST
	
	# mask repeats in flanking sequences 
	for {set x 1} {$x <= 3} {incr x} {
		if { $x == 2 } {continue} ; #Getting seq 2 in different manner
		set list [cg_validatesv_mask_seq $list_seq($x) $EVAL $MIN $fasta]
		set len_stretch${x} [lindex $list end]
		set list_seq${x}_mask [lreplace $list end end] ; #deleting stretch info from list
	}
	#puts $log "post1: [lindex $list_seq1_mask 0] [lindex $list_seq1_mask 1] [lindex $list_seq1_mask 2]" ; #TEST
	#puts $log "post3: [lindex $list_seq3_mask 0] [lindex $list_seq3_mask 1] [lindex $list_seq3_mask 2]" ; #TEST

	if {[llength $len_stretch1] > 0 || [llength $len_stretch3] > 0 } {
		puts "WARNING: There will be mononucleotide stretches in the amplicons"
		puts "Number nucleotides in such stretch surrounding \t left breakpoint: $len_stretch1"
		puts "\t \t \t \t \t \t right breakpoint: $len_stretch3"
		#puts $list_seq1_mask ; #TEST
	}
	# make seq2 from seq1 and seq3
	set list_seq2_mask [cg_validatesv_getSeqInv $list_seq1_mask $list_seq3_mask]

	#puts $log "post2: [lindex $list_seq2_mask 0] [lindex $list_seq2_mask 1] [lindex $list_seq2_mask 2]" ; #TEST

	#try to get the best primer pairs for the inversion breakpoints
	set primerDict1 [cg_validatesv_runPrimer3 $list_seq1_mask "leftRefSeq_${chr}_${patchstart}" 0 40]
	puts $log "PD1: $primerDict1"	;#TEST
	set len [dict get $primerDict1 PRIMER_PAIR_NUM_RETURNED]
	set i 0
	while {$i < $len} {
		set primerLF [list [dict get $primerDict1 PRIMER_LEFT_${i}_SEQUENCE]	\
			[dict get $primerDict1 PRIMER_LEFT_${i}_TM]  [dict get $primerDict1 PRIMER_LEFT_${i}_GC_PERCENT] ]
		set primerLR [list [dict get $primerDict1 PRIMER_RIGHT_${i}_SEQUENCE] \
			[dict get $primerDict1 PRIMER_RIGHT_${i}_TM]  [dict get $primerDict1 PRIMER_RIGHT_${i}_GC_PERCENT] ]
		set prod_sizeL [dict get $primerDict1 PRIMER_PAIR_${i}_PRODUCT_SIZE]
		#run ucsc_epcr on the 2 primers. There has to be 1 amplicon amplified
		set inv 0
		if {[catch {cg_validatesv_runEPCR $primerLF $primerLR $prod_sizeL $size $inv} out]} {
			puts $log "out1: $out"; #TEST
			incr i; continue
		}
		#picking the primer with no repeat
		if {$out != "repeat"} {set PP 1 ; break }
		incr i
	}
	#picking the best primer
	if {[info exists PrimerPair]} {
		set primerLF [lindex $PrimerPair 0]
		set beginLF [lindex $PrimerPair 1]
		set primerLR [lindex $PrimerPair 2]
		set beginLR [lindex $PrimerPair 3]
		set prod_sizeL [lindex $PrimerPair 4]	
	} else { 
		return 1	
	}
	if {[llength $PrimerPair] > 5} {set repeatL 1} ; #for later determination of location of repeats
	unset -nocomplain PrimerPair
	unset -nocomplain PP
	set primerDict2 [cg_validatesv_runPrimer3 $list_seq2_mask "leftInvSeq_${chr}_${patchstart}" [lindex $primerLF 0] "20"]
	puts $log "PD2: $primerDict2" ; #TEST
	set j 0
	while {$j < [dict get $primerDict2 PRIMER_RIGHT_NUM_RETURNED]} {
		set primerRF [list [dict get $primerDict2 PRIMER_RIGHT_${j}_SEQUENCE] \
			[dict get $primerDict2 PRIMER_RIGHT_${j}_TM] [dict get $primerDict2 PRIMER_RIGHT_${j}_GC_PERCENT] ]
		#run cindex_searchgenome on the 2 primers. There has to be 1 amplicon amplified
		#these are the inverted primers so the size of the inversion has to be brought into account
		set prod_size 0
		set inv 1
		if {[catch "cg_validatesv_runEPCR {$primerLF} {$primerRF} $prod_size $size $inv" out ]} {
			puts $log "out2: $out"
			incr j; continue
		}
		set primerDict3 [cg_validatesv_runPrimer3 $list_seq3_mask "leftInvSeq_${chr}_${patchstart}" [lindex $primerRF 0] "20"]
		puts $log "PR3: $primerDict3" ; #TEST
		set k 0
		
		while {$k < [dict get $primerDict3 PRIMER_RIGHT_NUM_RETURNED]} {	
			set primerRR [list [dict get $primerDict3 PRIMER_RIGHT_${k}_SEQUENCE] \
				[dict get $primerDict3 PRIMER_RIGHT_${k}_TM] [dict get $primerDict3 PRIMER_RIGHT_${k}_GC_PERCENT] ]
			set prod_sizeR [dict get $primerDict3 PRIMER_PAIR_${k}_PRODUCT_SIZE]
			#run epcr on the 2 primers. There has to be 1 amplicon amplified
			set inv 0
			if {[catch "cg_validatesv_runEPCR {$primerRF} {$primerRR} $prod_sizeR $size $inv" out ]} {
				puts $log "out3: $out"
				incr k; continue
			}
			if {$out != "repeat"} {
				set PP 1
				break
			}
			incr k
		}
		if {[info exists PP]} { break } ; #when there is already a primerpair found, stop searching
		incr j
	}
	if {[info exists PrimerPair]} {
		set primerRF [lindex $PrimerPair 0]
		set beginRF [lindex $PrimerPair 1]
		set primerRR [lindex $PrimerPair 2]
		set beginRR [lindex $PrimerPair 3]
		set prod_sizeR [lindex $PrimerPair 4]	
		set label "inv_${chr}_[string range $patchstart 0 [expr [string length $patchstart] - 7]]"
		#puts $log "primRF_j : $primerRF" ; #TEST
		#making output file of all the 4 primers and there info

		set primLF "${label}_LF [lindex $primerLF 0] {} {} {} {pg} {1} {Hs} $chr {} {} {} $beginLF [lindex $primerLF 1]	{} "
		set primLR	"${label}_LR [lindex $primerLR 0] {} {} {} {pg} {1} {Hs} $chr {} {} {} $beginLR [lindex $primerLR 1] {} $prod_sizeL "
		set primRF "${label}_RF [lindex $primerRF 0] {} {} {} {pg} {1} {Hs} $chr {} {} {} $beginRF [lindex $primerRF 1]	{}"
		set primRR "${label}_RR [lindex $primerRR 0] {} {} {} {pg} {1} {Hs} $chr {} {} {} $beginRR [lindex $primerRR 1] {} $prod_sizeR "

		if {[llength $PrimerPair] > 5} {set repeatR 1 }
		if {[info exists repeatL] && [info exists repeatR]} {
			puts "WARNING: Both primerpairs are located in a repeat"
		} elseif {[info exists repeatL]} {
			puts "WARNING: Primerpair around left breakpoint is located in a repeat"
		} elseif {[info exists repeatR]} {
			puts "WARNING: Primerpair around right breakpoint is located in a repeat"
		}
		return [list $primLF $primLR $primRF $primRR]
	}
	
	return 1
}

proc cg_validatesv args {
	global log maxnum max_amplicon prim_min MAX_SIZE archive extraseq cachedir 
	global searchGenomeDB
	# (ugly) globals
	set maxnum 1000
	set max_amplicon 1000
	set prim_min 25 ; # the minimal distance between 2 masked sequences 'closest' to breakpoint
	set MAX_SIZE 4000
	set archive may2009
	set extraseq 1000
	set cachedir [pwd]/cache
	file mkdir $cachedir
	set log [open testlog.txt a] ; #TEST
	#
	#
	# set options
	# -----------
	set READSIZE 360
	set EVAL 5
	set MIN 100
		
#	foreach {key value} $args {
#		switch -- $key \
#			"-h" - "--help" {cg_validatesv_help ; exit 0} \
#			"-f" "set file $value" \
#			"-r" "set READSIZE $value" \
#			"-m" "set MIN $value" \
#			"-o" "set file_out $value" \
#			"-e" "set EVAL $value" ;
#	}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-r - --readsize {
				set READSIZE $value
			}
			-m - --min {
				set MIN $value
			}
			-e - --eval {
				set EVAL $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 4) || ([llength $args] > 5)} {
		puts stderr "Wrong number of arguments"
		errorformat cg_validatesv
		exit 1
	}
	foreach {file file_out dbdir archive MAX_SIZE} $args break
	#
	set searchGenomeDB [lindex [glob $dbdir/genome_*.ssa] 0]
	if {[catch {gzopen $file} fileid]} {
		error "Could not open input file $file: $fileid"
	}
	#
	# get headers
	# -----------
	set headers [tsv_open $fileid]
	#
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
	#
	# setting header in outputfile
	# ----------------------------
	if {[catch {open $file_out w} fileid_out]} {
		puts "Could not open output file - $fileid_out"
		exit 1
	}
	#new_output
	set line_out "label sequence modification scale purification project pair species chromosome cyto target contig pos temperature mg size_amplicon "
	#set line_out "chr patchStart patchEnd primer sequenceF sizePrimerF TmF GC-contentF sequenceR sizePrimerR TmR GC-contentR sizeAmplicon"
	puts $fileid_out [join $line_out \t]
	#
	# Getting the primerpairs for each inversion
	# ------------------------------------------
	set in [gets $fileid]
	set inversion_count 1
	set fasta 1
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
		
		if {$size <= 400} {
			puts "The inversion is small enough to sequence the whole thing at once."
			set primerPairs [cg_validatesv_getOnePair $dbdir $chr $patchstart $breakpointL $breakpointR $size $EVAL $MIN $fasta]
		} else {
			set primerPairs [cg_validatesv_getPrimerPairs $dbdir $chr $patchSize $patchstart $size $breakpointL $breakpointR $READSIZE $EVAL $MIN $fasta]
		}
		if {$primerPairs == 1} {
			puts "No primer pairs were found for this inversion"
			puts " "
			# ask first for rerun scripts with other parameters
			set okvar 0
			while {!$okvar} {
				if {![info exists againLarge] } {
					puts "Rerun this inversion with opportunity for larger amplicon (1200nt)? (Y/N)"
					set LargerAmplicon [gets stdin]	
					if {$LargerAmplicon eq {Y}} {
						set okvar 1
						set againLarge 1
						break
					} elseif {$LargerAmplicon eq {N}} {
						set okvar 1
					} else {
						puts "$LargerAmplicon is not a valid answer, use Y or N."
					}
				}
				if {![info exists againFasta]} {
					puts "Rerun this inversion but don't mask with FASTA? (Y/N)"
					set noFasta [gets stdin]
					if {$noFasta eq {Y}} {
						set okvar 1
						set againFasta 1
						break
					} elseif {$noFasta eq {N}} {
						set okvar 1				
					} else {
						puts "$noFasta is not a valid answer, use Y or N."
					} 
				}
				#
				puts "Rerun this inversion but tolerate more primer hits in genome (but still only 1 amplicon!)? (Y/N)"
				set moreHits [gets stdin]
				if {$moreHits eq {Y}} {
					set okvar 1
					set againHits 1
				} elseif {$moreHits eq {N}} {
					set okvar 1	
					puts "Ok, No primers are found for this inversion"				
				} else {
					puts "$moreHits is not a valid answer, use Y or N."
				}				
			}
		} else {
			puts "Primerpairs for ${chr}_${patchstart} are found!"
			puts " "
			foreach primer $primerPairs {
				#new output
				#set line_out "$chr $patchstart [lindex $line $PATCHEND] $primer"
				set line_out $primer
				puts $fileid_out [join $line_out \t]
			}
			unset -nocomplain againLarge
			unset -nocomplain againFasta
			unset -nocomplain againHits
		}
		#rerun this inversion
		if {[info exists againLarge]} {
			tempfile clean
			set max_amplicon 1200
			continue
		}
		if {[info exists againFasta]} {
			tempfile clean
			set fasta 0
			continue
		}
		if {[info exists againHits]} {
			tempfile clean
			incr maxnum 5000
			continue
		}
		
		tempfile clean 
		set in [gets $fileid]
		#reset original values for following inversion
		set max_amplicon 1000
		set maxnum 1000
		set fasta 1
		
		incr inversion_count
	}
	close $fileid
	puts " ----------------------------------------------------------------------------------"
	puts " \t \t \t Program is finished! "
	puts " Output is written to: $file_out "
	puts ""
	puts " Please cite for use of Primer3: "
	puts " Rozen S, Skaletsky H (2000) Primer3 on the WWW for general users and for biologist programmers. "
	puts " In: Krawetz S, Misener S (eds) Bioinformatics Methods and Protocols: Methods in Molecular Biology. "
	puts " Humana Press, Totowa, NJ, pp 365-386   "
	puts " ----------------------------------------------------------------------------------"
	puts ""
	close $log ; #TEST
	return 0
}
