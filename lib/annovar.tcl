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
#  adapted Peter De Rijk
#
########################################################################

if 0 {

	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	append env(PATH) :/home/peter/dev/completegenomics/bin
	package require Extral
	cd /complgen/projects/dlb1
	set file dlb_compar.tsv
	set path /complgen/refseq/hg18

}



proc annovar {file resultfile path {build hg18}} {

	########################################
	#
	# Making a new input file from the old
	#
	########################################
	catch {close $fileid} ; catch {close $fileid_out}
	puts "Making the right input file...."
	if {[catch {open $file r} fileid]} {
		puts "Could not open input file - $fileid"
		exit 1
	}
	set temp $file.annovar
	if {[catch {open $temp w} fileid_out]} {
		puts "Could not open new inputfile - $fileid_out"
		exit 1
	}
	set poss [open_region $fileid header]
	foreach {chrom begin end} $poss break
	set refpos [lsearch $header reference]
	if {$refpos == -1} {
		set refpos [lsearch $header ref]
	}
	if {$refpos == -1} {
		error "No reference position found"
	}
	set type [lsearch $header type]
	set in [gets $fileid]
	set allel_place [lsearch -all $header "alt"]
	set next 1000000; set num 0
	while {![eof $fileid]} {
		incr num; if {$num >= $next} {puts $num; incr next 1000000}
		set in [gets $fileid]
		set line [split $in "\t"]
		set use_allels [split [lindex $line $allel_place] ,]
		set refallele [lindex $line $refpos]
		foreach allele $use_allels {
			if {$refallele ne $allele} {
				if {$refallele eq ""} {set refallele -}
				if {$allele eq ""} {set allele -}
				set b [lindex $line $begin]
				set e [lindex $line $end]
				if {$e == $b} {incr e}
				puts $fileid_out [lindex $line $chrom]\t$b\t$e\t$refallele\t$allele\t[lindex $line $type]
			} else {
				#discard allele
			}
		}	
	}
	close $fileid
	close $fileid_out

	###########################################
	#
	# Run Annovar
	#
	###########################################
	set dbs {refGene knownGene ensGene}
	foreach db $dbs {

		# annovar lopen
		puts "Running annovar for $db database ....."
		puts "this may take a while......"
		catch {exec $path/annotate_variation.pl --buildver $build -geneanno -zerostart -dbtype $db $temp $path/humandb} errmsg
		if {[regexp Error $errmsg] || ![file exists $temp.variant_function]} {
			error $errmsg
		}
		#
		catch {close $fileid}; catch {close $fileid_old}; catch {close $fileid_out};
		set typeA 0
		set nameA 1
		set chromA 2
		set beginA 3
		set endA 4
		set nucA 6
		set commentA 7
		foreach which_file {variant_function exonic_variant_function} {
			if {[file exists $temp.annot.${db}.$which_file ]} continue
			# data retrieval from annovar file
			puts "retrieving data from annovar file $db.$which_file....."
			if {[catch {open $temp.$which_file r} fileid]} {
				puts "Could not open .[lindex $which_file $i] inputfile - $fileid"
				exit 1
			}
			if {[catch {open $file r} fileid_old]} {
				puts "Could not open inputfile - $fileid_old"
				exit 1
			}
			if {[catch {open $temp.annot w} fileid_out]} {
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
			if {$which_file == "variant_function" } {
				set header "${db}_type ${db}_name"
				puts $fileid_out [join $header \t]
				set l 0
				set k 0
				while {![eof $fileid_old]} {
					if {[expr {$l%1000000}]==0} {	
						puts "reading line $l"
					}
					set begin_old [lindex $line_old $begin]
					set end_old [lindex $line_old $end]
					set type_old [lindex $line_old $type]
					if {$type_old eq "ins"} {incr end_old}
					if {[lindex $line1 $chromA] == [lindex $line_old $chrom] && [lindex $line1 $beginA] == $begin_old && [lindex $line1 $endA] == $end_old && [lindex $line1 $commentA] == $type_old} {
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
			} elseif {$which_file == "exonic_variant_function" } {
				set header "${db}_exonic_type ${db}_exonic_effect"
				puts $fileid_out [join $header \t]
				set l 0
				set k 0
				while {![eof $fileid_old]} {
					if {[expr {$l%1000000}]==0} {	
						puts "reading line $l"
					}
					set begin_old [lindex $line_old $begin]
					set end_old [lindex $line_old $end]
					set type_old [lindex $line_old $type]
					if {$type_old eq "ins"} {incr end_old}
					if {[lindex $line1 $chromA] == [lindex $line_old $chrom] && [lindex $line1 $beginA] == $begin_old && [lindex $line1 $endA] == $end_old && [lindex $line1 $commentA] == $type_old} {
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
			if {[catch "exec cp $temp.annot $temp.annot.${db}.$which_file " errmsg]} {
				puts "Copy annotated info to file failed - $errmsg"
			}
		}

	}

#	######################################
#	#
#	# Running sift in annovar
#	#
#	######################################
#	# annovar lopen
#	puts "Running annovar for avsift database ....."
#	puts "this may take a while......"
#	catch {exec $path/annotate_variation.pl --buildver $build -filter -sift_threshold 0 -zerostart -dbtype avsift $temp $path/humandb} errmsg
#	if {![file exists $temp.hg18_avsift_dropped]} {
#		puts "error getting sift data: $errmsg"
#	}
#
#	######
#	#
#	# Creating header for sorting
#	#
#	######
#	if {[catch "open $temp.header w" fileout_header]} {
#		puts "Could not open header file - $fileout_header"
#		exit 1
#	}	
#	set headers "avsift score chrom begin end ref allel comment"
#	puts $fileout_header [join $headers \t]
#	close $fileout_header
#	if {[catch "exec cat $temp.header $temp.hg18_avsift_dropped > $temp.hg18_avsift_dropped_header" errmsg]} {
#		puts "something went wrong while adding a header to input file - $errmsg"
#		exit 1
#	}
#	if {[catch "exec cg select -s {chrom end} $temp.hg18_avsift_dropped_header $temp.hg18_avsift_dropped_sort" errmsg]} {
#		puts "Something went wrong while sorting - $errmsg"
#		exit 1
#	}
#
#	######
#	#
#	# extract info from annovar files
#	#
#	###### 
#	# data retrieval 
#	puts "retrieving data from annovar files....."
#	if {[catch "open $temp.hg18_avsift_dropped_sort r" fileid]} {
#		puts "Could not open .hg18_avsift_dropped inputfile - $fileid"
#		exit 1
#	}
#	if {[catch "open $file r" fileid_old]} {
#		puts "Could not open old inputfile - $fileid_old"
#		exit 1
#	}
#	if {[catch "open $temp.annot w" fileid_out]} {
#		puts "Could not open outputfile - $fileid_out"
#		exit 1
#	}
#	set header "avsift_score avsift_minimalScore"
#	puts $fileid_out [join $header \t]
#	set scoreA 1
#	set chromA 2
#	set beginA 3
#	set endA 4
#	set nucA 6
#	set commentA 7
#	#skipping headers
#	set in_old [gets $fileid_old]
#	set in1 [gets $fileid]
#	#reading real lines
#	set in1 [gets $fileid]
#	set line1 [split $in1 "\t"]
#	set in_old [gets $fileid_old]
#	set line_old [split $in_old "\t"]
#	set l 0
#	set k 0
#	while {![eof $fileid_old]} {
#		if {[expr {$l%1000000}]==0} {	
#			puts "reading line $l"
#		}
#		if {[lindex $line1 $chromA] == [lindex $line_old $chrom] && [lindex $line1 $beginA] == [lindex $line_old $begin] && [lindex $line1 $endA] == [lindex $line_old $end] && [lindex $line1 $commentA] == [lindex $line_old $type]} {
#			set newscore "[lindex $line1 $nucA]:[lindex $line1 $scoreA]"
#			set in [gets $fileid]
#			set line [split $in "\t"]
#			set min_score [lindex $line1 $scoreA]
#			incr k
#			while {[lindex $line1 $beginA] == [lindex $line $beginA] && [lindex $line1 $endA] == [lindex $line $endA] && [lindex $line1 $commentA] == [lindex $line $commentA] }	{
#				lappend newscore "[lindex $line $nucA]:[lindex $line $scoreA]"
#				if {$min_score > [lindex $line $scoreA]} {
#					set min_score [lindex $line $scoreA]
#				}
#				set in [gets $fileid]
#				set line [split $in "\t"]
#				incr k
#			}
#			set n_score [join $newscore ";"]
#			set newline "$n_score $min_score"
#			puts $fileid_out [join $newline \t]
#			set line1 $line
#			set in_old [gets $fileid_old]
#			set line_old [split $in_old "\t"]
#			incr l
#		} else {	
#			set newline {" " " "}
#			puts $fileid_out [join $newline \t]
#			set in_old [gets $fileid_old]
#			set line_old [split $in_old "\t"]
#			incr l
#		}
#	}
#	puts "Done reading line: $l"
#	close $fileid
#	close $fileid_old
#	close $fileid_out

	set dbresults [glob $temp.annot.refGene.variant_function $temp.annot.refGene.exonic_variant_function $temp.annot.knownGene.variant_function $temp.annot.knownGene.exonic_variant_function $temp.annot.ensGene.variant_function $temp.annot.ensGene.exonic_variant_function]
	exec paste {*}$dbresults $temp.annot > $resultfile.temp
	file rename $resultfile.temp $resultfile
}
