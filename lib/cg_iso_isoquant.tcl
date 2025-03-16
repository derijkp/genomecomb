proc convert_isoquant_ambigcount {ambig umicount} {
	if {$ambig > 1} {
		if {$umicount > 1} {
			return [expr {1.0/($ambig * $umicount)}]
		} else {
			return [expr {1.0/$ambig}]
		}
	} elseif {$umicount > 1} {
		return [expr {1.0/$umicount}]
	} else {
		return 1
	}
	return $result
}

proc convert_isoquant_add_ambig {varVar ambig} {
	upvar $varVar var
	if {![info exists var]} {
		if {$ambig} {
			set var [expr {1.0/$ambig}]
		} else {
			set var 1
		}
	} else {
		if {$ambig} {
			set var [expr {$var + 1.0/$ambig}]
		} else {
			set var [expr {$var + 1}]
		}
	}
}

proc convert_isoquant_reggenedb {reftranscripts samregions refseq regreftranscriptsVar reggenedbVar} {
	upvar $regreftranscriptsVar regreftranscripts
	upvar $reggenedbVar reggenedb
	set regreftranscripts [tempfile].tsv
	set reggenedb [tempfile].gtf
	set regfile [tempfile].tsv
	if {$samregions eq ""} {
		set regreftranscripts $regreftranscripts[gzext $reftranscripts]
		mklink $reftranscripts $regreftranscripts
	} else {
		distrreg_reg2tsv $regfile $samregions $refseq
		cg regselect $reftranscripts $regfile > $regreftranscripts
	}
	# if empty
	set genecol [lindex [list_common {gene_id geneid gene} [cg select -h $regreftranscripts]] 0]
	cg_tsv2gtf -genecol $genecol -addgene 1 $regreftranscripts $reggenedb
}

proc convert_isoquant_add {varVar {count 1}} {
	upvar $varVar var
	if {![info exists var]} {
		set var $count
	} else {
		set var [expr {$var + $count}]
	}
}

proc iso_isoquant_add_gcounts {gcountaVar dgeneisos gene gambigcount} {
	upvar $gcountaVar gcounta
	if {[llength $dgeneisos]} {
		if {![info exists gcounta($gene,w)]} {
			set gcounta($gene,w) $gambigcount
		} else {
			set gcounta($gene,w) [expr {$gcounta($gene,w) + $gambigcount}]
		}
	}
	if {![info exists gcounta($gene,i)]} {
		set gcounta($gene,i) $gambigcount
	} else {
		set gcounta($gene,i) [expr {$gcounta($gene,i) + $gambigcount}]
	}
}

proc iso_isoquant_add_tcounts {tcountaVar iso assignment_type ambiguity polya covered_pct strictpct ambigcount {count 1}} {
	upvar $tcountaVar tcounta
	convert_isoquant_add tcounta($iso,t) $ambigcount
	if {$ambiguity <= 1 && $assignment_type ne "inconsistent"} {
		# unique
		convert_isoquant_add tcounta($iso,u) $ambigcount
		if {$covered_pct >= $strictpct} {
			# strict
			convert_isoquant_add tcounta($iso,s) $ambigcount
		}
	}
	if {$polya eq "True"} {
		convert_isoquant_add tcounta($iso,a) $ambigcount
		if {$ambiguity == 1 && $assignment_type ne "inconsistent"} {
			# unique
			convert_isoquant_add tcounta($iso,au) $ambigcount
			if {$covered_pct >= $strictpct} {
				# strict
				convert_isoquant_add tcounta($iso,as) $ambigcount
			}
		}
	}
}

proc iso_isoquant_add_gcount {cgcountaVar iso {count 1}} {
	upvar $cgcountaVar cgcounta
	if {![info exists cgcounta($iso)]} {
		set cgcounta($iso) $count
	} else {
		set cgcounta($iso) [expr {$cgcounta($iso) + $count}]
	}
}

proc convert_isoquant_intron_getgene_init {tocheckVar genelist} {
	upvar $tocheckVar tocheck
	set tocheck(genelist) $genelist
	set tocheck(glen) [llength $genelist]
	set tocheck(gpos) 0
	set tocheck(+) {}
	set tocheck(-) {}
	set tocheck(.) {}
}

proc convert_isoquant_intron_getgene {tocheckVar chr begin end strand unstranded} {
	upvar $tocheckVar tocheck
	set genelist $tocheck(genelist)
	set gpos $tocheck(gpos)
	set glen $tocheck(glen)
	# add to tocheck until we are sure to be past current read alignment
	if {$strand eq "."} {
		set stopa(-) 0 ; set stopa(+) 0
		set rstrands {+ -}
	} else {
		set stopa(-) 1 ; set stopa(+) 1
		set stopa($strand) 0
		set rstrands [list $strand]
	}
	while {$gpos < $glen} {
		set gline [lindex $genelist $gpos]
		foreach {gchr gbegin gend gstrand ggene ggeneid} $gline break
		# if {[info exists parenta($ggeneid)]} {incr gpos ; continue}
		set chrcomp [loc_compare $gchr $chr]
		# if $gchr > $chr -> break
		if {$chrcomp > 0} break
		# keep adding if not same chromosome or strand
		# only stop if get past current alignment (begin) on same chr and strand
		if {$unstranded || $strand eq "."} {
			if {$chrcomp == 0 && $gbegin >= $begin} {
				set stopa($gstrand) 1
				if {$stopa(-) && $stopa(+)} break
			}
		} else {
			if {$chrcomp == 0 && $gstrand eq $strand && $gbegin >= $begin} break
		}
		lappend tocheck($gstrand) $gline
		incr gpos
	}
	set tocheck(gpos) $gpos
	# find a match in tocheck, remove from tocheck what we can no longer match (if we come across it)
	foreach rstrand $rstrands {
		set remove {}
		set pos [llength $tocheck($rstrand)]
		set match 0
		while {[incr pos -1] >= 0} {
			set cline [lindex $tocheck($rstrand) $pos]
			foreach {cchr cbegin cend cstrand cgene cgeneid} $cline break
			if {$chr eq $cchr && $begin >= $cbegin && $end <= $cend} {
				set match 1
				break
			}
			if {[loc_lt $cchr $chr] || $begin >= $cend} {lappend remove $cline}
		}
		if {$match} break
		if {[llength $remove]} {set tocheck($rstrand) [list_lremove $tocheck($rstrand) $remove]}
	}
	if {$match} {
		if {$cgeneid ne ""} {set cgene $cgeneid}
		return $cgene
	} else {
		return .
	}
}

proc convert_isoquant {isodir destdir sample refseq reggenedb regreftranscripts root singlecell} {
	set unstranded 0
	set strictpct 90
	set read_assignmentsfile [gzfile $isodir/*.read_assignments.tsv]
	set targetisoformcountsfile $destdir/isoform_counts-${root}.tsv
	set targetgenecountsfile $destdir/gene_counts-${root}.tsv
	# info from known isos
	# --------------------
	# transcriptidsa gets the transcripts that are actually in the results (from readassignment file)
	# sizea will contain the size of the known transcripts
	# transcriptidsa: which transcript ids are actually in known file
	# transcriptsa: transcript (iso_name) to known name (transcript_id)
	# geneid2genea: translate geneid 2 gene
	# genebasica: basic gene info (chr begin end strand)
	# transcript2genea: get gene for transcript (oriname)
	# geneconva: original gene name to final genename
	unset -nocomplain sizea
	unset -nocomplain transcriptidsa
	unset -nocomplain transcriptsa
	unset -nocomplain geneid2genea
	unset -nocomplain genebasica
	unset -nocomplain transcript2genea
	unset -nocomplain geneconva
	array set transcriptidsa [split [cg select -hc 1 -g isoform_id $read_assignmentsfile] \n\t]
	catch {close $f}
	if {$regreftranscripts ne ""} {
		set f [gzopen $regreftranscripts]
		set header [tsv_open $f]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 11 12 13}]
		foreach {isopos genepos geneidpos} [lrange $poss end-2 end] break
		while 1 {
			if {[gets $f line] == -1} break
			set split [split $line \t]
			foreach {chr begin end strand starts ends oriname gene geneid} [list_sub $split $poss] break
			if {![info exists genebasica($geneid)]} {
				set genebasica($geneid) [list $chr $begin $end $strand]
			} else {
				set prevbegin [lindex $genebasica($geneid) 2]
				if {[lindex $genebasica($geneid) 0] ne "$chr"} {
					puts stderr "Chromosome does not match for transcripts of same gene $geneid: $line"
				} elseif {[lindex $genebasica($geneid) 3] ne "$strand"} {
					puts stderr "Strand does not match for transcripts of same gene $geneid: $line"
				} elseif {$begin > $prevbegin} {
					puts stderr "transcript not overlapping previous for same gene $geneid: $line"
					if {[expr {$begin - $prevbegin}] < 100000} {
						lset genebasica($geneid) 2 $end
					}
				} elseif {$end > $prevbegin} {
					lset genebasica($geneid) 2 $end
				}
			}
			set transcript2genea($oriname) $geneid
			if {$gene ne $geneid} {set geneid2genea($geneid) $gene}
			set transcript [iso_name $chr $strand $starts $ends size]
			set sizea($oriname) $size
			set transcriptsa($transcript) $oriname
		}
		close $f
	}
	# gene info from extended
	set file [gzfile $isodir/*.extended_annotation.gtf]
	if {[file exists $file]} {
		set f [open $file]
		while {1} {
			if {[gets $f line] == -1} break
			if {[string index $line 0] eq "#"} continue
			set line [split $line \t]
			foreach {chr src type begin end temp strand temp info} $line break
			if {$type ne "gene"} continue
			if {![regexp {gene_id "([^"]+)"} $info temp gene]} continue
			if {[regexp ^novel_gene $gene]} {
				set geneconva($gene) [gene_name $chr $strand $begin $end]
				set gene $geneconva($gene)
			}
			if {![info exists genebasica($gene)]} {
				set genebasica($gene) [list $chr $begin $end $strand]
			}
		}
		close $f
	}
	#
	#
	# info from model isos
	# --------------------
	# geneconva: to convert ori genename (for novel) to new
	catch {close $f}
	set mfile [gzfile $isodir/*.transcript_models.gtf]
	if {[file exists $mfile]} {
		set f [open $mfile]
		while {[gets $f line] != -1} {
			set line [split $line \t]
			foreach {chr src type begin end temp strand temp info} $line break
			if {$type ne "gene"} continue
			if {![regexp {gene_id "([^"]+)"} $info temp gene]} continue
			if {[regexp ^novel_gene $gene]} {
				set geneconva($gene) [gene_name $chr $strand $begin $end]
				set gene $geneconva($gene)
			}
			if {![info exists genebasica($gene)]} {
				set genebasica($gene) [list $chr $begin $end $strand]
			}
		}
		close $f
		#
		# get convert data from models file
		# converta: models that fully match known transcripts (oriname -> dbname)
		# outputa: all transcripts that must be output separately from known (oriname -> newname)
		# modelregiona: basic regions for transcripts -> needed for correctly assigning reads to novel transcripts when there are multiple alignment locations
		unset -nocomplain converta
		unset -nocomplain outputa
		unset -nocomplain modelregiona
		cg_gtf2tsv $mfile $destdir/transcripts_models-$sample.tsv.temp
		# check for and fix transcript out of gene area error
		# these were double assignments (one wrong), so wrong ones can be filtered out
		set f [gzopen $destdir/transcripts_models-$sample.tsv.temp]
		set header [tsv_open $f]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 11 12 13 7 8}]
		set o [wgzopen $destdir/transcripts_models-$sample.tsv.temp2]
		puts $o [join $header \t]
		while 1 {
			if {[gets $f line] == -1} break
			foreach {c b e s iso g gid es ee} [list_sub [split $line \t] $poss] break
			if {[info exists genebasica($g)]} {
				foreach {gc gb ge gs} $genebasica($g) break
				if {$s ne $gs || $c ne $gc || $e < $gb || $b >= $ge} {
					puts "skipping $iso ($g) because not in gene region $genebasica($g): $line"
					continue
				}
			}
			if {[regexp ^transcript $iso]} {
				set modelregiona($iso) [list $chr [lindex [split $es ,] 0] [lindex [split $ee ,] end] $strand]
			}
			puts $o $line
		}
		gzclose $o
		close $f
		file rename -force $destdir/transcripts_models-$sample.tsv.temp2 $destdir/transcripts_models-$sample.tsv
		#
		catch {close $f}
		set f [open $destdir/transcripts_models-$sample.tsv]
		set header [tsv_open $f comments]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 6 7 8 11 12 13}]
		foreach {isopos genepos geneidpos} [lrange $poss end-2 end] break
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {chr strand starts ends oriname gene geneid} [list_sub $line $poss] break
			set name [iso_name $chr $strand $starts $ends size]
			set sizea($oriname) $size
			if {[info exists geneconva($geneid)]} {set gene $geneconva($geneid)}
			if {[info exists geneconva($gene)]} {set gene $geneconva($gene)}
			if {[info exists transcriptsa($name)]} {
				set converta($oriname) $transcriptsa($name)
				if {![info exists transcriptidsa($transcriptsa($name))]} {
					set outputa($oriname) $transcriptsa($name)
				}
				set transcript2genea($transcriptsa($name)) $geneid
			} else {
				# set outputa($oriname) [list $name {*}[lrange $line 3 end]]
				set outputa($oriname) $name
				set transcript2genea($oriname) $geneid
			}
		}
		close $f
	}
	set genelist {}
	foreach gene [array names genebasica] {
		lappend genelist [list {*}$genebasica($gene) $gene]
	}
	set genelist [bsort $genelist]
	# join $genelist \n
	#
	# umi counts
	# ----------
	unset -nocomplain umicounta
	unset -nocomplain umicount_readsa
	catch {close $f}
	set f [gzopen $read_assignmentsfile]
	set prev {}
	while {[gets $f line] != -1} {
		if {[string index $line 0] ne {#}} break
		set prev $line
	}
	set header [split [string range $prev 1 end] \t]
	set readpos [lsearch $header read_id]
	while 1 {
		set line [split $line \t]
		set read [lindex $line $readpos]
		if {[regexp {^([A-Z]+)_([A-Z]+)#(.*)$} $read temp cellbarcode umi r]} {
			if {![info exists umicount_readsa($r)]} {
				incr umicounta($cellbarcode,$umi)
				set umicount_readsa($r) 1
			}
		}
		if {[gets $f line] == -1} break
	}
	close $f
	unset -nocomplain umicount_readsa
	#
	unset -nocomplain transcriptidsa
	#
	# read_assignments (models)
	# -------------------------
	# load model transcripts read info in read2isoa
	unset -nocomplain read2isoa
	set mfile [gzfile $isodir/*.transcript_model_reads.tsv]
	if {[file exists $mfile]} {
		catch {close $f} ; set f [gzopen $mfile]
		gets $f
		while {[gets $f line] != -1} {
			foreach {read transcript} $line break
			# only keep novel transcripts that are in the output models (so filter out known, *, and novel not in the output gtf)
			# replace these with *
			if {[info exists outputa($transcript)]} {
				list_addnew read2isoa($read) $transcript
			}
		}
		gzclose $f
	}
	#
	# check which reads have ambiguous mappings (from original read_assignmentsfile)
	unset -nocomplain ambiga
	foreach {read exons count} [cg select -hc 1 -g {read_id * exons *} $read_assignmentsfile] {
		dict set ambiga($read) $exons $count
	}
	foreach read [array names ambiga] {
		set ra [dict values $ambiga($read)]
		if {[llength $ra] <= 1 && [lindex $ra 0] <= 1} {
			unset ambiga($read)
		}
	}
	# foreach r [array names ambiga] {if {[llength $ambiga($r)] > 2} {puts [list set ambiga($r) $ambiga($r)]}}

	#
	# make read_assignments file
	# --------------------------
	# set target $destdir/read_assignments-${root}.tsv
	set temptarget $destdir/read_assignments-${root}.tsv.temp
	catch {close $f} ; catch {close $o}
	set target $destdir/read_assignments-${root}.tsv
	set f [gzopen $read_assignmentsfile]
	set prev {}
	while {[gets $f line] != -1} {
		if {[string index $line 0] ne {#}} break
		set prev $line
	}
	set header [split [string range $prev 1 end] \t]
	set newheader {read_id chromosome begin end strand exonStarts exonEnds aligned_size}
	lappend newheader {*}[list_remove $header read_id chromosome chr strand exons]
	set poss [list_cor $header $newheader]
	lset poss 1 [lsearch $header chr]
	set exonspos [lsearch $header exons]
	set chrpos [lsearch $newheader chromosome]
	set readpos [lsearch $newheader read_id]
	set isopos [lsearch $newheader isoform_id]
	set geneidpos [lsearch $newheader gene_id]
	set infopos [lsearch $newheader additional_info]
	set beginpos [lsearch $newheader begin]
	set endpos [lsearch $newheader end]
	set strandpos [lsearch $newheader strand]
	set exonStartspos [lsearch $newheader exonStarts]
	set exonEndspos [lsearch $newheader exonEnds]
	set assignment_typepos [lsearch $newheader assignment_type]
	set assignment_eventspos [lsearch $newheader assignment_events]
	set aligned_sizepos [lsearch $newheader aligned_size]
#	set nposs [list_cor $newheader {
#		chromosome begin end strand exonStarts exonEnds 
#		read_id isoform_id gene_id additional_info assignment_type aligned_size
#	}]
	lappend newheader ambiguity inconsistency covered_pct polya classification closest_known
	if {$singlecell} {
		lappend newheader cellbarcode umi umicount
	}
	#
	# set o [wgzopen $temptarget w {} 4]
	set o [open $temptarget w]
	puts $o [join $newheader \t]
	unset -nocomplain tcounta
	unset -nocomplain gcounta
	unset -nocomplain donea
	unset -nocomplain multirega
	for {set p 0} {$p <= 100} {incr p} {set pcta($p) 0}
	convert_isoquant_intron_getgene_init tocheck $genelist

	while 1 {
		set line [split $line \t]
		set exons [lindex $line $exonspos]
		set line [list_sub $line $poss]
		set read [lindex $line $readpos]
#if {$read eq "CATCGGGTCCCTCTCC_GTCAAAGCTGTC#cbff2f32-40d7-46ea-9128-7e6838ac2685"} {
#	error stop
#}
		if {$singlecell} {
			if {![regexp {([A-Z]+)_([A-Z]+)#} $read temp cellbarcode umi]} {
				set cellbarcode {} ; set umi {}
			}
		}
		set exonStarts {}
		set exonEnds {}
		set size 0
		foreach exon [split $exons ,] {
			foreach {s e} [split $exon -] break
			lappend exonStarts $s
			lappend exonEnds $e
			set size [expr {$size + ($e - $s)}]
		}
		set chr [lindex $line $chrpos]
		set strand [lindex $line $strandpos]
		set begin [lindex $exonStarts 0]
		set end [lindex $exonEnds end]
		set geneid [lindex $line $geneidpos]
		if {[info exists ambiga($read)]} {
			set temp [dict values $ambiga($read)]
			set ambig [lindex $temp 0]
			foreach v [lrange $temp 1 end] {
				incr ambig $v
			}
		} else {
			set ambig 0
		}
		set assignment_type [lindex $line $assignment_typepos]
		lset line $beginpos $begin
		lset line $endpos $end
		lset line $exonStartspos [join $exonStarts ,]
		lset line $exonEndspos [join $exonEnds ,]
		lset line $aligned_sizepos $size
		set closest_known [lindex $line $isopos]
		set info [lindex $line $infopos]
		set newinfo {}
		set polya {}
		set classification {}
		foreach temp [split $info \;] {
			set temp [string trim $temp]
			set pos [string first = $temp]
			if {$pos != -1} {
				set key [string range $temp 0 [expr {$pos-1}]]
				set value [string range $temp [expr {$pos+1}] end]
				if {$key eq "PolyA"} {
					set polya $value
				} elseif {$key eq "Classification"} {
					set classification $value
				} else {
					lappend newinfo $temp
				}
			} elseif {$temp ne ""} {
				lappend newinfo $temp
			}
			lset line $infopos $newinfo
		}
# if {$read eq "dd5ac6a9-2c67-4d86-a628-a0c74e8d6619"} {error selread}
		set assignment_events [lindex $line $assignment_eventspos]
		if {$assignment_type in "unique"} {
			set inconsistency 0
		} elseif {[regexp {major_exon_elongation|extra|alternative_structure|alt_donor_site|alt_acceptor_site|intron_alternation|migration|exclusive|skip|merge|gain|detach|terminal_exon_shift} $assignment_events]} {
			# significant inconsistencies
			set inconsistency 3
		} elseif {[regexp {intron_retention|unspliced_intron_retention|incomplete_intron_retention} $assignment_events]} {
			# intron retentions
			set inconsistency 2
		} elseif {[regexp {fake_terminal|alternative_polya|intron|exon} $assignment_events]} {
			# alignment artifacts, alternative transcription start / end
			set inconsistency 1
		} else {
			set inconsistency 0
		}
		if {[info exists read2isoa($read)]} {
			# if alternatives have no novel transcripts (because they e.g were not included in the final output),
			# remove and handle normally
			if {[lsearch -regexp $read2isoa($read) ^transcript] == -1} {
				unset read2isoa($read)
			}
		}
		if {[info exists read2isoa($read)]} {
			if {$assignment_type in "inconsistent unique_minor_difference intergenic" 
				|| $closest_known eq "."
			} {
				set knownmatch 0
			} else {
				if {$inconsistency > 1} {
					set knownmatch 0
				} else {
					set knownmatch 1
				}
			}
			if {[info exists ambiga($read)]} {
				if {[info exists multirega($read)]} {
					# multiple regions in readassignments, already encountered
				} else {
					set regions [dict keys $ambiga($read)]
					if {[llength $regions] > 1} {
						# multiple regions in readassignments, first encounter: try to find out which novels belong to which region
						# add isos to regions they overlap
						# this can still make some errors if we have a multilocal read where >1 are in the same region
						# overlapped by one iso
						# can only solve this fully by comparing full read exon structure iso and read
						# not doing this yet ..
						foreach iso $read2isoa($read) {
							foreach {mchr mstart mend mstrand} $modelregiona($iso) break
							set rpos 0
							foreach region $regions {
								set temp [split $region -]
								set rs [lindex $temp 0]
								set re [lindex $temp end]
								if {$rs >= $mstart && $re <= $mend} {
									dict lappend multirega($read) $region $iso
								}
								incr rpos
							}
						}
					}
				}
			}
			if {[info exists donea($exons,$read)]} {
				if {!$knownmatch} {
					# skip if already done: all (if more than one) knowns will be replaced by the novel ones
					if {[gets $f line] == -1} break
					continue
				} else {
					# new ones have been added already, proceed normally, adding this known one
					set isos [list $closest_known]
					set inconsistencylist [list $inconsistency]
				}
			} elseif {!$knownmatch} {
				# the readassignment is not a known, or is inconsistent with known -> replace with novel
				if {[info exists multirega($read)]} {
					if {[dict exists $multirega($read) $exons]} {
						set modelisos [dict get $multirega($read) $exons]
					} else {
						set modelisos {}
					}
				} else {
					set modelisos $read2isoa($read)
				}
				if {[llength $modelisos]} {
					# we have hits in models: use those, remove known hits (by setting donea)
					set isos [list_remdup $modelisos]
					set inconsistencylist [list_fill [llength $isos] 0]
					set ambig [llength $isos]
					if {$ambig > 1} {
						set assignment_type ambiguous_$assignment_type
					} else {
						set assignment_type unique_$assignment_type
					}
					lset line $assignment_typepos $assignment_type
					if {$ambig == 1} {set ambig 0}
					dict set ambiga($read) $exons $ambig
					# write only model hits
					set donea($exons,$read) 1
				} else {
					# not actually a replacement
					set isos [list $closest_known]
					set inconsistencylist [list $inconsistency]
				}
			} else {
				# the readassignment is known -> add novel ones to readassignment
				if {[info exists multirega($read)]} {
					if {[dict exists $multirega($read) $exons]} {
						set modelisos [dict get $multirega($read) $exons]
					} else {
						set modelisos {}
					}
				} else {
					set modelisos $read2isoa($read)
				}
				if {$ambig == 0} {set ambig 1}
				incr ambig [llength $modelisos]
				dict set ambiga($read) $exons $ambig
				lset line $assignment_typepos ambiguous
				set isos [list $closest_known {*}$modelisos]
				set inconsistencylist [list $inconsistency]
				lappend inconsistencylist {*}[list_fill [llength $modelisos] 0]
				set donea($exons,$read) 1
			}
		} else {
			set isos [list $closest_known]
			set inconsistencylist [list $inconsistency]
		}
		set ambigcount [convert_isoquant_ambigcount $ambig 1]
		if {$singlecell} {
			set umicount [get umicounta($cellbarcode,$umi) 1]
			if {$umicount > 1} {
				set count [expr {1.0/$umicount}]
			} else {
				set count 1
			}
		} else {
			set count 1
		}
		foreach iso $isos inconsistency $inconsistencylist {
			if {[info exists outputa($iso)]} {
				lset line $isopos $outputa($iso)
			} else {
				lset line $isopos $iso
			}
			if {$iso ne "."} {
				if {[info exists transcript2genea($iso)]} {
					set geneid $transcript2genea($iso)
					if {[info exists geneconva($geneid)]} {set geneid $geneconva($geneid)}
					lset line $geneidpos $geneid
				} else {
					set geneid [lindex $line $geneidpos]
					set transcript2genea($iso) $geneid
				}
			} else {
				set geneid [convert_isoquant_intron_getgene tocheck $chr $begin $end $strand $unstranded]
				lset line $geneidpos $geneid
			}
			if {[info exists sizea($iso)]} {
				set pct [expr {100*$size/$sizea($iso)}]
				if {$pct > 100} {set pct 100}
				incr pcta($pct)
			} else {
				set pct 0
			}
			if {$singlecell} {
				puts $o [join $line \t]\t$ambig\t$inconsistency\t$pct\t$polya\t$classification\t$closest_known\t$cellbarcode\t$umi\t$umicount
			} else {
				puts $o [join $line \t]\t$ambig\t$inconsistency\t$pct\t$polya\t$classification\t$closest_known
			}
			# counts
			if {$inconsistency < 2} {
				convert_isoquant_add tcounta($iso,t) $ambigcount
				if {!$ambig && $assignment_type ne "inconsistent"} {
					# unique
					convert_isoquant_add tcounta($iso,u) $count
					if {$pct >= $strictpct} {
						# strict
						convert_isoquant_add tcounta($iso,s) $count
					}
				}
				if {$polya eq "True"} {
					convert_isoquant_add tcounta($iso,a) $ambigcount
					if {!$ambig && $assignment_type ne "inconsistent"} {
						# unique
						convert_isoquant_add tcounta($iso,au) $count
						if {$pct >= $strictpct} {
							# strict
							convert_isoquant_add tcounta($iso,as) $count
						}
					}
				}
			}
			# don't try to add gambiguity here (per region)
			# we'll (have to) do it in the merge stage (to account for matches in other regions)
			convert_isoquant_add gcounta($geneid) $ambigcount
		}
		if {[gets $f line] == -1} break
	}

	gzclose $o
	gzclose $f
	cg select -overwrite 1 -s - $temptarget ${temptarget}2.zst
	file rename -force ${temptarget}2.zst $destdir/read_assignments-${root}.tsv.zst
	file delete $temptarget

	#
	# make isoform_counts file
	# ------------------------
	# read isoquant counts
	set f [gzopen [gzfile $isodir/*.transcript_counts.tsv]]
	gets $f
	foreach {iso value} [read $f] {
		set tcounta($iso,iq) $value
	}
	gzclose $f
	set mfile [gzfile $isodir/*.transcript_model_counts.tsv]
	if {[file exists $mfile]} {
		set f [gzopen $mfile]
		gets $f
		foreach {iso value} [read $f] {
			set tcounta($iso,iq) $value
		}
		gzclose $f
	}
	# open isoformcountsfile, write header and reference transcripts; novel transcripts will be added in next part
	set fields [list \
		counts_iqall-$root iq counts_weighed-$root t counts_unique-$root u counts_strict-$root s \
		counts_aweighed-$root a counts_aunique-$root au counts_astrict-$root as\
	]
	set o [iso_write_isoform_counts \
		$targetisoformcountsfile.temp \
		$regreftranscripts \
		tcounta \
		newheader \
		$fields \
	]
	#
	# models
	if {[file exists $destdir/transcripts_models-$sample.tsv]} {
		set f [open $destdir/transcripts_models-$sample.tsv]
		set header [tsv_open $f comments]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 11 12 13}]
		set remove [list_sub $header $poss]
		foreach {isopos genepos geneidpos} [lrange $poss 6 end] break
		lappend poss {*}[list_cor $header [lrange $newheader [llength $poss] end]]
	
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			set oriname [lindex $line $isopos]
			set geneid [lindex $line $geneidpos]
			if {![info exists outputa($oriname)]} continue
			set iq [get tcounta($oriname,iq) 0]
			set t [get tcounta($oriname,t) 0]
			if {$iq < 1 && $t < 1} continue
			if {[info exists converta($oriname)]} {
				set iso $converta($oriname)
				lset line $isopos $iso
				set category known
			} else {
				set gene [lindex $line $genepos]
				if {[regexp {novel} $gene]} {
					set category novel_gene
				} elseif {[regexp {\.nnic$} $oriname]} {
					set category novel_not_in_catalog
				} elseif {[regexp {\.nic$} $oriname]} {
					set category novel_in_catalog
				} else {
					set category novel
				}
			}
			lset line $isopos $outputa($oriname)
			if {[info exists geneconva($gene)]} {
				lset line $genepos $geneconva($gene)
				lset line $geneidpos $geneconva($gene)
			} else {
				lset line $genepos [get geneid2genea($geneid) $gene]
			}
			set line [list_sub $line $poss]
			# set genebasica($geneid) [list_sub $line {0 1 2 3}]
			if {![info exists genebasica($geneid)]} {
				set genebasica($geneid) [list_sub $line {0 1 2 3}]
			}
			puts $o [join $line \t]\t$category\t$sizea($oriname)\t$iq\t[format %.2f $t]\t[get tcounta($oriname,u) 0]\t[get tcounta($oriname,s) 0]\t[format %.2f [get tcounta($oriname,a) 0]]\t[get tcounta($oriname,au) 0]\t[get tcounta($oriname,as) 0]
		}
		close $f
	}
	gzclose $o
	cg select -overwrite 1 -s - $targetisoformcountsfile.temp $targetisoformcountsfile.temp2
	file delete $targetisoformcountsfile.temp
	file rename -force $targetisoformcountsfile.temp2 $targetisoformcountsfile
	#
	# make gene_counts file
	# ---------------------
	# file delete $targetisoformcountsfile
	unset -nocomplain gcounta(.)
	unset -nocomplain gcounta()
	set o [open $targetgenecountsfile.temp w]
	puts $o [join [list chromosome begin end strand gene geneid counts-${root}] \t]
	foreach origeneid [array names gcounta] {
		if {[info exists geneconva($origeneid)]} {
			set geneid $geneconva($origeneid)
			set gene $geneconva($origeneid)
		} else {
			set geneid $origeneid
			set gene [get geneid2genea($geneid) $origeneid]
		}
		puts $o [join $genebasica($geneid) \t]\t$gene\t$geneid\t$gcounta($origeneid)
	}
	close $o
	cg select -s - $targetgenecountsfile.temp $targetgenecountsfile.temp2
	file rename -force $targetgenecountsfile.temp2 $targetgenecountsfile
	file delete $targetgenecountsfile.temp
}

proc iso_isoquant_mergeresults_old_monmerge {isofiles genefiles readfiles strictpct sample {analysisname isoquant}} {
	set o [open isoform_counts-${root}.tsv.temp w]
	foreach refisofile $isofiles {
		set f [open $refisofile]
		set header [tsv_open $f comments]
		if {[gets $f line] == -1} {
			close $f
		} else {
			close $f
			break
		}
	}
	puts $o $comments[join $header \t]
	unset -nocomplain donea
	foreach isofile $isofiles {
		set f [open $isofile]
		set cheader [tsv_open $f comments]
		if {$cheader ne $header} {
			# if file is empty, diff header does not matter -> skip
			if {[gets $f line] == -1} continue
			error "header of $isofile differs from header of $refisofile"
		}
		set pos [lsearch $cheader name]
		while {[gets $f line] != -1} {
			set split [split $line \t]
			set name [lindex $split $pos]
			if {[info exists donea($name)]} {
				if {[regexp {^(transcript_[0-9]+)(.*)$} $name temp pre post]} {
					set num 2
					while 1 {
						set test $pre-$num$post
						if {![info exists donea($test)]} break
						incr num
					}
					set name $test
				} else {
					set num 2
					while 1 {
						set test $name-$num
						if {![info exists donea($test)]} break
						incr num
					}
					set name $test
				}
				lset split $pos $name
				puts $o [join $split \t]
			} else {
				puts $o $line
			}
			set donea($name) 1
		}
		close $f
	}
	close $o
	cg select -s - isoform_counts-${root}.tsv.temp isoform_counts-${root}.tsv.temp2
	file rename -force isoform_counts-${root}.tsv.temp2 isoform_counts-${root}.tsv
	file delete isoform_counts-${root}.tsv.temp
	#
	cg cat {*}$genefiles > gene_counts-${root}.tsv.temp
	file rename -force gene_counts-${root}.tsv.temp gene_counts-${root}.tsv
	#
	cg cat {*}$readfiles | cg zst > read_assignments-${root}.tsv.temp.zst
	file rename -force read_assignments-${root}.tsv.temp.zst read_assignments-${root}.tsv.zst
	# 
	cg select -overwrite 1 -g all -gc {sum(count*-*)} isoform_counts-${root}.tsv | cg select -rf all > totalcounts-${root}.tsv.temp
	file rename -force totalcounts-${root}.tsv.temp totalcounts-${root}.tsv
}

proc iso_isoquant_mergeresults {isofiles genefiles readfiles strictpct sample root {analysisname isoquant} {addumis 0}} {
	set tempreads [tempfile].tsv.zst
	set tempreads2 [tempfile].tsv.zst
	# set tempreads temp_sort_read_assignments-${root}.tsv.zst
	# set tempreads2 temp_addsort_read_assignments-${root}.tsv.zst
	# need to concat all and sort on read_id to find gambiguity
	# (alignments to different genes can be in different regions)
	putslog "Merging readfiles"
	cg cat -m 1 {*}$readfiles | cg select -s read_id | cg zst -c 1 > $tempreads

	catch {close $f} ; catch {close $o}
	putslog "gathering (bulk) counts"
	# gather counts in memory (tcounta,gcounta)
	# also write readassignment with added gambiguity in $tempreads2, to be sorted later
	unset -nocomplain tcounta
	unset -nocomplain gcounta
	unset -nocomplain utcounta
	unset -nocomplain ugcounta
	catch {close $f} ; catch {close $o}
	set f [gzopen $tempreads]
	set header [tsv_open $f comment]
	set poss [list_cor $header {assignment_type isoform_id covered_pct polya inconsistency}]
	set chrpos [lsearch $header chromosome]
	set readpos [lsearch $header read_id]
	set inconsistencypos [lsearch $header inconsistency]
	set sizepos [lsearch $header aligned_size]
	set transcriptpos [lsearch $header isoform_id]
	set isopos [lsearch $header isoform_id]
	set genepos [lsearch $header gene_id]
	set ambpos [lsearch $header ambiguity]
	set umicountpos [lsearch $header umicount]
	set pctpos [lsearch $header covered_pct]
	set o [wgzopen $tempreads2]
	puts $o $comment[join $header \t]\tgambiguity
	set prevgrouper {}
	set prevumi {}
	unset -nocomplain todoa
	set todo {}
	set read 0
	set nr 0

	#catch {close $odbg}
	# set odbg [open ~/tmp/genedebug.tsv w]
	# puts $odbg [join [list read_id gene gambigcount wcount icount dgenes dgeneisos genes] \t]
	set umi {}
	while 1 {
		incr nr
		if {![expr $nr%1000000]} {puts $nr}
		if {$read == -1} break
		set read [gets $f line]
		set line [split $line \t]
		set read_id [lindex $line $readpos]
		if {$read_id eq ""} continue
		if {$addumis} {
			if {[regexp {^([A-ZA-Z_]+)#(.*)$} $read_id temp umi r]} {
				set grouper $umi
			} else {
				set grouper $read_id
			}
		} else {
			set grouper $read_id
		}
		if {$grouper ne $prevgrouper} {
			set reads [array names todoa]
			set umicount [llength $reads]
			foreach read $reads {
				set todo $todoa($read)
				set remtodo {}
				set ambiguity [llength $todo]
				if {$ambiguity == 1} {
					set ambigcount 1
					set gambigcount 1
					set gambiguity 1
					set genes [list_subindex $todo $genepos]
					set isos [list_subindex $todo $isopos]
					set inconsistency [lindex $todo 0 $inconsistencypos]
					set ca([lindex $genes 0]) {1 1}
				} elseif {$ambiguity > 1} {
					set genes [list_subindex $todo $genepos]
					set isos [list_subindex $todo $isopos]
					set dgenes [list_remdup $genes]
					if {[llength $dgenes] > 1} {
						# filter out shorter alignments
						set sizes [list_subindex $todo $sizepos]
						set dsizes [list_remdup $sizes]
						if {[llength $dsizes] > 1} {
							set cutoff [expr {[lmath_max $dsizes]-100}]
							set temp {}
							foreach tline $todo {
								set size [lindex $tline $sizepos]
								if {$size >= $cutoff} {
									lappend temp $tline
								} else {
									lappend remtodo $tline
								}
							}
							set todo $temp
							set genes [list_subindex $todo $genepos]
							set isos [list_subindex $todo $isopos]
						}
						if {[llength $todo] > 1} {
							# filter out none-gene if there are gene
							set nogene [list_find $genes .]
							if {[llength $nogene] > 0 && [llength $nogene] < [llength $todo]} {
								lappend remtodo {*}[list_sub $todo $nogene]
								set todo [list_sub $todo -exclude $nogene]
								set genes [list_subindex $todo $genepos]
								set isos [list_subindex $todo $isopos]
							}
						}
					}
					if {[llength $todo] > 1} {
						# filter out non-transcript if there are transcript
						set transcripts [list_subindex $todo $transcriptpos]
						set notranscript [list_find $transcripts .]
						if {[llength $notranscript] > 0 && [llength $notranscript] < [llength $todo]} {
							lappend remtodo {*}[list_sub $todo $notranscript]
							set todo [list_sub $todo -exclude $notranscript]
							set genes [list_subindex $todo $genepos]
							set isos [list_subindex $todo $isopos]
						}
					}
					if {[llength $todo] > 1} {
						# filter out inconsistent (if there are more consistent)
						set incs [list_subindex $todo $inconsistencypos]
						set incs [lsort -integer [list_remdup $incs]]
						set inconsistency [lindex $incs 0]
						if {[llength $incs] > 1 && [llength [list_remove $incs 0 1]]} {
							if {$inconsistency == 0} {set set {0 1}} else {set set [list $inconsistency]}
							set temp {}
							foreach l $todo {
								set linc [lindex $l $inconsistencypos]
								if {$linc in $set} {
									lappend temp $l
								}
							}
							set todo $temp
							set ambiguity [llength $todo]
							set genes [list_subindex $todo $genepos]
							set isos [list_subindex $todo $isopos]
						}
					}
					set ambiguity [llength $todo]
					unset -nocomplain ca
					set dgenes [list_remdup $genes]
					foreach gene $dgenes {
						set ambiguity [llength [list_find $genes $gene]]
						set ambigcount [expr {1.0/$ambiguity}]
						set ca($gene) [list $ambiguity $ambigcount]
					}
					set gambiguity [llength $dgenes]
					set gambigcount [expr {1.0/$gambiguity}]
				}
				if {[llength $todo]} {
					set dgenes [list_remdup $genes]
					foreach gene $dgenes {
						set dgeneisos [list_remdup [list_remove [list_sub $isos [list_find $genes $gene]] .]]
						iso_isoquant_add_gcounts gcounta $dgeneisos $gene $gambigcount
						if {$addumis} {
							iso_isoquant_add_gcounts ugcounta $dgeneisos $gene [expr {$gambigcount/$umicount}]
						}
						# puts $odbg [join [list $read_id $gene $gambigcount [get gcounta($gene,w) 0] [get gcounta($gene,i) 0] $dgenes $dgeneisos $genes] \t]
						# flush $odbg
					}
					foreach l $todo {
						# add gamibiguity and write to (temp) read_assignments
						set gene [lindex $l $genepos]
						foreach {ambiguity ambigcount} $ca($gene) break
						lset l $ambpos $ambiguity
						if {$addumis} {
							lset l $umicountpos $umicount
						}
						lappend l $gambiguity
						puts $o [join $l \t]
						# count
						foreach {assignment_type iso covered_pct polya inconsistency} [list_sub $l $poss] break
						if {$inconsistency < 2 && $iso ne "."} {
							iso_isoquant_add_tcounts tcounta $iso $assignment_type $ambiguity $polya $covered_pct $strictpct $ambigcount
							if {$addumis} {
								iso_isoquant_add_tcounts utcounta $iso $assignment_type $ambiguity $polya $covered_pct $strictpct [expr {$ambigcount/$umicount}]
							}
						}
					}
				}
				# add "removed" hits (because there are better, e.g. non-transcript if there are transcript matches) 
				# with ambiguity and gambiguity set to 0
				foreach l $remtodo {
					lset l $ambpos 0
					if {$addumis} {
						lset l $umicountpos $umicount
					}
					lappend l 0
					puts $o [join $l \t]
				}
			}
			# set todo {}
			# set prevread $read_id
			unset -nocomplain todoa
			set prevgrouper $grouper
		}
		lappend todoa($read_id) $line
	}

	close $o

	putslog "write (bulk) isoform_counts -> isoform_counts-${root}.tsv"
	# write (bulk) isoform_counts
	set o [open isoform_counts-${root}.tsv.temp w]
	foreach refisofile $isofiles {
		set f [open $refisofile]
		set header [tsv_open $f comments]
		if {[gets $f line] == -1} {
			close $f
		} else {
			close $f
			break
		}
	}
	if {$addumis} {
		set temp [lrange $header end-5 end]
		set temp [list_regsub ^ $temp umi]
		puts $o $comments[join [list {*}$header {*}$temp] \t]
	} else {
		puts $o $comments[join $header \t]
	}
	unset -nocomplain donea
	foreach isofile $isofiles {
		set f [open $isofile]
		set cheader [tsv_open $f comments]
		if {$cheader ne $header} {
			if {[gets $f line] == -1} {
				# different header because empty file -> just skip
				continue
			}
			error "header of $isofile differs from header of $refisofile"
		}
		set isopos [lsearch $cheader transcript]
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set iso [lindex $line $isopos]
			if {![info exists tcounta($iso,t)]} continue
			if {$addumis} {
				puts $o [join [lrange $line 0 end-6] \t]\t[format %.2f [get tcounta($iso,t) 0]]\t[get tcounta($iso,u) 0]\t[get tcounta($iso,s) 0]\t[format %.2f [get tcounta($iso,a) 0]]\t[get tcounta($iso,au) 0]\t[get tcounta($iso,as) 0]\t[format %.2f [get utcounta($iso,t) 0]]\t[get utcounta($iso,u) 0]\t[get utcounta($iso,s) 0]\t[format %.2f [get utcounta($iso,a) 0]]\t[get utcounta($iso,au) 0]\t[get utcounta($iso,as) 0]
			} else {
				puts $o [join [lrange $line 0 end-6] \t]\t[format %.2f [get tcounta($iso,t) 0]]\t[get tcounta($iso,u) 0]\t[get tcounta($iso,s) 0]\t[format %.2f [get tcounta($iso,a) 0]]\t[get tcounta($iso,au) 0]\t[get tcounta($iso,as) 0]
			}
		}
		close $f
	}
	close $o
	cg select -s - isoform_counts-${root}.tsv.temp isoform_counts-${root}.tsv.temp2
	file rename -force isoform_counts-${root}.tsv.temp2 isoform_counts-${root}.tsv
	file delete isoform_counts-${root}.tsv.temp

	putslog "write (bulk) gene_counts -> gene_counts-${root}.tsv"
	# write (bulk) gene_counts
	set o [open gene_counts-${root}.tsv.temp w]
	set f [open [lindex $genefiles 0]]
	set header [tsv_open $f comments]
	close $f
	if {$addumis} {
		puts $o $comments[join $header \t]\tnicounts-$root\tumicounts-$root\tuminicounts-$root
	} else {
		puts $o $comments[join $header \t]\tnicounts-$root
	}
	unset -nocomplain donea
	foreach genefile $genefiles {
		set f [open $genefile]
		set cheader [tsv_open $f comments]
		if {$cheader ne $header} {
			error "header of $genefile differs from header of [lindex $genefiles 0]"
		}
		set genepos [lsearch $cheader geneid]
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set gene [lindex $line $genepos]
			if {![info exists gcounta($gene,w)] && ![info exists gcounta($gene,i)]} continue
			if {$addumis} {
				puts $o [join [lrange $line 0 end-1] \t]\t[format %.2f [get gcounta($gene,i) 0]]\t[format %.2f [get gcounta($gene,w) 0]]\t[format %.2f [get ugcounta($gene,i) 0]]\t[format %.2f [get ugcounta($gene,w) 0]]
			} else {
				puts $o [join [lrange $line 0 end-1] \t]\t[format %.2f [get gcounta($gene,i) 0]]\t[format %.2f [get gcounta($gene,w) 0]]
			}
		}
		close $f
	}
	close $o
	cg select -s - gene_counts-${root}.tsv.temp gene_counts-${root}.tsv.temp2
	file rename -force gene_counts-${root}.tsv.temp2 gene_counts-${root}.tsv
	file delete gene_counts-${root}.tsv.temp

	putslog "sort read_assignment file -> read_assignments-${root}.tsv.zst"
	# sort read_assignment file
	cg select -s - $tempreads2 read_assignments-${root}.tsv.temp.zst
	file rename -force read_assignments-${root}.tsv.temp.zst read_assignments-${root}.tsv.zst

	putslog "totalcounts totalcounts-${root}.tsv"
	# totalcounts
	cg select -overwrite 1 -g all -gc {sum(count*-*)} isoform_counts-${root}.tsv | cg select -rf all > totalcounts-${root}.tsv.temp
	file rename -force totalcounts-${root}.tsv.temp totalcounts-${root}.tsv
	cg select -overwrite 1 -g all -gc {sum(umicount*-*)} isoform_counts-${root}.tsv | cg select -rf all > totalumicounts-${root}.tsv.temp
	file rename -force totalumicounts-${root}.tsv.temp totalumicounts-${root}.tsv
}

proc iso_isoquant_sc_ambigcount {ambiguity gambiguity umicount} {
	global iso_isoquant_cachea
	if {[info exists iso_isoquant_cachea($ambiguity,$gambiguity,$umicount)]} {
		return $iso_isoquant_cachea($ambiguity,$gambiguity,$umicount)
	}
	if {$ambiguity > 1} {
		if {$umicount > 1} {
			set ambigcount [expr {1.0/($ambiguity * $umicount)}]
		} else {
			set ambigcount [expr {1.0/$ambiguity}]
		}
	} else {
		if {$umicount > 1} {
			set ambigcount [expr {1.0/$umicount}]
		} else {
			set ambigcount 1
		}
	}
	if {$gambiguity > 1} {
		set gambigcount [expr {$ambigcount/double($gambiguity)}]
	} else {
		set gambigcount $ambigcount
	}
	set result [list $ambigcount $gambigcount]
	set iso_isoquant_cachea($ambiguity,$gambiguity,$umicount) $result
	return $result
}

proc iso_isoquant_sc_counts {genefile isofile readfile target target2 strictpct reads_per_cell_file} {

	#
	# get cellbarcodes used (and cell counts) from reads_per_cell_file
	set useintronic 1
	set unstranded 0
	putslog "get cellbarcodes from $reads_per_cell_file"
	set temp [lrange [string trim [file_read $reads_per_cell_file]] 2 end]
	unset -nocomplain cellsa
	array set allcellsa $temp
	set cells [list_unmerge $temp]
	#
	# see which novel are intronic (and will be added to "parent" gene)
	# put all in todo, put intronic in parenta (says in which genes intron)
	putslog "match intronic reads/novel genes to parent gene"
	unset -nocomplain matcha
	unset -nocomplain parenta
	catch {close $f}
	set f [gzopen $genefile]
	set gheader [tsv_open $f comments]
	set name [file tail [file dir $genefile]]
	set gposs [lrange [tsv_basicfields $gheader 3] 1 end]
	lappend gposs [lsearch $gheader strand]
	set geneidpos [lsearch $gheader geneid]
	lappend gposs $geneidpos
	set strandpos [lindex $gposs 2]
	set todo {}
	unset -nocomplain tocheck
	set tocheck(+) {}
	set tocheck(-) {}
	set tocheck(.) {}
	while 1 {
		set read [gets $f line]
		if {$read == -1} break
		set sline [split $line \t]
		foreach {begin end strand geneid} [list_sub $sline $gposs] break
		lappend todo $sline
		# match unstranded here for intronic read addition
		# add "novel" intronic from each strand
		if {$unstranded} {
			lset sline $strandpos .
		}
		set strand .
		if {[regexp ^novelg_ $geneid]} {
			set remove {}
			set pos [llength $tocheck($strand)]
			set match 0
			while {[incr pos -1] >= 0} {
				set cline [lindex $tocheck($strand) $pos]
				foreach {cbegin cend cstrand cgeneid} [list_sub $cline $gposs] break
				if {$begin >= $cbegin && $end <= $cend} {
					set match 1
					set temp [lindex $sline $geneidpos]
					lappend matcha($cgeneid) $temp
					set parenta($temp) $cgeneid
					break
				}
				if {$begin >= $cend} {lappend remove $cline}
			}
			if {[llength $remove]} {set tocheck($strand) [list_lremove $tocheck($strand) $remove]}
			# don't add to tocheck, only checking against known genes
		} else {
			lappend tocheck($strand) $sline
		}
	}
	close $f

	#
	# get counts from readfile, load in memory
	unset -nocomplain cellsa
	unset -nocomplain genea
	unset -nocomplain ctcounta
	unset -nocomplain cgcounta
	unset -nocomplain ireadsa
	unset -nocomplain ireadstargeta
	if {$useintronic} {
		set gpos 0
		set glen [llength $todo]
		set tocheck(+) {}
		set tocheck(-) {}
	}
	catch {close $f}
	set f [gzopen $readfile]
	set header [tsv_open $f]
	set poss [list_cor $header {chromosome begin end strand isoform_id gene_id ambiguity inconsistency assignment_type covered_pct polya cellbarcode umi umicount gambiguity}]
	set readidpos [lsearch $header read_id]
	putslog "get counts from $readfile"
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		# set read_id [lindex $line $readidpos]
		foreach {chr begin end strand isoform_id gene_id ambiguity inconsistency assignment_type covered_pct polya cellbarcode umi umicount gambiguity} [list_sub $line $poss] break
		if {$gambiguity == 0} continue
		if {![info exists allcellsa($cellbarcode)]} {
			continue
		}
		incr cellsa($cellbarcode)
		foreach {ambigcount gambigcount} [iso_isoquant_sc_ambigcount $ambiguity $gambiguity $umicount] break
		if {$isoform_id ne "."} {
			if {$inconsistency < 2} {
				iso_isoquant_add_tcounts ctcounta $isoform_id,$cellbarcode $assignment_type $ambiguity $polya $covered_pct $strictpct $ambigcount $umicount
			}
			if {![info exists cgcounta($gene_id,$cellbarcode,w)]} {
				lappend genea($gene_id) $cellbarcode
			}
			iso_isoquant_add_gcount cgcounta $gene_id,$cellbarcode,w $gambigcount
			if {$useintronic} {
				if {[info exists parenta($gene_id)]} {
					set usegene $parenta($gene_id)
				} else {
					set usegene $gene_id
				}
				if {![info exists cgcounta($usegene,$cellbarcode,i)]} {
					lappend genea($usegene) $cellbarcode
				}
				iso_isoquant_add_gcount cgcounta $usegene,$cellbarcode,i $gambigcount
				iso_isoquant_add_gcount cgcounta $usegene,$cellbarcode,m $ambigcount
				if {$gambiguity <= 1} {
					iso_isoquant_add_gcount cgcounta $usegene,$cellbarcode,u $ambigcount
				}
			}
		} elseif {$useintronic} {
			if {$gene_id ne "."} {
				if {![info exists cgcounta($gene_id,$cellbarcode,w)]} {
					lappend genea($gene_id) $cellbarcode
				}
				iso_isoquant_add_gcount cgcounta $gene_id,$cellbarcode,i $gambigcount
				iso_isoquant_add_gcount cgcounta $gene_id,$cellbarcode,m $ambigcount
				if {$gambiguity <= 1} {
					iso_isoquant_add_gcount cgcounta $gene_id,$cellbarcode,u $ambigcount
				}
			}
		}
	}
	close $f

	#
	putslog "write sc_gene_counts"
	analysisinfo_write $genefile $target
	# write sc_gene_counts
	catch {gzclose $o2}
	set o2 [wgzopen $target.temp.zst]
	# end-2 to remove the counts
	puts $o2 [join [lrange $gheader 0 end-2] \t]\tcell\tcount\tnicount\tmaxcount\tuniquecount
	# set rcells [list_common $cells [array names cellsa]]
	set len [llength $todo]
	set num 0
	foreach line $todo {
		# puts "[incr num]/$len $name"
		set gene [lindex $line $geneidpos]
		if {![info exists genea($gene)]} continue
		set rcells [list_remdup $genea($gene)]
		# end-2 to remove the counts
		set temp [join [lrange $line 0 end-2] \t]
		# set resultline $temp
		foreach cell $rcells {
			# default count includes introns
			if {[info exists cgcounta($gene,$cell,i)]} {
				set count $cgcounta($gene,$cell,i)
			} else {
				set count 0
			}
			# nicount without the introns
			if {[info exists cgcounta($gene,$cell,w)]} {
				set nicount $cgcounta($gene,$cell,w)
			} else {
				set nicount 0
			}
			if {[info exists cgcounta($gene,$cell,m)]} {
				set maxcount $cgcounta($gene,$cell,m)
			} else {
				set maxcount 0
			}
			if {[info exists cgcounta($gene,$cell,u)]} {
				set ucount $cgcounta($gene,$cell,u)
			} else {
				set ucount 0
			}
			# append resultline \t$count
			if {$count > 0 || $nicount > 0 || $maxcount > 0} {
				set count [format %.2f $count]
				set nicount [format %.2f $nicount]
				puts $o2 $temp\t$cell\t$count\t$nicount\t$maxcount\t$ucount
			}
		}
	}
	# close $o
	gzclose $o2
	cg select -overwrite 1 -s - $target.temp.zst $target.temp2.zst
	file rename -force $target.temp2.zst $target
	file delete $target.temp.zst
	# check
	# cg select -f {{weight=1.0/($umicount * $gambiguity)}} -g gene_id -gc {sum(weight)} $readfile
	#

	# make sc transcript counts
	putslog "write sc transcript counts"
	analysisinfo_write $isofile $target2
	catch {close $f} catch {close $o}
	set f [gzopen $isofile]
	set header [tsv_open $f comments]
	set newheader [list_sub $header -exclude [list_find -regexp $header ^counts]]
	set poss [list_cor $header $newheader]
	set isopos [lsearch $header transcript]
	set genepos [lsearch $header geneid]
	# set name [file tail $isofile]
	set o [wgzopen $target2.temp.zst]
	lappend newheader cell counts_weighed counts_unique counts_strict counts_aweighed counts_aunique counts_astrict
	puts $o [join $newheader \t]
	set num 0
	while 1 {
		if {[gets $f line] == -1} break
		set sline [split $line \t]
		set sline [list_sub $sline $poss]
		set iso [lindex $sline $isopos]
		set gene [lindex $sline $genepos]
		if {![info exists genea($gene)]} continue
		set rcells [list_remdup $genea($gene)]
		# puts "[incr num] $iso"
		foreach cell $rcells {
			if {![info exists ctcounta($iso,$cell,t)] || $ctcounta($iso,$cell,t) == 0} continue
			puts $o [join $sline \t]\t$cell\t[format %.2f $ctcounta($iso,$cell,t)]\t[get ctcounta($iso,$cell,u) 0]\t[get ctcounta($iso,$cell,s) 0]\t[format %.2f [get ctcounta($iso,$cell,a) 0]]\t[get ctcounta($iso,$cell,au) 0]\t[get ctcounta($iso,$cell,as) 0]
		}
	}
	gzclose $o
	gzclose $f
	file rename -force $target2.temp.zst $target2
}

proc iso_isoquant_job {args} {
	# putslog [list iso_isoquant_job {*}$args]
	upvar job_logdir job_logdir
	global appdir
	set cmdline [clean_cmdline cg iso_isoquant {*}$args]
	set preset nanopore
	set distrreg chr
	set refseq {}
	set skips {}
	set reftranscripts {}
	set sqanti 1
	set threads 8
	set data_type nanopore
	set quantification all
	set gene_quantification {}
	set options {}
	set strictpct 90
	set analysisname isoquant
	set datatype nanopore
	set model_construction_strategy default_ont
	set splice_correction_strategy default_ont
	set cleanup 1
	set resultfile {}
	set regions {}
	# set skipregions {chrM M}
	set organelles {}
	set skipregions {}
	set options {}
	set singlecell 0
	set addumis 0
	cg_options iso_isoquant args {
		-preset {
			set preset $value
			if {[regexp ^sc $value]} {
				set singlecell 1
				set value [string range $value 2 end]
			}
			if {$value in {nanopore ont {}}} {
			} elseif {$value in "sens ont_sens nanopore_sens"} {
				set datatype nanopore
				set splice_correction_strategy default_ont
				set model_construction_strategy sensitive_ont
			} elseif {$value in "all"} {
				set datatype nanopore
				set splice_correction_strategy default_ont
				set model_construction_strategy all
			} elseif {$value in "pacbio"} {
				set datatype pacbio
				set splice_correction_strategy default_pacbio
				set model_construction_strategy default_pacbio
			} elseif {$value in "pacbiosens"} {
				set datatype pacbio
				set splice_correction_strategy default_pacbio
				set model_construction_strategy sensitive_pacbio
			} elseif {$value in "pacbioall"} {
				set datatype pacbio
				set splice_correction_strategy default_pacbio
				set model_construction_strategy all
			} elseif {$value in "assembly"} {
				set datatype assembly
				set splice_correction_strategy assembly
				set model_construction_strategy assembly
			} else {
				error "isoquant preset $value not supported, must be one of: nanopore sens sc scnanopore all pacbio pacbiosens assembly"
			}
			if {$preset eq "nanopore"} {
				set analysisname isoquant
			} else {
				set analysisname isoquant_$preset
			}
		}
		-addumis {
			set addumis $value
		}
		-refseq {
			set refseq $value
		}
		-organelles {
			set organelles $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-distrreg {
			if {$value eq "regionfile"} {
				set distrreg regionfile
			} else {
				set distrreg [distrreg_checkvalue $value]
			}
		}
		-transcript_quantification - -quantification {
			set quantification $value
		}
		-gene_quantification {
			set gene_quantification $value
		}
		-data_type {
			set data_type $value
		}
		-splice_correction_strategy {
			set splice_correction_strategy $value
		}
		-model_construction_strategy {
			set model_construction_strategy $value
		}
		-matching_strategy {
			lappend options --matching_strategy $value
		}
		-threads {
			set threads $value
		}
		-skip {
			lappend skips -skip $value
		}
		-regions {
			set regions $value
		}
		-skipregions {
			set skipregions $value
		}
		-cleanup {
			if {$value ni {0 1}} {error "wrong option $value for -cleanup: must be 0 or 1"}
			set cleanup $value
		}
	} {bam resultfile} 1 2
	set bam [file_absolute $bam]
	set bamindex $bam.[indexext $bam]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root ${analysisname}-[file_rootname $bam]
		set resultfile [file dir $bam]/isoform_counts-$root.tsv
		set sample [file tail [file dir $bam]]
	} else {
		set resultfile [file_absolute $resultfile]
		set root [file_rootname $resultfile]
	}
	set resultfile [file_absolute $resultfile]
	set sampledir [file dir $resultfile]
	if {$gene_quantification eq ""} {set gene_quantification $quantification}
	set reftranscripts [get_transcriptsfile_tsv $reftranscripts $refseq]
	# hanging problems with threads, run single
	set threads 1
	# set iso_isoquantdir [findiso_isoquant]

	#
	# make job log file
	job_logfile $sampledir/iso_${root} $sampledir $cmdline \
		{*}[versions iso_isoquant dbdir zstd os]

	# analysis per sample
	if {$regions eq ""} {
		set regions [list_remove [distrreg_regs $distrreg $refseq g] unaligned]
	}
	cd $sampledir
	set sample [file tail $sampledir]
	if {$bam eq ""} {
		set bam [lindex [jobglob map-sminimap*.bam map-sminimap*.cram map-*.bam map-*.cram] 0]
	}
	if {![llength $regions]} {
		set regions {{}}
	}
	set workdir [file_absolute $root.temp]
	shadow_mkdir $workdir
	set mainskips [list \
		$resultfile \
		$sampledir/gene_counts-${root}.tsv \
		$sampledir/read_assignments-${root}.tsv.zst \
		$sampledir/totalcounts-${root}.tsv \
	]
	foreach region $regions {
		if {[regions_skip $region $skipregions]} continue
		set regdir $workdir/$root-$region
		set regionskips [list \
			$regdir/isoform_counts-${root}.tsv \
			$regdir/gene_counts-${root}.tsv \
			$regdir/read_assignments-${root}.tsv.zst \
		]
		if {$reftranscripts eq "none"} {set depreftranscripts ""} else {set depreftranscripts $reftranscripts}
		if {[regions_organelle $refseq $organelles $region]} {
			# dont need regionsskip here, as these are the direct results of iso_organelle_job
			iso_organelle_job \
				-refseq $refseq \
				-reftranscripts $depreftranscripts \
				-region $region \
				-skip $mainskips \
				$bam \
				$regdir/isoform_counts-${root}.tsv
			continue
		}
		job isoquant-${root}-$region -mem 6G -cores $threads -skip $mainskips -skip $regionskips -deps {
			$bam $bamindex $refseq $depreftranscripts
		} -targets {
			$regdir/00_regali
		} -vars {
			bam refseq regdir region reftranscripts threads root sample workdir
			options data_type quantification gene_quantification analysisname distrreg
			model_construction_strategy splice_correction_strategy matching_strategy
		} -code {
			if {[file exists $regdir.temp]} {file delete -force $regdir.temp}
			mkdir $regdir.temp
			analysisinfo_write $bam $regdir.temp/00_regali \
				analysis $root sample $sample \
				isocaller_reftranscripts [file tail $reftranscripts] \
				isocaller_distrreg $distrreg \
				isocaller_data_type $data_type \
				isocaller_model_construction_strategy $model_construction_strategy \
				isocaller_splice_correction_strategy $splice_correction_strategy \
				isocaller_transcript_quantification $quantification \
				isocaller_gene_quantification $gene_quantification \
				isocaller_reference $refseq \
				isocaller $analysisname isocaller_version [version isoquant]
			# region bamfile
			set tempbam [tempfile].bam
			# set tempbam $regdir.temp/regali.bam
			if {$region ne ""} {
				set samregions [samregions $region $refseq]
				catch_exec samtools view -h -b -1 $bam {*}$samregions > $tempbam
			} else {
				set samregions {}
				# use absolute link, as this may be in a shadowdir (testing version $regdir.temp/regali.bam)
				# and relative links caanot cross that
				mklink $bam $tempbam 1
			}
			if {![catch {catch_exec samtools view $tempbam | head -1} out]} {
				# only one read aligned -> skip running isoquant
				file mkdir $regdir.temp
				file mkdir $regdir.temp/00_regali
				file_write $regdir.temp/not_enough_reads ""
			} else {
				catch_exec samtools index $tempbam
				# region gene file
				if {$reftranscripts ne "none"} {
					convert_isoquant_reggenedb $reftranscripts $samregions $refseq regreftranscripts reggenedb
					set tempgenedb [tempfile].db
					set emptyref [tsv_empty $reggenedb]
					if {$emptyref} {
						set tempfile [tempfile]
						file_write $tempfile [deindent {
							name	gene	chromosome	begin	end	strand	exonStarts	exonEnds	exonCount	source	transcript_name	gene_id	gene_name	name2
							dummyname	dummygene	dummychr	1	1	+	1	1	1	dummysource	dummysource	dummytranscript	dummygene	dummyname2
						}]\n
						cg_tsv2gtf -genecol gene_id -addgene 1 $tempfile $reggenedb
					}
					exec isoquant3_gtf2db --complete_genedb --input $reggenedb --output $tempgenedb
					lappend options --genedb $tempgenedb
				}
				file delete -force $regdir.temp/.params $regdir.temp/00_regali $regdir.temp/OUT $regdir.temp/isoquant.log 
				exec isoquant3 \
					--data_type $data_type \
					--model_construction_strategy $model_construction_strategy \
					--splice_correction_strategy $splice_correction_strategy \
					--transcript_quantification $quantification \
					--gene_quantification $gene_quantification \
					--threads $threads \
					--reference $refseq \
					--bam $tempbam \
					--keep_tmp \
					{*}$options \
					-o $regdir.temp 2>@ stderr >@ stdout
				file delete $tempbam
			}
			if {![file exists $regdir.temp/00_regali] && [file exists $regdir.temp/OUT]} {
				file rename $regdir.temp/OUT $regdir.temp/00_regali
			}
			file delete -force $regdir
			file rename $regdir.temp $regdir
		}
		job isquant_convert-${root}-$region -cores 1 -mem 5g -skip $mainskips -deps {
			$regdir/00_regali $refseq $depreftranscripts
		} -targets {
			$regdir/isoform_counts-${root}.tsv
			$regdir/gene_counts-${root}.tsv
			$regdir/read_assignments-${root}.tsv.zst
		} -vars {
			sample refseq regdir region reftranscripts root analysisname cleanup singlecell
		} -code {
			analysisinfo_write $regdir/00_regali $regdir/isoform_counts-${root}.tsv
			analysisinfo_write $regdir/00_regali $regdir/gene_counts-${root}.tsv
			analysisinfo_write $regdir/00_regali $regdir/read_assignments-${root}.tsv
			set isodir $regdir/00_regali
			set destdir $regdir
			if {[file exists $regdir/not_enough_reads]} {
				set header {chromosome begin end strand exonStarts exonEnds transcript gene geneid}
				lappend header category size \
					counts_iqall-$root counts_weighed-$root counts_unique-$root counts_strict-$root \
					counts_aweighed-$root counts_aunique-$root counts_astrict-$root
				file_write $regdir/isoform_counts-${root}.tsv [join $header \t]\n
				file_write $regdir/gene_counts-${root}.tsv \
					[join [list chromosome begin end strand gene geneid counts-${root}] \t]\n
				file_write $regdir/read_assignments-${root}.tsv [join {
					read_id chromosome begin end strand exonStarts exonEnds aligned_size
					ambiguity inconsistency covered_pct polya classification closest_known
				} \t]\n
				cg zst $regdir/read_assignments-${root}.tsv
			} else {
				if {$reftranscripts ne "none"} {
					set samregions [samregions $region $refseq]
					convert_isoquant_reggenedb $reftranscripts $samregions $refseq regreftranscripts reggenedb
					convert_isoquant $isodir $destdir $sample $refseq $reggenedb $regreftranscripts $root $singlecell
				} else {
					convert_isoquant $isodir $destdir $sample $refseq "" "" $root $singlecell
				}
			}
			if {$cleanup} {
				cg_shadow_delete $isodir
			}
		}
	}

	#
	# combine results from separately calculated regions
	set isofiles {}
	set genefiles {}
	set readfiles {}
	set missing {}
	foreach region $regions {
		if {[regions_skip $region $skipregions]} continue
		set regdir $workdir/$root-$region
		set isofile $regdir/isoform_counts-${root}.tsv
		set genefile $regdir/gene_counts-${root}.tsv
		set readfile $regdir/read_assignments-${root}.tsv.zst
		lappend isofiles $isofile
		lappend genefiles $genefile
		lappend readfiles $readfile
		if {![file exists $isofile]} {
			puts "$isofile missing"
			lappend missing $isofile
		}
		if {![file exists $genefile]} {
			puts "$genefile missing"
			lappend missing $genefile
		}
		if {![file exists $readfile]} {
			puts "$readfile missing"
			lappend missing $readfile
		}
	}
	job isoquant_join-$root -cores 1 -deps [list {*}$isofiles {*}$genefiles {*}$readfiles] -targets {
		$resultfile
		gene_counts-${root}.tsv
		read_assignments-${root}.tsv.zst
		totalcounts-${root}.tsv
	} -vars {
		isofiles genefiles readfiles sample refseq regdir region reftranscripts genedb root analysisname
		strictpct workdir singlecell resultfile skipregions addumis
	} -code {
		analysisinfo_write [lindex $isofiles 0] $resultfile isocaller_skipregions $skipregions
		analysisinfo_write [lindex $genefiles 0] gene_counts-${root}.tsv isocaller_skipregions $skipregions
		analysisinfo_write [lindex $readfiles 0] read_assignments-${root}.tsv isocaller_skipregions $skipregions
		iso_isoquant_mergeresults $isofiles $genefiles $readfiles $strictpct $sample $root $analysisname $addumis
		cg tsv2gtf $resultfile isoforms-${root}.gtf
		# if {!$singlecell} {catch {shadow_delete $workdir}}
	}
	if {$singlecell} {
		# parallelize sc analysis
		# will redistribute the merge at bulk level (which added e.g. gambiguity)
		mkdir $workdir/gene_counts
		set sc_gene_counts_files [distrreg_job -refseq $refseq \
			-skip [list sc_gene_counts_raw-${root}.tsv.zst sc_isoform_counts_raw-${root}.tsv.zst] \
			gene_counts-${root}.tsv $workdir/gene_counts/gene_counts- .tsv $regions \
		]
		mkdir $workdir/isoform_counts
		set sc_iso_counts_files [distrreg_job -refseq $refseq \
			-skip [list sc_gene_counts_raw-${root}.tsv.zst sc_isoform_counts_raw-${root}.tsv.zst] \
			$resultfile $workdir/isoform_counts/isoform_counts- .tsv $regions \
		]
		mkdir $workdir/readassignments
		# distributed readassignments are only made for making the sc_gene and sc_isoforms results
		# so those are the skip targets here
		distrreg_job \
			-skip [list sc_gene_counts_raw-${root}.tsv.zst sc_isoform_counts_raw-${root}.tsv.zst] \
			-refseq $refseq \
			read_assignments-${root}.tsv.zst $workdir/readassignments/readassignments- .tsv $regions
		set sc_gene_counts_files {}
		set sc_iso_counts_files {}
		foreach region $regions {
			if {[regions_skip $region $skipregions]} continue
			set readfile $workdir/readassignments/readassignments-$region.tsv
			set genefile $workdir/gene_counts/gene_counts-$region.tsv
			set isofile $workdir/isoform_counts/isoform_counts-$region.tsv
			mkdir $workdir/sc_gene_counts_long
			mkdir $workdir/sc_isoform_counts_long
			set target $workdir/sc_gene_counts_long/sc_gene_counts_long-$root-$region.tsv.zst
			set target2 $workdir/sc_isoform_counts_long/sc_isoform_counts_long-$root-$region.tsv.zst
			lappend sc_gene_counts_files $target
			lappend sc_iso_counts_files $target2
			job isquant_sc_count-$region-$root \
			-skip [list sc_gene_counts_raw-${root}.tsv.zst sc_isoform_counts_raw-${root}.tsv.zst] \
			-deps {
				$readfile $genefile $sampledir/reads_per_cell_raw.tsv
			} -targets {
				$target $target2
			} -vars {
				readfile genefile strictpct isofile sampledir
			} -code {
				set reads_per_cell_file $sampledir/reads_per_cell_raw.tsv
				iso_isoquant_sc_counts $genefile $isofile $readfile $target $target2 $strictpct $reads_per_cell_file
			}
		}
		job isoquant_sc_join_gene-$root -cores 1 \
		-deps [list {*}$sc_gene_counts_files {*}$sc_iso_counts_files] \
		-targets {
			sc_gene_counts_raw-${root}.tsv.zst
		} -vars {
			sc_gene_counts_files root
		} -code {
			analysisinfo_write [lindex $sc_gene_counts_files 0] sc_gene_counts_raw-${root}.tsv
			cg cat {*}$sc_gene_counts_files | cg zst > sc_gene_counts_raw-${root}.tsv.zst.temp
			file rename sc_gene_counts_raw-${root}.tsv.zst.temp sc_gene_counts_raw-${root}.tsv.zst
		}
		set sc_result sc_isoform_counts_raw-${root}.tsv.zst
		job isoquant_sc_join_iso-$root -cores 1 -deps $sc_iso_counts_files -targets {
			$sc_result
		} -vars {
			sc_iso_counts_files root sc_result
		} -code {
			analysisinfo_write [lindex $sc_iso_counts_files 0] $sc_result
			cg cat {*}$sc_iso_counts_files | cg zst > $sc_result.temp
			file rename $sc_result.temp $sc_result
		}
	}
}

proc cg_iso_isoquant {args} {
	set args [job_init {*}$args]
	iso_isoquant_job {*}$args
	job_wait
}


