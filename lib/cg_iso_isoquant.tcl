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
	set f [open $regreftranscripts]
	tsv_open $f
	set read [gets $f line]
	close $f
	if {$read == -1} {
		set header [cg select -h $reftranscripts]
		catch {close $o}
		set o [open $regreftranscripts w]
		puts $o [join $header \t]
		set line [list_fill [llength $header] {}]
		set poss [list_sub [tsv_basicfields $header 9 0] {0 1 2 6 7 8}]
		foreach pos $poss value {chr1 1 2 + 1, 2,} {
			if {$pos == -1} continue
			lset line $pos $value
		}
		foreach {field value} {
			name dummy1 gene dummy2 cdsStart 1 cdsEnd 1 exonCount 1 
			source dummy3 transcript_name dummy4 gene_id dummy5 gene_name dummy6
			name2 dummyname2
		} {
			set pos [lsearch $header $field]
			if {$pos == -1} continue
			lset line $pos $value
		}
		puts $o [join $line \t]
		close $o
	}
	set genecol [lindex [list_common {gene_id geneid gene} [cg select -h $regreftranscripts]] 0]
	cg_tsv2gtf -genecol $genecol -addgene 1 $regreftranscripts $reggenedb
	set tempgenedb [tempfile].db
}

proc convert_isoquant_add {varVar {count 1}} {
	upvar $varVar var
	if {![info exists var]} {
		set var $count
	} else {
		set var [expr {$var + $count}]
	}
}

proc convert_isoquant {isodir destdir sample refseq reggenedb regreftranscripts root singlecell} {
	set strictpct 90
	set read_assignmentsfile [gzfile $isodir/*.read_assignments.tsv]
	set targetisoformcountsfile $destdir/isoform_counts-${root}.tsv
	set targetgenecountsfile $destdir/gene_counts-${root}.tsv
	set targetreadassignmentsfile $destdir/read_assignments-${root}.tsv
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
				if {[lindex $genebasica($geneid) 0] ne "$chr"} {
					puts stderr "Chromosome does not match for transcripts of same gene $geneid: $line"
				} elseif {[lindex $genebasica($geneid) 3] ne "$strand"} {
					puts stderr "Strand does not match for transcripts of same gene $geneid: $line"
				} elseif {$begin > [lindex $genebasica($geneid) 2]} {
					puts stderr "transcript not overlapping previous for same gene $geneid: $line"
				} elseif {$end > [lindex $genebasica($geneid) 2]} {
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
	}

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
		unset -nocomplain converta
		unset -nocomplain outputa
		cg_gtf2tsv $mfile $destdir/transcripts_models-$sample.tsv
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

	#
	# umi counts
	# ----------
	unset -nocomplain umicounta
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
		if {[regexp {([A-Z]+)_([A-Z]+)#} $read temp cellbarcode umi]} {
			incr umicounta($umi)
		}
		if {[gets $f line] == -1} break
	}


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
			if {[info exists outputa($transcript)]} {
				lappend read2isoa($read) $transcript
			}
		}
		gzclose $f
	}

	# check which reads have ambiguous mappings (from original read_assignmentsfile)
	unset -nocomplain ambiga
	array set ambiga [split [cg select -hc 1 -g read_id $read_assignmentsfile | cg select -q {$count > 1}] \n\t]

	#
	# make read_assignments file
	# --------------------------
	# set temptarget $destdir/read_assignments-${root}.tsv.temp
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
	set genepos [lsearch $newheader gene_id]
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
	unset -nocomplain ctcounta
	unset -nocomplain cgcounta
	unset -nocomplain donea
	for {set p 0} {$p <= 100} {incr p} {set pcta($p) 0}
	while 1 {
		set line [split $line \t]
		set exons [lindex $line $exonspos]
		set line [list_sub $line $poss]
		set read [lindex $line $readpos]
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
		set gene [lindex $line $genepos]
		if {[info exists ambiga($read)]} {set ambig $ambiga($read)} else {set ambig 0}
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
			if {$assignment_type in "inconsistent unique_minor_difference intergenic" 
				|| $closest_known ne "."
			} {
				set knownmatch 0
			} else {
				if {$inconsistency > 1} {
					set knownmatch 0
				} else {
					set knownmatch 1
				}
			}
			if {!$knownmatch} {
				if {[info exists donea($read)]} {
					if {[gets $f line] == -1} break
					continue
				}
				set modelisos {}
				set knownisos {}
				foreach iso $read2isoa($read) {
					if {![info exists converta($iso)]} {
						lappend modelisos $iso
					} else {
						lappend knownisos $converta($iso)
					}
				}
				if {![llength $modelisos]} {
					# known hits only (no new model) -> keep original known hits
					# keep original readassignment if model was scrapped from output
					if {![info exists outputa($iso)]} {
						set isos [list $closest_known]
					} else {
						set isos [list $converta($iso)]
					}
					set inconsistencylist [list $inconsistency]
				} else {
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
					set ambiga($read) $ambig
					# write only model hits
					set donea($read) 1
				}
			} else {
				if {[info exists donea($read)]} {
					# new ones have been added already, proceed normally
					set isos [list $closest_known]
					set inconsistencylist [list $inconsistency]
				} else {
					set modelisos {}
					set knownisos {}
					foreach iso $read2isoa($read) {
						if {![info exists converta($iso)]} {
							lappend modelisos $iso
						} else {
							lappend knownisos $converta($iso)
						}
					}
					if {[llength $modelisos]} {
						if {$ambig == 0} {set ambig 1}
						incr ambig [llength $modelisos]
						set ambiga($read) $ambig
						lset line $assignment_typepos ambiguous
						set isos [list $closest_known {*}$modelisos]
						set inconsistencylist [list $inconsistency]
						lappend inconsistencylist {*}[list_fill [llength $modelisos] 0]
					} else {
						# no model hits -> proceed normally
						# keep original readassignment if model was scrapped from output
						if {![info exists outputa($iso)]} {
							set isos [list $closest_known]
						} else {
							set isos [list $converta($iso)]
						}
						set inconsistencylist [list $inconsistency]
					}
					set donea($read) 1
				}
			}
		} else {
			set isos [list $closest_known]
			set inconsistencylist [list $inconsistency]
		}
		set umicount [get umicounta($umi) 1]
		set ambigcount [convert_isoquant_ambigcount $ambig $umicount]
		if {$umicount > 1} {
			set count [expr {1.0/$umicount}]
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
				} else {
					set geneid [lindex $line $genepos]
					set transcript2genea($iso) $geneid
				}
			} else {
				set geneid .
			}
			if {[info exists geneconva($geneid)]} {
				set gene $geneconva($geneid)
			} elseif {[info exists geneid2genea($geneid)]} {
				set gene $geneid2genea($geneid)
			} else {
				set gene $geneid
			}
			lset line $genepos $gene
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
			convert_isoquant_add gcounta($geneid) $ambigcount
#			if {$singlecell} {
#				if {$inconsistency < 2} {
#					convert_isoquant_add ctcounta($iso,$cellbarcode,t) $ambigcount
#					if {!$ambig && $assignment_type ne "inconsistent"} {
#						# unique
#						convert_isoquant_add ctcounta($iso,$cellbarcode,u) $count
#						if {$pct >= $strictpct} {
#							# strict
#							convert_isoquant_add ctcounta($iso,$cellbarcode,s) $count
#						}
#					}
#					if {$polya eq "True"} {
#						convert_isoquant_add ctcounta($iso,$cellbarcode,a) $ambigcount
#						if {!$ambig && $assignment_type ne "inconsistent"} {
#							# unique
#							convert_isoquant_add ctcounta($iso,$cellbarcode,au) $count
#							if {$pct >= $strictpct} {
#								# strict
#								convert_isoquant_add ctcounta($iso,$cellbarcode,as) $count
#							}
#						}
#					}
#				}
#				convert_isoquant_add cgcounta($geneid,$cellbarcode) $ambigcount
#			}
		}
		if {[gets $f line] == -1} break
	}

	gzclose $o
	gzclose $f
	cg select -overwrite 1 -s - $temptarget ${temptarget}2
	file rename -force ${temptarget}2 $targetreadassignmentsfile
	file delete $temptarget

	# make isoform_counts file
	# ------------------------
	# read isoquant counts
	unset -nocomplain countsa
	set f [gzopen [gzfile $isodir/*.transcript_counts.tsv]]
	gets $f
	array set countsa [read $f]
	gzclose $f
	set mfile [gzfile $isodir/*.transcript_model_counts.tsv]
	if {[file exists $mfile]} {
		set f [gzopen $mfile]
		gets $f
		array set countsa [read $f]
		gzclose $f
	}

	set basefields {chromosome begin end strand exonStarts exonEnds transcript gene geneid}
	# lappend remove {*}$basefields
	set newheader $basefields
	if {$regreftranscripts ne ""} {
		catch {close $f}
		set f [gzopen $regreftranscripts]
		set header [tsv_open $f comments]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 11 12 13}]
		set remove [list_sub $header $poss]
		foreach {isopos genepos geneidpos} [lrange $poss 6 end] break
		lappend remove {*}[list_sub $header [lrange $poss 6 end]]
		set left [list_lremove $header $remove]
		lappend poss {*}[list_cor $header $left]
		set temp [list_common $left $basefields]
		foreach field $temp {
			set pos [lsearch $left $field]
			lset left $pos ${field}_ori
		}
		lappend newheader {*}$left
	}
	#
	if {$comments eq ""} {
		set comments [deindent {
			#filetype	tsv/transcriptsfile
			#fileversion	0.99
			#fields	table
			#fields	field	number	type	description
			#fields	name	1	String	Name of transcript (usually transcript_id from GTF)
			#fields	chromosome	1	String	Chromosome name
			#fields	strand	1	String	+ or - for strand
			#fields	begin	1	Integer	Transcription start position
			#fields	end	1	Integer	Transcription end position
			#fields	exonCount	1	Integer	Number of exons
			#fields	exonStarts	E	Integer	Exon start positions
			#fields	transcript	1	String	transcript id
			#fields	gene	1	String	gene name
			#fields	geneid	1	String	gene id
			#fields	cdsStart	1	Integer	Coding region start
			#fields	cdsEnd	1	Integer	Coding region end
			#fields	type	1	String	type of element
			#fields	exonEnds	E	Integer	Exon end positions
			#fields	exonFrames	E	Integer	Exon frame offsets {0,1,2}
			#fields	score	1	Float	Score
			#fields	name2	1	String	Alternate name (e.g. gene_id from GTF)
			#fields	length	1	Integer	isoform length
			#fields	exons	1	Integer	Number of exons
			#fields	category	1	String	one of the isoform categories (known, novel_in_catalog, novel_not_in_catalog, intergenic)
			#fields	associated_gene	1	String	the reference gene name
			#fields	associated_transcript	1	String	the reference transcript name
			#fields	counts	1	Integer	Number of reads mapping to isoform
		}]
	}
	catch {close $o}
	set o [open $targetisoformcountsfile.temp w]
	puts $o $comments
	puts $o [join $newheader \t]\tcategory\tsize\tcounts_iqall-$sample\tcounts_weighed-${root}\tcounts_unique-${root}\tcounts_strict-${root}\tcounts_aweighed-${root}\tcounts_aunique-${root}\tcounts_astrict-${root}
	if {$regreftranscripts ne ""} {
		set keepoutputa [array get outputa]
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			set iso [lindex $line $isopos]
			if {$iso eq ""} continue
			if {[info exists outputa($iso)]} {
				unset outputa($iso)
			}
			set t [get tcounta($iso,t) 0]
			if {[get countsa($iso) 0] == 0 && $t == 0} continue
			set geneid [lindex $line $geneidpos]
			if {[info exists geneid2genea($geneid)]} {
				lset line $genepos $geneid2genea($geneid)
			}
			set line [list_sub $line $poss]
			puts $o [join $line \t]\tknown\t$sizea($iso)\t[get countsa($iso)]\t[format %.2f $t]\t[get tcounta($iso,u) 0]\t[get tcounta($iso,s) 0]\t[format %.2f [get tcounta($iso,a) 0]]\t[get tcounta($iso,au) 0]\t[get tcounta($iso,as) 0]
		}
		catch {close $f}
	}
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
			if {[info exists geneid2genea($geneid)]} {
				lset line $genepos $geneid2genea($geneid)
			}
			if {![info exists outputa($oriname)]} continue
			set t [get tcounta($oriname,t) 0]
			if {[get countsa($oriname) 0] < 1 && $t < 1} continue
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
			}
			set line [list_sub $line $poss]
			# set genebasica($geneid) [list_sub $line {0 1 2 3}]
			if {![info exists genebasica($geneid)]} {
				set genebasica($geneid) [list_sub $line {0 1 2 3}]
			}
			puts $o [join $line \t]\t$category\t$sizea($oriname)\t[get countsa($oriname)]\t[format %.2f $t]\t[get tcounta($oriname,u) 0]\t[get tcounta($oriname,s) 0]\t[format %.2f [get tcounta($oriname,a) 0]]\t[get tcounta($oriname,au) 0]\t[get tcounta($oriname,as) 0]
		}
		close $f
	}
	close $o
	cg select -overwrite 1 -s - $targetisoformcountsfile.temp $targetisoformcountsfile.temp2
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

proc iso_isoquant_mergeresults {isofiles genefiles readfiles strictpct sample root {analysisname isoquant} {singlecell 0}} {
	set tempreads [tempfile]
	set tempreads2 [tempfile]
	cg cat -m 1 {*}$readfiles | cg select -s read_id > $tempreads

	unset -nocomplain tcounta
	unset -nocomplain gcounta
	unset -nocomplain ctcounta
	unset -nocomplain cgcounta
	unset -nocomplain cella
	catch {close $f} ; catch {close $o}
	set f [open $tempreads]
	set header [tsv_open $f comment]
	set poss [list_cor $header {assignment_type isoform_id covered_pct polya cellbarcode umicount}]
	set readpos [lsearch $header read_id]
	set inconsistencypos [lsearch $header inconsistency]
	set isopos [lsearch $header isoform_id]
	set genepos [lsearch $header gene_id]
	set ambpos [lsearch $header ambiguity]
	set pctpos [lsearch $header covered_pct]
	set o [wgzopen $tempreads2]
	puts $o $comment[join $header \t]\tgambiguity
	set prevread {}
	set todo {}
	set read 0
	set nr 0

	# gather counts in memory (tcounta,gcounta), can get very big for singlecell (ctcounta,cgcounta)
	while 1 {
		incr nr
		if {![expr $nr%1000]} {puts $nr}
		if {$read == -1} break
		set read [gets $f line]
		set line [split $line \t]
		set read_id [lindex $line $readpos]
		if {$read_id eq ""} continue
		if {$read_id ne $prevread} {
			set ambiguity [llength $todo]
			if {$ambiguity > 1} {
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
				}
				set genes [list_remdup [list_subindex $todo $genepos]]
				set gambiguity [llength $genes]
				set ambigcount [expr {1.0/$ambiguity}]
				set gambigcount [expr {1.0/$gambiguity}]
			} elseif {$ambiguity == 1} {
				set ambigcount 1
				set gambigcount 1
				set gambiguity 1
				set genes [list_remdup [list_subindex $todo $genepos]]
				set inconsistency [lindex $todo 0 $inconsistencypos]
			}
			if {$ambiguity > 0} {
				foreach gene $genes {
					if {![info exists gcounta($gene)]} {
						set gcounta($gene) $gambigcount
					} else {
						set gcounta($gene) [expr {$gcounta($gene) + $gambigcount}]
					}
				}
				foreach l $todo {
					lset l $ambpos $ambiguity
					lappend l $gambiguity
					puts $o [join $l \t]
					# count
					foreach {assignment_type iso pct polya cellbarcode umicount} [list_sub $l $poss] break
					if {$inconsistency < 2} {
						convert_isoquant_add tcounta($iso,t) $ambigcount
						if {$ambiguity == 1 && $assignment_type ne "inconsistent"} {
							# unique
							convert_isoquant_add tcounta($iso,u)
							if {$pct >= $strictpct} {
								# strict
								convert_isoquant_add tcounta($iso,s)
							}
						}
						if {$polya eq "True"} {
							convert_isoquant_add tcounta($iso,a) $ambigcount
							if {$ambiguity == 1 && $assignment_type ne "inconsistent"} {
								# unique
								convert_isoquant_add tcounta($iso,au)
								if {$pct >= $strictpct} {
									# strict
									convert_isoquant_add tcounta($iso,as)
								}
							}
						}
					}
				}
				if {$singlecell} {
					if {$umicount <=1} {
						set count 1
					} else {
						set count [expr {1.0/$umicount}]
						set gambigcount [expr {double($gambigcount)/$umicount}]
						set ambigcount [expr {double($ambigcount)/$umicount}]
					}
					incr cella($cellbarcode)
					foreach gene $genes {
						if {![info exists cgcounta($gene)]} {
							set cgcounta($gene,$cellbarcode) $gambigcount
						} else {
							set cgcounta($gene,$cellbarcode) [expr {$cgcounta($gene) + $gambigcount}]
						}
					}
					foreach l $todo {
						lset l $ambpos $ambiguity
						lappend l $gambiguity
						puts $o [join $l \t]
						# count
						foreach {assignment_type iso pct polya cellbarcode umicount} [list_sub $l $poss] break
						if {$inconsistency < 2} {
							convert_isoquant_add ctcounta($iso,$cellbarcode,t) $ambigcount
							if {$ambiguity == 1 && $assignment_type ne "inconsistent"} {
								# unique
								convert_isoquant_add ctcounta($iso,$cellbarcode,u) $count
								if {$pct >= $strictpct} {
									# strict
									convert_isoquant_add ctcounta($iso,$cellbarcode,s) $count
								}
							}
							if {$polya eq "True"} {
								convert_isoquant_add ctcounta($iso,$cellbarcode,a) $ambigcount
								if {$ambiguity == 1 && $assignment_type ne "inconsistent"} {
									# unique
									convert_isoquant_add ctcounta($iso,$cellbarcode,au) $count
									if {$pct >= $strictpct} {
										# strict
										convert_isoquant_add ctcounta($iso,$cellbarcode,as) $count
									}
								}
							}
						}
					}
				}
			}
			set todo {}
			set prevread $read_id
		}
		lappend todo $line
	}
	close $o

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
	puts $o $comments[join $header \t]
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
			set t [get tcounta($iso,t) 0]
			puts $o [join [lrange $line 0 end-6] \t]\t[format %.2f $t]\t[get tcounta($iso,u) 0]\t[get tcounta($iso,s) 0]\t[format %.2f [get tcounta($iso,a) 0]]\t[get tcounta($iso,au) 0]\t[get tcounta($iso,as) 0]
		}
		close $f
	}
	close $o
	cg select -s - isoform_counts-${root}.tsv.temp isoform_counts-${root}.tsv.temp2
	file rename -force isoform_counts-${root}.tsv.temp2 isoform_counts-${root}.tsv
	file delete isoform_counts-${root}.tsv.temp

	# write (bulk) gene_counts
	set o [open gene_counts-${root}.tsv.temp w]
	set f [open [lindex $genefiles 0]]
	set header [tsv_open $f comments]
	close $f
	puts $o $comments[join $header \t]
	unset -nocomplain donea
	foreach genefile $genefiles {
		set f [open $genefile]
		set cheader [tsv_open $f comments]
		if {$cheader ne $header} {
			error "header of $genefile differs from header of [lindex $genefiles 0]"
		}
		set genepos [lsearch $cheader gene]
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set gene [lindex $line $genepos]
			if {![info exists gcounta($gene)]} continue
			puts $o [join [lrange $line 0 end-1] \t]\t[format %.2f $gcounta($gene)]
		}
		close $f
	}
	close $o
	cg select -s - gene_counts-${root}.tsv.temp gene_counts-${root}.tsv.temp2
	file rename -force gene_counts-${root}.tsv.temp2 gene_counts-${root}.tsv
	file delete gene_counts-${root}.tsv.temp

	if {$singlecell} {
		# write sc_gene_counts
		# --------------------
		# get cells (filter out cells with < 200 reads)
		unset -nocomplain cella()
		set list {}
		foreach cell [array names cella] {
			if {$cell eq ""} continue
			if {$cella($cell) < 200} {
				unset cella($cell)
				continue
			}
			lappend list [list $cell $cella($cell)]
		}
		set list [lsort -index 1 -integer -decreasing $list]
		set cells [list_subindex $list 0]
		#
		# make sc_gene_counts
		catch {close $f}
		set f [open [lindex $genefiles 0]]
		set header [tsv_open $f comments]
		close $f
		catch {close $o} ; catch {close $o2}
		set o [wgzopen sc_gene_counts-${root}.tsv.temp.zst]
		set resultheader [join [lrange $header 0 end-1] \t]
		foreach cell $cells {
			append resultheader \tcounts-$cell
		}
		puts $o $comments$resultheader
		set o2 [wgzopen sc_gene_counts_long-${root}.tsv.temp.zt]
		puts $o2 [join [lrange $header 0 end-1] \t]\tcell\tcount
		unset -nocomplain donea
		set genefiles [bsort $genefiles]
		foreach genefile $genefiles {
			set name [file tail [file dir $genefile]]
			set f [open $genefile]
			set cheader [tsv_open $f comments]
			if {$cheader ne $header} {
				error "header of $genefile differs from header of [lindex $genefiles 0]"
			}
			set poss [lrange [tsv_basicfields $header 3] 1 end]
			lappend poss [lsearch $cheader strand]
			lappend poss [lsearch $cheader gene]
			unset -nocomplain matcha
			set todo {}
			set tocheck(+) {}
			set tocheck(-) {}
			while 1 {
				set read [gets $f line]
				if {$read == -1} break
				set sline [split $line \t]
				foreach {begin end strand gene} [list_sub $sline $poss] break
				if {[regexp ^novelg_ $gene]} {
					set remove {}
					set pos [llength $tocheck($strand)]
					set match 0
					while {[incr pos -1] >= 0} {
						set cline [lindex $tocheck($strand) $pos]
						foreach {cbegin cend cstrand cgene} [list_sub $cline $poss] break
						if {$begin >= $cbegin && $end <= $cend} {
							set match 1
							lappend matcha($cgene) [lindex $sline $genepos]
							break
						}
						if {$begin >= $cend} {lappend remove $cline}
					}
					if {[llength $remove]} {set tocheck($strand) [list_lremove $tocheck($strand) $remove]}
					if {!$match} {
						lappend todo $sline
						# lappend tocheck($strand) $sline
					}
				} else {
					lappend todo $sline
					lappend tocheck($strand) $sline
				}
			}
			close $f
			set len [llength $todo]
			set num 0
			foreach line $todo {
				puts "[incr num]/$len $name"
				set gene [lindex $line $genepos]
				if {![info exists gcounta($gene)]} continue
				set temp [join [lrange $line 0 end-1] \t]
				set resultline $temp
				set total 0
				foreach cell $cells {
					if {[info exists cgcounta($gene,$cell)]} {
						set count $cgcounta($gene,$cell)
					} else {
						set count 0
					}
					if {[info exists matcha($gene)]} {
						foreach g $matcha($gene) {
							if {[info exists cgcounta($g,$cell)]} {
								set count [expr {$count + $cgcounta($g,$cell)}]
							}
						}
					}
					set total [expr {$total + $count}]
					set count [format %.2f $count]
					append resultline \t$count
					if {$count > 0} {
						puts $o2 $temp\t$cell\t$count
					}
				}
				if {$total < 1} continue
				puts $o $resultline
			}
			flush $o
			flush $o2
		}
		close $o
		cg select -s - sc_gene_counts-${root}.tsv.temp.zst sc_gene_counts-${root}.tsv.temp2.gz
		file rename -force sc_gene_counts-${root}.tsv.temp2.gz sc_gene_counts-${root}.tsv.gz
		file delete sc_gene_counts-${root}.tsv.temp
		
	}

	# sort read_assignment file
	cg select -s - $tempreads2 read_assignments-${root}.tsv.temp.zst
	file rename -force read_assignments-${root}.tsv.temp.zst read_assignments-${root}.tsv.zst

	# totalcounts
	cg select -overwrite 1 -g all -gc {sum(count*-*)} isoform_counts-${root}.tsv | cg select -rf all > totalcounts-${root}.tsv.temp
	file rename -force totalcounts-${root}.tsv.temp totalcounts-${root}.tsv

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
	set skipregions {chrM M}
	set options {}
	set singlecell 0
	cg_options iso_isoquant args {
		-preset {
			set preset $value
			if {[regexp ^sc] $value} {
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
		-refseq {
			set refseq $value
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
			set cleanup $value
		}
	} {bam resultfile} 1 2
	set bam [file_absolute $bam]
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
	if {$reftranscripts eq ""} {
		set reftranscripts [ref_tsvtranscripts $refseq]
	} elseif {$reftranscripts eq "none"} {
	} elseif {[file extension $reftranscripts] eq ".gtf"} {
		set reftranscripts [file_absolute $reftranscripts]
		set reftranscriptstsv [file root $reftranscripts].tsv
		if {![file exists $reftranscriptstsv]} {
			cg_gtf2tsv $reftranscripts $reftranscriptstsv.temp
			cg select -overwrite 1 -s - $reftranscriptstsv.temp $reftranscriptstsv.temp2
			file delete $reftranscriptstsv.temp
			file rename -force $reftranscriptstsv.temp2 $reftranscriptstsv
		}
		set reftranscripts $reftranscriptstsv
	} else {
		set reftranscripts [file_absolute $reftranscripts]
	}
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
	set bam [lindex [jobglob map-sminimap*.bam map-*.bam] 0]
	if {![llength $regions]} {
		set regions {{}}
	}
	# set workdir [file_absolute $root.temp]
	set workdir [shadow_workdir [file_absolute $root]]
	set mainskips [list \
		isoform_counts-${root}.tsv \
		gene_counts-${root}.tsv \
		read_assignments-${root}.tsv \
		totalcounts-${root}.tsv \
	]
	foreach region $regions {
		if {[regions_skip $region $skipregions]} continue
		set regdir $workdir/$root-$region
		set regionskips [list \
			$regdir/isoform_counts-${root}.tsv \
			$regdir/gene_counts-${root}.tsv \
			$regdir/read_assignments-${root}.tsv \
		]
		if {$reftranscripts eq "none"} {set depreftranscripts ""} else {set depreftranscripts $reftranscripts}
		job isquant-${root}-$region -mem 15G -cores $threads -skip $mainskips -skip $regionskips -deps {
			$bam $bam.bai $refseq $depreftranscripts
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
				isocaller $analysisname isocaller_version [version isoquant]
			# region bamfile
			# set tempbam [tempfile].bam
			set tempbam $regdir.temp/regali.bam
			if {$region ne ""} {
				set samregions [samregions $region $refseq]
				exec samtools view -h -b -1 $bam {*}$samregions > $tempbam
			} else {
				set samregions {}
				mklink $bam $tempbam
			}
			if {![catch {exec samtools view $tempbam | head -1} out]} {
				# only one read aligned -> skip running isoquant
				file mkdir $regdir.temp
				file mkdir $regdir.temp/00_regali
				file_write $regdir.temp/not_enough_reads ""
			} else {
				exec samtools index $tempbam
				# region gene file
				if {$reftranscripts ne "none"} {
					convert_isoquant_reggenedb $reftranscripts $samregions $refseq regreftranscripts reggenedb
					set tempgenedb [tempfile].db
					exec isoquant3_gtf2db --complete_genedb --input $reggenedb --output $tempgenedb
					lappend options --genedb $tempgenedb
				}
				file delete -force $regdir.temp/.params $regdir.temp/00_regali $regdir.temp/isoquant.log 
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
			file delete -force $regdir
			file rename $regdir.temp $regdir
		}
		job isquant_convert-${root}-$region -cores 1 -mem 10g -skip $mainskips -deps {
			$regdir/00_regali $refseq $depreftranscripts
		} -targets {
			$regdir/isoform_counts-${root}.tsv
			$regdir/gene_counts-${root}.tsv
			$regdir/read_assignments-${root}.tsv
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
			} else {
				if {$reftranscripts ne "none"} {
					set samregions [samregions $region $refseq]
					convert_isoquant_reggenedb $reftranscripts $samregions $refseq regreftranscripts reggenedb
					convert_isoquant $isodir $destdir $sample $refseq $reggenedb $regreftranscripts $root $singlecell
				} else {
					convert_isoquant $isodir $destdir $sample $refseq "" "" $root
				}
			}
			if {$cleanup} {
				file delete -force $isodir
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
		set readfile $regdir/read_assignments-${root}.tsv
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
		isoform_counts-${root}.tsv
		gene_counts-${root}.tsv
		read_assignments-${root}.tsv
		totalcounts-${root}.tsv
	} -vars {
		isofiles genefiles readfiles sample refseq regdir region reftranscripts genedb root analysisname
		strictpct workdir
	} -code {
		analysisinfo_write [lindex $isofiles 0] isoform_counts-${root}.tsv
		analysisinfo_write [lindex $genefiles 0] gene_counts-${root}.tsv
		analysisinfo_write [lindex $readfiles 0] read_assignments-${root}.tsv
		iso_isoquant_mergeresults $isofiles $genefiles $readfiles $strictpct $sample $root $analysisname $singlecell
		cg tsv2gtf isoform_counts-${root}.tsv isoforms-${root}.gtf.gz
		catch {shadow_delete $workdir}
	}
}

proc cg_iso_isoquant {args} {
	set args [job_init {*}$args]
	iso_isoquant_job {*}$args
	job_wait
}
