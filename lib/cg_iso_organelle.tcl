
proc iso_organelle_removeribo {overlaps {genes {MT-RNR1 MT-RNR2 MT-TF}}} {
	set ps {}
	set p 0
	set do 0
	foreach temp [list_subindex $overlaps 4] {
		if {$temp in $genes} {
			lappend ps $p
		} else {
			set do 1
		}
		incr p
	}
	if {$do} {
		set overlaps [list_sub $overlaps -exclude $ps]
	}
	return $overlaps
}

# a simplified caller for organelles (which give too many problems with normals callers because of different structure)
# It only counts known genes based on overlap. 
# It just picks the isoform/gene with the largest overlap, not taking into account that some organelles show some splicing/alt transcripts
proc iso_organelle_job {args} {
	upvar job_logdir job_logdir
	global appdir
	set cmdline [clean_cmdline cg iso_organelle {*}$args]
	set refseq {}
	set skips {}
	set reftranscripts {}
	set resultfile {}
	set region {}
	set options {}
	set singlecell 0
	set strictpct 90
	cg_options iso_organelle args {
		-refseq {
			set refseq $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-skip {
			lappend skips -skip $value
		}
		-region {
			set region $value
		}
	} {bam resultfile} 1 2
	if {$region eq ""} {
		error "region cannot be left empty for iso_organelle"
	}
	set bam [file_absolute $bam]
	set refseq [refseq $refseq]
	set reftranscripts [get_transcriptsfile_tsv $reftranscripts $refseq]
	if {$reftranscripts eq "none"} {
		error "reftranscripts none not allowed for iso_organelle"
	}
	if {$resultfile eq ""} {
		set root ${analysisname}-[file_rootname $bam]
		set resultfile [file dir $bam]/isoform_counts-$root.tsv
		set sample [file tail [file dir $bam]]
	} else {
		set resultfile [file_absolute $resultfile]
		set root [file_rootname $resultfile]
	}
	set analysisname iso_organelle
	set resultfile [file_absolute $resultfile]
	set regdir [file dir $resultfile]
	mkdir $regdir
	job iso_organelle-${root}-$region -mem 15G {*}$skips -deps {
		$bam $bam.bai $refseq $reftranscripts
	} -targets {
		$regdir/isoform_counts-${root}.tsv
		$regdir/gene_counts-${root}.tsv
		$regdir/read_assignments-${root}.tsv.zst
	} -vars {
		bam refseq regdir region reftranscripts threads root sample
		options analysisname strictpct
	} -code {
		if {[file exists $regdir.temp]} {file delete -force $regdir.temp}
		set info [list analysis $root sample $sample \
			isocaller_reference $refseq \
			isocaller_reftranscripts [file tail $reftranscripts] \
			isocaller $analysisname isocaller_version [version genomecomb]]
		analysisinfo_write $bam $regdir/isoform_counts-${root}.tsv {*}$info
		analysisinfo_write $bam $regdir/gene_counts-${root}.tsv {*}$info
		analysisinfo_write $bam $regdir/read_assignments-${root}.tsv {*}$info

		# read genes/isoforms
		set regreftranscripts [tempfile].tsv
		if {$region ne ""} {
			cg select -rc 1 -q "region(\"$region\")" $reftranscripts $regreftranscripts
		} else {
			cg select -rc 1 $reftranscripts $regreftranscripts
		}
		catch {gzclose $f}
		set f [gzopen $regreftranscripts]
		set header [tsv_open $f]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 11 12 13}]
		# list_sub $header $poss fields: chrom start end strand exonStarts exonEnds name gene geneid
		foreach {isopos genepos geneidpos} [lrange $poss end-2 end] break
		set tocheck {}
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {chrom start end strand exonStarts exonEnds transcript gene geneid} [list_sub $line $poss] break
			if {[regexp , [string trim $exonStarts ,]] || [regexp , [string trim $exonEnds ,]]} {
				puts stderr "warning: $transcript ($geneid $chrom:$start-$end $exonStarts $exonEnds) has splicing (which is not really supported in iso_organelle)"
			}
			lappend tocheck [list $start $end $strand $transcript $gene $geneid [expr {$end-$start}]]
		}
		gzclose $f

		set regtsvali [tempfile].tsv.zst
		set regtsvali $regdir/map-$region.sam.tsv.zst
		if {$region ne ""} {
			set samregions [samregions $region $refseq]
			exec samtools view -h $bam {*}$samregions | cg sam2tsv | cg zst -c 1 > $regtsvali
		} else {
			set samregions {}
			cg sam2tsv $bam | cg zst -c 1 > $regtsvali
		}
		# umicount
		catch {gzclose $f}
		set f [gzopen $regtsvali] 
		set header [tsv_open $f]
		set readpos [lsearch $header qname]
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			set read [lindex $line $readpos]
			if {[regexp {^([A-Z]+)_([A-Z]+)#(.*)$} $read temp cellbarcode umi r]} {
				if {![info exists umicount_readsa($r)]} {
					incr umicounta($cellbarcode,$umi)
					set umicount_readsa($r) 1
				}
			}
		}
		gzclose $f
		unset -nocomplain umicount_readsa

		# count
		catch {gzclode $ot} ; catch {gzclode $og} ; catch {gzclode $or}
		set ot [wgzopen $regdir/isoform_counts-${root}.tsv]
		puts $ot [join [list \
			chromosome begin end strand exonStarts exonEnds transcript gene geneid score bin \
			cdsStart cdsEnd exonCount name2 cdsStartStat cdsEndStat exonFrames ROW category size \
			counts_iqall-$root counts_weighed-$root counts_unique-$root counts_strict-$root counts_aweighed-$root \
		] \t]
		set og [wgzopen $regdir/gene_counts-${root}.tsv]
		puts $og [join {
			chromosome begin end strand gene geneid
		} \t]
		set or [wgzopen $regdir/read_assignments-${root}.tsv.zst]
		puts $or [join {
			read_id chromosome begin end strand exonStarts exonEnds aligned_size isoform_id 
			gene_id assignment_type assignment_events additional_info ambiguity inconsistency 
			covered_pct polya classification closest_known cellbarcode umi umicount
		} \t]
		catch {gzclose $f}
		set f [gzopen $regtsvali]
		set header [tsv_open $f]
		set poss [list_cor $header {chromosome begin end strand qname qstart qend mapquality cigar other seq}]
		unset -nocomplain tcounta
		unset -nocomplain gcounta
		set nr 0

		while 1 {
			incr nr ; if {![expr {$nr%10000}]} {puts $nr}
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {chromosome begin end strand read_id qstart qend mapquality cigar other seq} [list_sub $line $poss] break
			set exonStarts $begin ; set exonEnds $end ; set aligned_size [expr {$end-$begin}]
			set isoform_id . ; set gene_id . 
			set assignment_type noninformative ; set assignment_events . ; set additional_info .
			set ambiguity 0 ; set inconsistency 0 ;	set covered_pct 0 ; set polya {} 
			set classification intergenic ; set closest_known .
			set cellbarcode {} ; set umi {} ; set umicount 1
			if {[regexp {^([A-Z]+)_([A-Z]+)#(.*)$} $read_id temp cellbarcode umi r]} {
				set umicount $umicounta($cellbarcode,$umi)
			}
			if {[regexp N $cigar]} {
				set matches [cigar2exons $cigar $begin]
				set exonStarts [list_unmerge $matches 1 exonEnds]
				set exonStarts [join $exonStarts ,]
				set exonEnds [join $exonEnds ,]
				# putsvars exonStarts exonEnds
			} else {
				set matches [list $begin $end]
			}
			set overlaps {}
			set lmatches [llength $matches]
			set polya [iso_polya $seq $qstart $qend]
			foreach {begin end} $matches {
				foreach cline $tocheck {
					foreach {cbegin cend cstrand} $cline break
					if {$cend < $begin} continue
					if {$cbegin >= $end} break
					set covered_pct [format %.2f [expr {100.0*(min($cend,$end) - max($cbegin,$begin))/($cend - $cbegin)}]]
					# putsvars cbegin begin end cend
					if {[expr {abs($cbegin - $begin)}] <= 3 && [expr {abs($cend - $end)}] <= 3} {
						set assignment_events mono_exon_match
					} elseif {$begin >= [expr {$cbegin-3}] && $end <= [expr {$cend+3}]} {
						set assignment_events mono_exon_enclosed
					} else {
						if {$covered_pct < 5} continue
						set assignment_events mono_exon_overlap
					}
					set cpolya 0
					if {$polya == -1} {
						if {$cstrand == "-"} {set cpolya 1}
						set polya 0
					}
					lappend cline $assignment_events $covered_pct $cpolya
					lappend overlaps $cline
				}
			}
			set loverlaps [llength $overlaps]
			if {$loverlaps > 0} {
				if {$polya == 1 && $cstrand == "+"} {
					lset overlaps end 9 1
				}
				if {$loverlaps > 1} {
#					unset -nocomplain a
#					foreach overlap $overlaps {
#						set transcript [lindex $overlap 3]
#						lappend a($transcript) $overlap
#					}
#					set temp {}
#					foreach transcript [array names a] {
#						if {[llength $a($transcript)] > 1} {
#							lappend temp [lindex [lsort -index 9 -real $a($transcript)] end]
#						} else {
#							lappend temp [lindex $a($transcript) 0]
#						}
#					}
#					set overlaps $temp
					# remove ribosomal genes (human specific naming) from overlaps
					# if there are other matches
					set overlaps [iso_organelle_removeribo $overlaps {MT-RNR1 MT-RNR2 MT-TF}]
					# take the one with the best pct coverage (is in col 9)
					if {$loverlaps > 1} {
						set overlaps [list [lindex [lsort -index 9 -real $overlaps] end]]
					}
				}
				foreach overlap $overlaps {
					foreach {
						cbegin cend cstrand isoform_id gene gene_id size
						assignment_events covered_pct polya
					} $overlap break
					set closest_known $isoform_id
					set assignment_type unique
					if {$loverlaps > 2} {set assignment_events mono_exon_overlap}
					set classification mono_exon_match
					puts $or [join [list \
						$read_id $chromosome $begin $end $strand $exonStarts $exonEnds $aligned_size $isoform_id \
						$gene_id $assignment_type $assignment_events $additional_info $ambiguity $inconsistency \
						$covered_pct $polya $classification $closest_known $cellbarcode $umi $umicount \
					] \t]
					incr tcounta($isoform_id,t)
					if {$polya} {incr tcounta($isoform_id,t)}
					if {$assignment_events in "mono_exon_match mono_exon_enclosed"} {
						incr tcounta($isoform_id,u)
						if {$polya} {incr tcounta($isoform_id,u)}
						if {$covered_pct >= $strictpct} {
							incr tcounta($isoform_id,s)
							if {$polya} {incr tcounta($isoform_id,s)}
						}
					}
					incr gcounta($gene_id)
				}
			} else {
				# join [lrange $tocheck 0 5] \n
				puts $or [join [list \
					$read_id $chromosome $begin $end $strand $exonStarts $exonEnds $aligned_size $isoform_id \
					$gene_id $assignment_type $assignment_events $additional_info $ambiguity $inconsistency \
					$covered_pct $polya $classification $closest_known $cellbarcode $umi $umicount \
				] \t]
			}
		}

		gzclose $or
		# putsvars matches
		# putsvars overlaps

		#
		# make gene_counts and isoform_counts files
		# -----------------------------------------
		set o [iso_write_isoform_counts \
			$regdir/isoform_counts-${root}.tsv.temp \
			$regreftranscripts \
			tcounta \
			[list \
				counts_iqall-$root iq counts_weighed-$root t counts_unique-$root u counts_strict-$root s \
				counts_aweighed-$root a counts_aunique-$root au counts_astrict-$root as\
			]
		]
		gzclose $o
		file rename -force $regdir/isoform_counts-${root}.tsv.temp $regdir/isoform_counts-${root}.tsv

		#
		# make gene_counts file
		# ---------------------
		# make genebasica
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
		# write targetisoformcountsfile
		set targetgenecountsfile $regdir/gene_counts-${root}.tsv
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
}
