proc sc_demultiplex_read_dmfile {dmfile dmaVar} {
	upvar $dmaVar dma
	unset -nocomplain dma
	set f [gzopen $dmfile]
	set header [tsv_open $f]
	set cellpos [lindex [list_remove [list_cor $header {cell barcode callbarcode}] -1] 0]
	if {$cellpos eq ""} {
		error "no field \"cell\" or \"barcode\" found in $dmfile"
	}
	set dmsamplepos [lindex [list_remove [list_cor $header {sample donor_id subsample}] -1] 0]
	if {$dmsamplepos eq ""} {
		error "no field \"sample\" or \"donor_id\", or \"subsample\" found in $dmfile"
	}
	set poss [list $cellpos $dmsamplepos]
	while {[gets $f line] != -1} {
		foreach {cell dmsample} [list_sub [split $line \t] $poss] break
		regsub -- {-1$} $cell {} cell
		set dma($cell) $dmsample
	}
	gzclose $f
}

proc sc2bulk {scgenefile target} {
	set f [gzopen $scgenefile]
	set header [tsv_open $f]
	set countposs [list_find -regexp $header count]
	set cellpos [lindex [list_remove [list_cor $header {cell cellbarcode barcode}] -1] 0]
	set rmposs [list $cellpos {*}$countposs]
	set o [wgzopen $target]
	puts $o [join [list {*}[list_sub $header -exclude $rmposs] {*}[list_sub $header $countposs]] \t]
	set previd {}
	set ::tcl_precision 2
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set id [list_sub $line -exclude $rmposs]
		if {$id ne $previd} {
			if {$previd ne ""} {
#				set out [join $previd \t]
#				foreach v $counts {
#					append out \t[format %.2f $v]
#				}
#				puts $o $out
				puts $o [join $previd \t]\t[join $counts \t]
			}
			set counts [list_sub $line $countposs]
			set previd $id
		} else {
			set newcounts {}
			foreach pv $counts v [list_sub $line $countposs] {
				lappend newcounts [expr {$pv+$v}]
			}
			set counts $newcounts
		}
	}
	set ::tcl_precision 0
	gzclose $f
	gzclose $o
}

proc sc_demultiplex_job {sampledir dmfile refseq {destVar {}}} {
	upvar job_logdir job_logdir
	set sampledir [file_absolute $sampledir]
	set sample [file tail $sampledir]
	if {$destVar ne ""} {
		upvar $destVar dest
	} else {
		set f [gzopen $dmfile]
		set header [tsv_open $f]
		set field [lindex [list_common $header {sample donor_id subsample}] 0]
		gzclose $f
		set dmsamples [lrange [list_subindex [split [string trim [cg select -g $field $dmfile]] \n] 0] 1 end]
		unset -nocomplain dest
		foreach dmsample $dmsamples {
			set dest($dmsample) [file dir $sampledir]/dm_${dmsample}__$sample
			lappend dmsamplesdone $dest($dmsample)
			file mkdir $dest($dmsample)
		}
	}

	set dmsamples [array names dest]
	set files [jobgzfiles $sampledir/sc_gene_*.tsv $sampledir/sc_isoform_*.tsv \
		$sampledir/sc_cellinfo_*.tsv $sampledir/sc_group*.tsv \
		$sampledir/read_assignments-isoquant_sc-*.tsv $sampledir/map-*.bam $sampledir/map-*.cram]
	foreach file $files {
		set targets {}
		set tail [file tail $file]
		set gzext [gzext $tail]
		set gzroot [gzroot $tail]
		set fullext [file ext $gzroot]$gzext
		set temp [file root [gzroot $tail]]
		set pre [join [lrange [split $temp -] 0 end-1] -]
		foreach dmsample $dmsamples {
			set destsample [file tail $dest($dmsample)]
			set destfile $dest($dmsample)/$pre-$destsample$fullext
			lappend targets $destfile
		}
		job demultiplex-[file tail $file] -deps {
			$file $refseq
		} -targets $targets -vars {
			file dmfile dmsamples refseq
		} -code {
			sc_demultiplex_read_dmfile $dmfile dma
			set tail [file tail $file]
			set gzext [gzext $tail]
			set ext [file extension $tail]
			if {$ext in ".cram .bam"} {set bam 1} else {set bam 0}
			set comments {}
			if {$bam} {
				set f [open "| cg sam2tsv $file"]
				set header [tsv_open $f comments]
				set qnamepos [lsearch $header qname]
			} else {
				set f [gzopen $file]
				set header [tsv_open $f comments]
				set cellpos [lindex [list_remove [list_cor $header {cell barcode cellbarcode}] -1] 0]
				if {$cellpos eq ""} {
					error "no field \"cell\" or \"barcode\" found in $file"
				}
			}
			unset -nocomplain desta
			foreach dmsample $dmsamples destfile $targets {
				if {$bam} {
					set desta($dmsample) [open "| cg tsv2sam -refseq $refseq -outformat cram - $destfile.temp$ext" w]
				} else {
					set desta($dmsample) [wgzopen $destfile.temp$gzext]
				}
				puts $desta($dmsample) $comments[join $header \t]
			}
			if {$bam} {
				while {[gets $f line] != -1} {
					set qname [lindex [split $line \t] $qnamepos]
					if {![regexp {([^_]+)_[^#]+#} $qname temp cell]} continue
					if {![info exists dma($cell)]} continue
					puts $desta($dma($cell)) $line
				}
			} else {
				while {[gets $f line] != -1} {
					set cell [lindex [split $line \t] $cellpos]
					if {![info exists dma($cell)]} continue
					puts $desta($dma($cell)) $line
				}
			}

			foreach dmsample $dmsamples destfile $targets {
				gzclose $desta($dmsample)
				unset desta($dmsample)
				if {$bam} {
					file rename -force $destfile.temp$ext $destfile
					exec samtools index $destfile
				} else {
					file rename -force $destfile.temp$gzext $destfile
				}
			}
		}
	}
	# pre
	foreach dmsample $dmsamples {
		set sampledir $dest($dmsample)
		set isocallers {}
		foreach scgenefile [bsort [jobgzfiles $sampledir/sc_gene_counts_filtered-*.tsv]] {
			lappend isocallers [lindex [split [file tail $scgenefile] -] 1]
		}
		set isocallers [list_remdup $isocallers]
		set sc_celltypers {}
		foreach groupfile [bsort [jobgzfiles $sampledir/sc_group-*.tsv]] {
			lappend sc_celltypers [lindex [split [file tail $groupfile] -] 1]
		}
		set sc_celltypers [list_remdup $sc_celltypers]
		foreach isocaller $isocallers {
			set scgenefile [jobgzfile $sampledir/sc_gene_counts_raw-$isocaller-*.tsv $sampledir/sc_gene_counts_filtered-$isocaller-*.tsv]
			set scisoformfile [jobgzfile $sampledir/sc_isoform_counts_raw-$isocaller-*.tsv $sampledir/sc_isoform_counts_filtered-$isocaller-*.tsv]
			set tail [file tail $scgenefile]
			if {![regsub {^sc_gene_counts_[^-]+-} $tail gene_counts- tail]} {
				regsub {^sc_} $tail {} tail
			}
			set target $sampledir/$tail
			sc2bulk $scgenefile $target
			set tail [file tail $scisoformfile]
			if {![regsub {^sc_isoform_counts_[^-]+-} $tail isoform_counts- tail]} {
				regsub {^sc_} $tail {} tail
			}
			set target $sampledir/$tail
			sc2bulk $scisoformfile $target

			set scgenefile [jobgzfile $sampledir/sc_gene_counts_filtered-$isocaller-*.tsv]
			set scisoformfile [jobgzfile $sampledir/sc_isoform_counts_filtered-$isocaller-*.tsv]
			foreach sc_celltyper $sc_celltypers {
				set groupfile [gzfile $sampledir/sc_group-$sc_celltyper-$isocaller-*.tsv]
				sc_pseudobulk_job $scgenefile $scisoformfile $groupfile
			}
		}
	}
}

proc cg_sc_demultiplex {args} {
	set args [job_init {*}$args]
	sc_demultiplex_job {*}$args
	job_wait
}
