proc sc_demultiplex_read_dmfile {dmfile dmaVar} {
	upvar $dmaVar dma
	unset -nocomplain dma
	set f [gzopen $dmfile]
	set header [tsv_open $f]
	set cellpos [lindex [list_remove [list_cor $header {cell barcode callbarcode}] -1] 0]
	if {$cellpos eq ""} {
		error "no field \"cell\" or \"barcode\" found in $dmfile"
	}
	set dmsamplepos [lindex [list_remove [list_cor $header {sample donor_id}] -1] 0]
	if {$dmsamplepos eq ""} {
		error "no field \"sample\" or \"donor_id\" found in $dmfile"
	}
	set poss [list $cellpos $dmsamplepos]
	while {[gets $f line] != -1} {
		foreach {cell dmsample} [list_sub [split $line \t] $poss] break
		regsub -- {-1$} $cell {} cell
		set dma($cell) $dmsample
	}
	gzclose $f
}

proc sc_demultiplex_job {sampledir dmfile refseq destVar} {
	upvar job_logdir job_logdir
	set sampledir [file_absolute $sampledir]
	upvar $destVar dest
	set dmsamples [array names dest]
	set files [gzfiles $sampledir/sc_gene_*.tsv $sampledir/sc_isoform_*.tsv \
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
				} else {
					file rename -force $destfile.temp$gzext $destfile
				}
			}
		}

	}
}

proc cg_sc_demultiplex {args} {
	set args [job_init {*}$args]
	sc_demultiplex_job {*}$args
	job_wait
}
