proc pb_combine_cat_colinfo {result files} {
	file copy [lindex $files 0] $result.temp
	set o [open $result.temp a]
	foreach file [lrange $files 1 end] {
		set f [open $file]
		set header [tsv_open $f]
		set pos [lsearch $header fieldtype]
		while {[gets $f line] != -1} {
			set type [lindex [split $line \t] $pos]
			if {$type ne "data"} continue
			puts $o $line
		}
		close $f
	}
	close $o
	file rename -force $result.temp $result
}

proc pb_combine_job {projectdir sc_celltyper {iso_match {}}} {
	upvar job_logdir job_logdir
	# combined analysis
	cd $projectdir
	mkdir compar
	set exproot [file tail $projectdir]
	if {$sc_celltyper eq "*"} {
		set root $exproot
	} else {
		set root ${sc_celltyper}-$exproot
	}
	set pbisoformfiles [bsort [jobglob -checkcompressed 1 samples/*/pb_isoform_counts-${sc_celltyper}-*.tsv]]
	set pbisoformcolinfofiles [jobglob -checkcompressed 1 samples/*/pb_isoform_counts-${sc_celltyper}-*.tsv.colinfo]
	if {[llength $pbisoformfiles]} {
		job pb_compar-pb_isoform_counts-$root \
		-deps $pbisoformfiles \
		-targets {
			compar/pb_isoform_counts-$root.tsv.zst
		} -vars {
			pbisoformfiles exproot root sc_celltyper iso_match
		} -code {
			set isoformcounts compar/pb_isoform_counts-$root.tsv
			cg multitranscript -match $iso_match $isoformcounts.temp.zst {*}$pbisoformfiles
			file rename -force $isoformcounts.temp.analysisinfo $isoformcounts.analysisinfo
			file rename -force $isoformcounts.temp.zst $isoformcounts.zst
		}
		job pb_compar-pb_isoform_counts_colinfo-$root \
		-deps $pbisoformcolinfofiles \
		-targets {
			compar/pb_isoform_counts-$root.tsv.colinfo
		} -vars {
			pbisoformcolinfofiles exproot root sc_celltyper iso_match
		} -code {
			set isoformcounts_colinfo compar/pb_isoform_counts-$root.tsv.colinfo
			pb_combine_cat_colinfo $isoformcounts_colinfo $pbisoformcolinfofiles
		}
	}
	set pbgenefiles [bsort [jobglob -checkcompressed 1 samples/*/pb_gene_counts-${sc_celltyper}-*.tsv]]
	set pbgenecolinfofiles [jobglob samples/*/pb_gene_counts-${sc_celltyper}-*.tsv.colinfo]
	if {[llength $pbgenefiles]} {
		job pb_compar-gene_counts-$root \
		-deps $pbgenefiles \
		-targets {
			compar/pb_gene_counts-$root.tsv.zst
		} -vars {
			pbgenefiles exproot root sc_celltyper iso_match
		} -code {
			set genecounts compar/pb_gene_counts-$root.tsv
			cg_multigene -match $iso_match $genecounts.temp.zst {*}$pbgenefiles
			file rename -force $genecounts.temp.analysisinfo $genecounts.analysisinfo
			file rename -force $genecounts.temp.zst $genecounts.zst
		}
		job pb_compar-pb_gene_counts_colinfo-$root \
		-deps $pbgenecolinfofiles \
		-targets {
			compar/pb_gene_counts-$root.tsv.colinfo
		} -vars {
			pbgenecolinfofiles exproot root sc_celltyper iso_match
		} -code {
			set genecounts_colinfo compar/pb_gene_counts-$root.tsv.colinfo
			pb_combine_cat_colinfo $genecounts_colinfo $pbgenecolinfofiles
		}
	}
}

