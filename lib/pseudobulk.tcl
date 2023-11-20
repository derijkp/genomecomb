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
	set pbisoformfiles [bsort [jobglob samples/*/pb_isoform_counts-${sc_celltyper}-*.tsv]]	
	set pbisoformfiles [list_lremove $pbisoformfiles [jobglob samples/*/pb_isoform_counts-${sc_celltyper}-*.colinfo.tsv]]	
	if {[llength $pbisoformfiles]} {
		job pb_compar-pb_isoform_counts-$root \
		-deps $pbisoformfiles \
		-targets {
			compar/pb_isoform_counts-$root.tsv
		} -vars {
			pbisoformfiles exproot root sc_celltyper iso_match
		} -code {
			set isoformcounts compar/pb_isoform_counts-$root.tsv
			cg multitranscript -match $iso_match $isoformcounts.temp.zst {*}$pbisoformfiles
			file rename -force $isoformcounts.temp.zst $isoformcounts.zst
		}
	}
	set pbgenefiles [bsort [jobglob samples/*/pb_gene_counts-${sc_celltyper}-*.tsv]]
	set pbgenefiles [list_lremove $pbgenefiles [jobglob samples/*/pb_gene_counts-${sc_celltyper}-*.colinfo.tsv]]	
	if {[llength $pbgenefiles]} {
		job pb_compar-gene_counts-$root \
		-deps $pbgenefiles \
		-targets {
			compar/gene_counts-$root.tsv
		} -vars {
			pbgenefiles exproot root sc_celltyper iso_match
		} -code {
			set genecounts compar/pb_gene_counts-$root.tsv
			cg_multigene -match $iso_match $genecounts.temp.zst {*}$pbgenefiles
			file rename -force $genecounts.temp.zst $genecounts.zst
		}
	}
}

