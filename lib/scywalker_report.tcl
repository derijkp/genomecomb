proc scywalker_report_job {sampledir} {
	upvar job_logdir job_logdir
	set sample [file tail $sampledir]
	set bam [lindex [jobglob $sampledir/map-sminimap2_splice-*.bam $sampledir/map-sminimap2*.bam] 0]
	set target $sampledir/reports/report_scywalker-$sample.html
	job scywalker_report_$sample -optional 1 -deps {
		$bam
		$sampledir/sc_cellinfo_raw-isoquant_sc-sminimap2_splice-$sample.tsv.zst
		$sampledir/read_assignments-isoquant_sc-sminimap2_splice-$sample.tsv.zst
		$sampledir/sc_celltype_umap-scsorter-isoquant_sc-sminimap2_splice-$sample.png
		$sampledir/sc_gene_counts_filtered-isoquant_sc-sminimap2_splice-$sample.tsv.zst
	} -targets {
		$target
	} -vars {
		sampledir sample bam
	} -code {
		exec cramino $bam > $sampledir/cramino-output.tsv
		exec scywalker-report -o $target.temp $sampledir
		file rename -force $target.temp $target
		mklink $sampledir/reports/report_scywalker-$sample.html $sampledir/report_scywalker-$sample.html
	}
}
