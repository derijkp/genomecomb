
proc sc_downstream {sampledir sc_filters sc_celltyping} {
	set sampledir [file_absolute $sampledir]
	set scgenefiles [jobgzfiles $sampledir/sc_gene_counts_raw-*.tsv]
	foreach scgenefile $scgenefiles {
		set scisoformfile [file dir $scgenefile]/[regsub ^sc_gene_counts_raw- [file tail $scgenefile] sc_isoform_counts_raw-]
		foreach sc_filter $sc_filters {
			cg_sc_filter_$sc_filter $scgenefile $scisoformfile
		}
	}
}
