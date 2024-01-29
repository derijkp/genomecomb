proc preset_scywalker {} {
	return {
		split 1
		paired 0
		clip 0
		removeduplicates 0
		aligners minimap2_splice
		realign 0
		varcallers {}
		svcallers {}
		methcallers {}
		counters {}
		isocallers isoquant_sc
		singlecell ontr10x
		sc_umisize 12
		sc_barcodesize 16
		sc_adaptorseq CTACACGACGCTCTTCCGATCT
		sc_filters default
		sc_celltypers {}
		distrreg g5000000
		iso_match novel
		reports {fastqstats singlecell flagstat_reads samstats histodepth hsmetrics vars covered histo}
	}
}

proc preset_ont {} {
	return {
		split 1
		paired 0
		clip 0
		removeduplicates 0
		aligners minimap2
		realign 0
		varcallers {clair3}
		svcallers {sniffles npinv cuteSV}
		methcallers {}
		counters {}
		isocallers {}
		singlecell {}
		distrreg chr
		reports {fastqstats flagstat_reads samstats histodepth hsmetrics vars covered histo}
	}
}

proc preset_ontr {} {
	return {
		split 1
		paired 0
		clip 0
		removeduplicates 0
		aligners minimap2_splice
		realign 0
		varcallers {clair3}
		svcallers {}
		methcallers {}
		counters {}
		isocallers {isoquant}
		iso_joint {isoquant}
		singlecell {}
		distrreg chr
		reports {-fastqc predictgender}
	}
}

proc preset_srs {} {
	# default removeduplicates to {} because then it set to 1 or 0 according to amplicons parameter
	return {
		split 1
		paired 1
		clip 1
		removeduplicates {}
		aligners bwa
		realign 0
		varcallers {gatkh strelka}
		svcallers {manta lumpy}
		methcallers {}
		counters {}
		isocallers {}
		iso_joint {}
		singlecell {}
		distrreg 30000000
		reports {basic}
	}
}

proc preset_rseq {} {
	# default removeduplicates to {} because then it set to 1 or 0 according to amplicons parameter
	return {
		split 1
		paired 1
		clip 1
		removeduplicates {}
		aligners star_2p
		realign 0
		varcallers {gatkh strelka}
		svcallers {}
		methcallers {}
		counters {qorts}
		isocallers {}
		iso_joint {}
		singlecell {}
		distrreg chr
		reports {fastqstats fastqc flagstat_reads histodepth hsmetrics vars covered histo}
	}
}

proc presets {} {
	set list [list_regsub -all ^preset_ [command_list preset_*] {}]
}