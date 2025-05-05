proc AnnotSV {args} {
	set refseq hg38
	set hpo {}
	cg_options AnnotSV args {
		-refseq {
			set value [file tail $value]
			if {$value in "hg38 GRCh38"} {
				set refseq GRCh38
			} elseif {$value in "hg37 GRCh37"} {
				set refseq GRCh37
			} elseif {$value in "mm10"} {
				set refseq mm10
			} elseif {$value in "mm9"} {
				set refseq mm9
			} else {
				error "usupported refseq ($value) for AnnotSV: must be one of hg38, hg37, mm10, mm9"
			}
			set refseq $value
		}
		-hpo {
		}
	} {svfile resultfile annotsvresultfile} 3 3
	set tempfile [tempfile].bed
	cg select -f {chromosome begin end {type=toupper($type)}} -sh /dev/null $svfile $tempfile
	set options {}
	if {$hpo ne ""} {
		lappend options -hpo $hpo
	}
	exec AnnotSV -SVinputFile $tempfile -outputFile $annotsvresultfile.temp -svtBEDcol 4 -genomeBuild $refseq {*}$options
	file rename -force $annotsvresultfile.temp.tsv $annotsvresultfile

	catch {gzclose $f} ; catch {gzclose $fa} ;  catch {gzclose $o} ;

	set f [gzopen $svfile]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 6 0]
	if {-1 in [lrange $poss 0 3]} {
		error "some necessary fields not found in $svfile: [list_sub {chromsome begin end type} [list_find [lrange $poss 0 3] -1]]"
	}
	set fa [open [list | cg select -s {SV_chrom SV_start SV_end SV_type} -q {$Annotation_mode eq "full"} $annotsvresultfile]]
	set aheader [tsv_open $fa]
	set addfields [list_common $aheader {
		Gene_name
		Closest_left
		Closest_right
		Gene_count
		RE_gene
		P_gain_phen
		P_gain_hpo
		P_gain_source
		P_gain_coord
		P_loss_phen
		P_loss_hpo
		P_loss_source
		P_loss_coord
		P_ins_phen
		P_ins_hpo
		P_ins_source
		P_ins_coord
		po_P_gain_phen
		po_P_gain_hpo
		po_P_gain_source
		po_P_gain_coord
		po_P_gain_percent
		po_P_loss_phen
		po_P_loss_hpo
		po_P_loss_source
		po_P_loss_coord
		po_P_loss_percent
		P_snvindel_nb
		P_snvindel_phen
		B_gain_source
		B_gain_coord
		B_gain_AFmax
		B_loss_source
		B_loss_coord
		B_loss_AFmax
		B_ins_source
		B_ins_coord
		B_ins_AFmax
		B_inv_source
		B_inv_coord
		B_inv_AFmax
		po_B_gain_allG_source
		po_B_gain_allG_coord
		po_B_gain_someG_source
		po_B_gain_someG_coord
		po_B_loss_allG_source
		po_B_loss_allG_coord
		po_B_loss_someG_source
		po_B_loss_someG_coord
		GC_content_left
		GC_content_right
		Repeat_coord_left
		Repeat_type_left
		Repeat_coord_right
		Repeat_type_right
		Gap_left
		Gap_right
		SegDup_left
		SegDup_right
		ENCODE_blacklist_left
		ENCODE_blacklist_characteristics_left
		ENCODE_blacklist_right
		ENCODE_blacklist_characteristics_right
		ACMG
		HI
		TS
		DDD_HI_percent
		ExAC_delZ
		ExAC_dupZ
		ExAC_cnvZ
		ExAC_synZ
		ExAC_misZ
		GenCC_disease
		GenCC_moi
		GenCC_classification
		GenCC_pmid
		NCBI_gene_ID
		OMIM_ID
		OMIM_phenotype
		OMIM_inheritance
		OMIM_morbid
		OMIM_morbid_candidate
		LOEUF_bin
		GnomAD_pLI
		ExAC_pLI
		PhenoGenius_score
		PhenoGenius_phenotype
		PhenoGenius_specificity
		Exomiser_gene_pheno_score
		Human_pheno_evidence
		Mouse_pheno_evidence
		Fish_pheno_evidence
		AnnotSV_ranking_score
		AnnotSV_ranking_criteria
		ACMG_class
	}]
	set aposs [list_cor $aheader {SV_chrom SV_start SV_end SV_type}]
	set newheader $header
	foreach field $addfields {
		regsub ^AnnotSV_ $field {} field
		lappend newheader AnnotSV_$field
	}
	set empty [list_fill [llength $addfields] {}]
	#
	set o [wgzopen $resultfile]
	puts $o $comment[join $newheader \t]
	while 1 {
		if {[gets $f line] == -1} break
		foreach {chr begin end type} [list_sub [split $line \t] $poss] break
		set add 1
		while 1 {
			if {[gets $fa aline] == -1} break
			set aline [split $aline \t]
			foreach {achr abegin aend atype} [list_sub $aline $aposs] break
			incr abegin -1
			if {$chr eq $achr && $begin == $abegin && $end == $aend && [string toupper $type] eq [string toupper $atype]} {
				set add 1
				break
			} else {
				error "annotation not found"
			}
		}
		if {$add} {
			lappend line {*}$aline
		} else {
			lappend line {*}$empty
		}
		puts $o [join $line \t]
	}


	gzclose $o
	gzclose $f
	gzclose $fa
}
