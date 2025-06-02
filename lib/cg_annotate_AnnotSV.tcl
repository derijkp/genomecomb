proc AnnotSV_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg annotate_AnnotSV {*}$args]
	set refseq hg38
	set hpo {}
	set distrreg chr
	cg_options AnnotSV args {
		-refseq {
			set refseq $value
			set ref [file tail [refdir $value]]
			if {$ref in "hg38 GRCh38"} {
				set ref GRCh38
			} elseif {$ref in "hg19 hg37 GRCh37"} {
				set ref GRCh37
			} elseif {$ref in "mm39"} {
				set ref mm39
			} elseif {$ref in "mm10"} {
				set ref mm10
			} elseif {$ref in "mm9"} {
				set ref mm9
			} else {
				error "usupported ref ([file tail $refseq]) for AnnotSV: must be one of hg38, hg19, hg37, GRCh38, GRCh37"
			}
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-hpo {
			set value $value
		}
	} {svfile resultfile annotsvresultfile} 3 3 {
		add (extra) annotation to SV using AnnotSV
	}
	set svfile [file_absolute $svfile]
	set resultfile [file_absolute $resultfile]
	set annotsvresultfile [file_absolute $annotsvresultfile]
	# logfile
	job_logfile [file dir $resultfile]/annotate_AnnotSV_[file tail $resultfile] [file dir $resultfile] $cmdline \
		{*}[versions os AnnotSV]
	# logdir
	if {![info exists job_logdir]} {
		set_job_logdir [file dir $resultfile]/log_jobs
	}
	set resultname [file tail $resultfile]
	set workdir [shadow_workdir $resultfile]
	job_cleanup_add_shadow $workdir
	if {$distrreg == "0"} {
		job AnnotSV-[file tail $annotsvresultfile] -mem 16G -deps {
			$svfile
		} -targets {
			$annotsvresultfile
		} -vars {
			svfile annotsvresultfile hpo ref
		} -code {
			set tempfile [tempfile].bed
			# convert to the bed format as used by AnnotSV (ins and bnd end gets, incorrectly, +1)
			cg select -sh /dev/null \
				-f {chromosome begin {end=if($type in "ins bnd",$end+1,$end)} {type=toupper($type)}} \
				-q {$type in "del ins dup inv bnd" and not($chromosome regexp "_")} \
				$svfile $tempfile
			set options {}
			if {$hpo ne ""} {
				lappend options -hpo $hpo
			}
			exec AnnotSV -SVinputFile $tempfile -outputFile $annotsvresultfile.temp -svtBEDcol 4 -genomeBuild $ref {*}$options
			if {[gziscompressed $annotsvresultfile]} {
				# compress already takes care of temp and rename
				compress $annotsvresultfile.temp.tsv $annotsvresultfile
				file delete $annotsvresultfile.temp.tsv
			} else {
				file rename -force $annotsvresultfile.temp.tsv $annotsvresultfile
			}
		}
	} else {
		set chromosomes [distrreg_regs chr $refseq]
		set distrsrcs [distrreg_job -skip [list $annotsvresultfile] -refseq $refseq $svfile $workdir/svfile.part- .tsv $chromosomes]
		set todo {}
		foreach chromosome $chromosomes src $distrsrcs {
			if {[regexp _ $chromosome]} continue
			lappend todo $workdir/result.$chromosome
			job AnnotSV-$resultname-$chromosome -mem 4G -skip {
				$annotsvresultfile
			} -deps {
				$workdir/svfile.part-$chromosome.tsv
			} -targets {
				$workdir/result.$chromosome
			} -vars {
				workdir chromosome hpo ref
			} -code {
				set tempfile $workdir/svfile.part-$chromosome.bed
				cg select -overwrite 1 -sh /dev/null \
					-f {chromosome begin {end=if($type in "ins bnd",$end+1,$end)} {type=toupper($type)}} \
					-q {$type in "del ins dup inv bnd"} \
					$workdir/svfile.part-$chromosome.tsv $tempfile
				if {[file size $tempfile] != 0} {
					set options {}
					if {$hpo ne ""} {
						lappend options -hpo $hpo
					}
					# AnnotSV gave intermittent errors in testing larger data sets
					# to avoid/reduce these, do a retry
					if {[catch {
						exec AnnotSV -SVinputFile $tempfile -outputFile $workdir/result.$chromosome.temp -svtBEDcol 4 -genomeBuild $ref {*}$options
					} m]} {
						# really trying, a third retry
						if {[catch {
							exec AnnotSV -SVinputFile $tempfile -outputFile $workdir/result.$chromosome.temp -svtBEDcol 4 -genomeBuild $ref {*}$options
						} m]} {
							puts stderr "warning: retrying AnnotSV after error: $m"
							exec AnnotSV -SVinputFile $tempfile -outputFile $workdir/result.$chromosome.temp -svtBEDcol 4 -genomeBuild $ref {*}$options
						}
					}
				} else {
					file_write $workdir/result.$chromosome.temp.tsv {}
				}
				file rename -force $workdir/result.$chromosome.temp.tsv $workdir/result.$chromosome
				file delete $tempfile
			}
		}
		job AnnotSV-$resultname -deps $todo -targets {
			$annotsvresultfile
		} -vars {
			annotsvresultfile svfile todo
		} -code {
			analysisinfo_write $svfile $annotsvresultfile AnnotSV_version [version AnnotSV]
			set temp [filetemp $annotsvresultfile]
			cg cat -c 0 {*}$todo {*}[compresspipe $annotsvresultfile] > $temp
			file rename -- $temp $annotsvresultfile
		}
	}

	job annotate_AnnotSV-$resultname -deps {
		$svfile $annotsvresultfile
	} -targets {
		$resultfile
	} -vars {
		workdir svfile annotsvresultfile resultfile
	} -code {
		analysisinfo_write $svfile $resultfile AnnotSV_version [version AnnotSV]
		catch {gzclose $f} ; catch {gzclose $fa} ;  catch {gzclose $o} ;
		set f [gzopen $svfile]
		set header [tsv_open $f comment]
		set poss [tsv_basicfields $header 6 0]
		if {-1 in [lrange $poss 0 3]} {
			error "some necessary fields not found in $svfile: [list_sub {chromsome begin end type} [list_find [lrange $poss 0 3] -1]]"
		}
		cg select \
			-f {
				chromosome=$SV_chrom
				{begin=$SV_start - 1}
				{end=if($SV_type in "INS BND",$SV_end - 1,$SV_end)}
				{type="[string tolower $SV_type]"}
				*
			} -q {
				$Annotation_mode eq "full"
			} \
			$annotsvresultfile \
		| cg select -s - \
			-rf {SV_chrom SV_start SV_end SV_type} \
		 	> $workdir/temp.tsv
		set fa [open $workdir/temp.tsv]
		set aheader [tsv_open $fa]
		set addfields [list_common $aheader {
			Gene_name
			Closest_left Closest_right
			Gene_count RE_gene
			P_gain_phen P_gain_hpo P_gain_source P_gain_coord P_loss_phen
			P_loss_hpo P_loss_source P_loss_coord P_ins_phen P_ins_hpo P_ins_source P_ins_coord
			po_P_gain_phen po_P_gain_hpo po_P_gain_source po_P_gain_coord po_P_gain_percent
			po_P_loss_phen po_P_loss_hpo po_P_loss_source po_P_loss_coord po_P_loss_percent
			P_snvindel_nb P_snvindel_phen
			B_gain_source B_gain_coord B_gain_AFmax B_loss_source B_loss_coord B_loss_AFmax
			B_ins_source B_ins_coord B_ins_AFmax B_inv_source B_inv_coord B_inv_AFmax
			po_B_gain_allG_source po_B_gain_allG_coord po_B_gain_someG_source po_B_gain_someG_coord
			po_B_loss_allG_source po_B_loss_allG_coord po_B_loss_someG_source po_B_loss_someG_coord
			GC_content_left GC_content_right
			Repeat_coord_left Repeat_type_left Repeat_coord_right Repeat_type_right Gap_left Gap_right SegDup_left SegDup_right
			ENCODE_blacklist_left ENCODE_blacklist_characteristics_left ENCODE_blacklist_right ENCODE_blacklist_characteristics_right
			ACMG HI TS DDD_HI_percent
			ExAC_delZ ExAC_dupZ ExAC_cnvZ ExAC_synZ ExAC_misZ
			GenCC_disease GenCC_moi GenCC_classification GenCC_pmid
			NCBI_gene_ID OMIM_ID OMIM_phenotype OMIM_inheritance OMIM_morbid OMIM_morbid_candidate
			LOEUF_bin GnomAD_pLI ExAC_pLI
			PhenoGenius_score PhenoGenius_phenotype PhenoGenius_specificity Exomiser_gene_pheno_score
			Human_pheno_evidence Mouse_pheno_evidence Fish_pheno_evidence
			AnnotSV_ranking_score AnnotSV_ranking_criteria ACMG_class
		}]
		set addfields [list_lremove $aheader {
			chromosome begin end type
			AnnotSV_ID SV_chrom SV_start SV_end SV_length SV_type Samples_ID Annotation_mode
			CytoBand Tx Tx_version Tx_start Tx_end Overlapped_tx_length
			Overlapped_CDS_length Overlapped_CDS_percent Frameshift Exon_count
			Location Location2 Dist_nearest_SS Nearest_SS_type Intersect_start Intersect_end ROW
		}]
		set aposs [list_cor $aheader {chromosome begin end type}]
		set newheader $header
		foreach field $addfields {
			regsub ^AnnotSV_ $field {} field
			lappend newheader AnnotSV_$field
		}
		set empty [join [list_fill [llength $addfields] {}] \t]
		set afposs [list_cor $aheader $addfields]
		#
		set tempresultfile $workdir/finalresult.temp[gzext $resultfile]
		set o [wgzopen $tempresultfile]
		puts $o $comment[join $newheader \t]
		while 1 {
			if {[gets $f line] == -1} break
			# we use line later
			set sline [split $line \t]
# putsvars sline
			foreach {chr begin end type} [list_sub $sline $poss] break
			set add 0
			if {$type in "del ins dup inv bnd" && ![regexp _ $chr]} {
				while 1 {
					if {[gets $fa aline] == -1} break
					set aline [split $aline \t]
# putsvars aline
					foreach {achr abegin aend atype} [list_sub $aline $aposs] break
					if {$chr eq $achr && $begin == $abegin && $end == $aend && [string toupper $type] eq [string toupper $atype]} {
						set add 1
						break
					} elseif {$chr ne $achr && [bsort [list $achr $chr]] ne [list $achr $chr]} {
						# error if $achr > $chr
						error "annotation not found"
					} elseif {$chr eq $achr} {
						# if same chr, test begin and end
						if {$abegin > $begin || ($abegin == $begin && $aend > $end)} {
							error "annotation not found"
						}
					}
				}
			}
			if {$add} {
				append line \t[join [list_sub $aline $afposs] \t]
			} else {
				append line \t$empty
			}
			puts $o $line
		}
		gzclose $o
		gzclose $f
		gzclose $fa
		file rename -force $tempresultfile $resultfile
	}

}

proc cg_annotate_AnnotSV args {
	set args [job_init {*}$args]
	AnnotSV_job {*}$args
	job_wait
}
