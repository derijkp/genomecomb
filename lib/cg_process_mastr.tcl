
proc make_alternative_compar_job {experiment {destdir {}} {varcaller gatk}} {
	upvar job_logdir job_logdir
	if {$destdir eq ""} {set destdir [pwd]}
	job altcompar-$experiment -deps {
		$destdir/compar/annot_compar-$experiment.tsv
	} -vars {
		varcaller
	} -targets [list $destdir/compar/annot_compar_${varcaller}-${experiment}.tsv.zst $destdir/compar/annot_compar_${varcaller}-${experiment}_long.tsv.zst] \
	-code {
		set target1 [lindex $targets 0]
		set target2 [lindex $targets 1]
		##remove samtools analysis & move specific annotfields forwards		
		set cfields [cg select -h $dep]
		set cfields [list_sub $cfields -exclude [list_find -regexp $cfields -sam-]]
		set fields [list_common {chromosome begin end type ref alt amplicons dbnsfp_SIFT_score dbnsfp_Polyphen2_HDIV_score dbnsfp_Polyphen2_HDIV_pred dbnsfp_Polyphen2_HVAR_score dbnsfp_Polyphen2_HVAR_pred snp138_name 1000gCEU refGene_impact refGene_gene refGene_descr dbnsfp_MutationTaster_score dbnsfp_MutationTaster_pred} $cfields]
		lappend fields *
		lappend fields [string_change {log2_allele_ratio-gatk-crsbwa-*=if(llen(${alleledepth-gatk-crsbwa-*})>1, log10(lindex(${alleledepth-gatk-crsbwa-*},0))/log10(2) - log10(lindex(${alleledepth-gatk-crsbwa-*},1))/log10(2), 0)} [list gatk $varcaller]]
		set compress [compresspipe $target1]
		set temp [filetemp_ext $target1]
		exec cg select -rf {*-sam-*} $dep | cg select -f $fields {*}$compress > $temp
		file rename -force $temp $target1
		##depivot compar file
		set compress [compresspipe $target2]
		set temp [filetemp_ext $target2]
		if {$compress eq ""} {
			cg long $target1 $temp
		} else {
			exec cg long $target1 {*}$compress > $temp
		}
		file rename -force $temp $target2
		if {$compress ne ""} {cg_zindex $target}
	}
}

proc analysis_complete_job {experiment {destdir {}} {varcaller gatk}} {
	upvar job_logdir job_logdir
	if {$destdir eq ""} {set destdir [pwd]}
	job analysis_complete-$experiment -deps [list $destdir/coverage_${experiment}_avg.tsv $destdir/coverage_${experiment}_frac_above_20.tsv $destdir/compar/annot_compar_${varcaller}-${experiment}_long.tsv $destdir/${experiment}.html] \
	-targets {$destdir/analysis_complete} -vars destdir -code {
		file delete $destdir/analysis_running
		exec touch $target
	}
}

proc generate_coverage_report_job {experiment regfile histofiles {destdir {}}} {
	upvar job_logdir job_logdir
	if {$destdir eq ""} {set destdir [pwd]}
	job coverage_report-$experiment -deps [list $regfile {*}$histofiles] \
	-targets [list $destdir/coverage_${experiment}_avg.tsv $destdir/coverage_${experiment}_frac_above_20.tsv ] -code {
		set regfile [tempfile]
		cg regcollapse $dep1 > $regfile
		set oheader {name chr begin end}
		set names {}
		set regionlist [split [cg select -sh /dev/null -f {name chromosome begin end} $regfile] \n]
		foreach line $regionlist {
			set line [split $line \t]
			set avga($line) $line
			set fraca($line) $line
		}
		set histofiles [lrange $deps 1 end]
		foreach file $histofiles {
			if {![file exists $file]} continue
			set file [file_absolute $file]
			set temp [file split $file]
			if {[lindex $temp end-1] eq "reports"} {
				lappend oheader [lindex $temp end-2]
			} else {
				# old arrangment
				lappend oheader [lindex $temp end-1]
			}
			set f [open $file]
			set header [tsv_open $f]
			set poss [list_cor $header {name avg size {r<1} {r1<5} {r5<10} {r10<20}}]
			foreach region $regionlist {
				if {[gets $f line] == -1} {error "file $file too short"}
				set line [list_sub [split $line \t] $poss]
				foreach {name avg size r1 r5 r10 r20} $line break
				set regionname [lindex [split $region \t] 0]
				if {$name ne $regionname} {error "wrong name (order) in file $file"}
				set frac20 [format %.4g [expr {1-double($r1+$r5+$r10+$r20)/$size}]]
				lappend avga($region) $avg
				lappend fraca($region) $frac20
			}
			close $f
		}
		set o [open [lindex $targets 0] w]
		puts $o [join $oheader \t]
		foreach region $regionlist {
			puts $o $region\t[join $avga($region) \t]
		}
		close $o
		set o [open [lindex $targets 1] w]
		puts $o [join $oheader \t]
		foreach region $regionlist {
			puts $o $region\t[join $fraca($region) \t]
		}
		close $o
	}
}

# Needs R to be installed together with some R packages:
# install.packages("rmarkdown") ; install.packages("stringr") ; install.packages("dplyr") ; install.packages("tidyr") ; install.packages("googleVis") ; install.packages("DT")
proc generate_html_report_job {experiment {destdir {}}} {
	upvar job_logdir job_logdir
	if {$destdir eq ""} {set destdir [pwd]}
	job html_report-$experiment -deps {
		$destdir/compar/compar-${experiment}.tsv
		$destdir/coverage_${experiment}_avg.tsv
		$destdir/coverage_${experiment}_frac_above_20.tsv
		($destdir/reports/report_stats-$experiment.tsv)
		($destdir/report_stats-$experiment.tsv)
		($destdir/demultiplex_stats.tsv)
	} -targets {$destdir/$experiment.html} -vars {experiment destdir} -code {
		set keepdir [pwd]
		cd $destdir
		cg select -overwrite 1 -g analysis -gc {sequenced {v} count} $dep $destdir/compar/summary-compar-${experiment}.tsv
		set rmd $::appdir/res/mastrreport.Rmd
		set chartjs $::appdir/res/displayChartHistogram.js
		set cmd [string_change {library(rmarkdown); library(stringr); mastrdir="@destdir@"; local_jsapi="@chartjs@"; mastr <- str_replace(mastrdir,".*/([^/]*)","\\1"); render("@rmd@", output_file=paste(mastr,"html.temp",sep="."), output_dir = mastrdir)} [list @destdir@ $destdir @rmd@ $rmd @chartjs@ $chartjs]]
		exec [findR] -e $cmd >@ stdout 2>@ stderr
		file rename -force $target.temp $target
		file delete $destdir/compar/summary-compar-${experiment}.tsv
		cd $keepdir
	}
}

proc cg_process_mastr {args} {
	error "This command has been depricated, use cg process_project instead"
}

proc cg_process_mastrdesign {args} {
	error "This command has been depricated, use cg process_project instead"
}
