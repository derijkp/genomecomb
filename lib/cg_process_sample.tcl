#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc process_compare_checkfield {file field} {
	set f [open $file]
	set header [tsv_open $f]
	close $f
	inlist $header $field
}

proc process_indexcompress {file} {
	set ext [file extension $file]
	if {$ext eq ".gz"} {
		gunzip $file
		set file [file root $file]
	}
	set f [gzopen $file]
	set header [tsv_open $f]
	foreach field {offset end1} {
		set fpos [lsearch $header $field]
		if {$fpos != -1} break
	}
	if {$fpos == -1} {error "no column offset, end1 in file $file"}
	if {([gets $f] eq "") && [eof $f]} return
	gzclose $f
	tsv_index $field $file
	putslog "Compressing $file"
	if {![inlist {.rz .bgz} $ext]} {
		exec bgzip -c $file > $file.gz.temp
		file rename -force -- $file.gz.temp $file.gz
		file delete $file
	}
}

proc cg_process_indexcompress {args} {
	global scriptname action
	if {[llength $args] != 1} {
		error "format is: $scriptname $action file\n - makes index, and compresses to bgzip"
	}
	foreach {file} $args break
	process_indexcompress $file
}

proc process_sample_cgi_job {workdir split} {
	set sample [file tail $workdir]
	set keepdir [pwd]
	cd $workdir
	set_job_logdir $workdir/log_jobs
	set chromosomes {}
	set files [jobglob -checkcompressed 1 $workdir/ori/ASM/REF/coverage*.tsv]
	if {[llength $files]} {
		foreach file $files {
			lappend chromosomes [chr_clip [lindex [split [file tail $file] -] 1]]
		}
		set chromosomes [bsort [list_remdup $chromosomes]]
	} else {
		set files [jobglob -checkcompressed 1 $workdir/coverage/coverage*.tsv]
		foreach file $files {
			lappend chromosomes [chr_clip [lindex [split [file tail $file] -] end-1]]
		}
		set chromosomes [bsort [list_remdup $chromosomes]]
	}
	# start from CGI data
	# convert overage files to bcol first (will be used to add coverage and refscore to vars)
	set files [bsort [glob -nocomplain ori/ASM/REF/coverage*-chr*-*]]
	if {[llength $files]} {
		# this will only work if ori/ASM/REF/coverage*-chr* already exist from the start
		# maybe later make more flexible
		set chrs {}
		set header {}
		foreach file $files {
			if {[catch {regexp {^ori/ASM/REF/coverage.*-chr([^-]*)-.*$} $file temp chr}]} {error "Could not find chromosome in filename $file"}
			set chr [chr_clip $chr]
			lappend chrs $chr
			set f [gzopen  $file]
			set nheader [tsv_open $f]
			gzclose $f
			if {![llength $header]} {
				set header $nheader
			} elseif {$nheader ne $header} {
				error "file $file has a different header"
			}
		}
		if {![llength [list_common $header {uniqueSequenceCoverage coverage}]]} {
			error "No coverage/uniqueSequenceCoverage field found in $file"
		}
		set fields [list_remove $header offset]
		set poss [list_cor $header $fields]
		foreach posfield {offset pos} {
			set offsetpos [lsearch $header $posfield]
			if {$offsetpos != -1}  break
		}
		if {$offsetpos == -1} {
			error "No position/offset field found in $file"
		}
		foreach field $fields {
			if {$field eq "uniqueSequenceCoverage"} {set outfield coverage} else {set outfield $field}
			file mkdir bcolall
			set finaltarget bcolall/$outfield-cg-cg-$sample.bcol
			set tomerge {}
			set tomergebins {}
			foreach file $files chr $chrs {
				set target bcolall/$outfield-cg-cg-$chr-$sample.bcol
				lappend tomerge $target
				lappend tomergebins $target.bin
				job cg_coverage-cg-$sample-$outfield-$chr-$sample -checkcompressed 1 -deps $file \
				   -vars {sample chr field posfield} \
				   -skip [list $finaltarget $finaltarget.bin.zst $finaltarget.bin.zst.zsti] \
				   -targets {$target $target.bin} -code {
					# make coverage files
					set file $dep
					exec {*}[gzcat $file] $file | cg bcol make -n $chr -co 0 -p $posfield -t i $target $field >@ stdout 2>@ stderr
				}
			}
			job cg_coverage-cg-${sample}_merge-$outfield-$sample \
			    -checkcompressed 1 \
			    -vars {tomerge tomergebins finaltarget} \
			    -deps [list {*}$tomerge {*}$tomergebins] \
			    -rmtargets [list {*}$tomerge {*}$tomergebins] \
			    -targets {$finaltarget $finaltarget.bin.zst $finaltarget.bin.zst.zsti} -code {
				cg cat -c f {*}$tomerge > $finaltarget.temp
				exec cat {*}$tomergebins > $finaltarget.bin.temp
				zst -keep 0 -i 1 -o $finaltarget.bin.zst $finaltarget.bin.temp
				file delete $finaltarget.bin.temp
				file rename $finaltarget.temp $finaltarget
				file delete {*}$deps
			}
		}
	}
	# variants
	job cg_svar-$sample -checkcompressed 1 -optional 1 -deps {ori/ASM/var-*-ASM*.tsv} -targets {svar-$sample.tsv} \
	-skip [list var-cg-cg-$sample.tsv sreg-cg-cg-$sample.tsv reg_refcons-$sample.tsv reg_nocall-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv] \
	-code {
		set varfile $dep
		set f [gzopen $varfile]
		set info {}
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] ne "#"} break
			lappend info $line
		}
		catch {gzclose $f}
		if {[file exists info.txt]} {
			set test [split [file_read info.txt] \n]
			if {$info ne $test} {
				error "$target1 already has info.txt that differs from data in the source $dir"
			}
		}
		file_write info.txt [join $info \n]
		putslog "Sort var file ($varfile)"
		cg select -s "chromosome begin end varType" $varfile $target.temp
		file rename -force -- $target.temp $target
	}
	# we no longer convert the CG gene info. better annotated by genomecomb itself later
#	job cg_sgene-$sample -checkcompressed 1 -optional 1 -deps {ori/ASM/gene-*-ASM*.tsv} -targets {sgene-$sample.tsv} \
#	-skip [list var-cg-cg-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv] -code {
#		set genefile $dep
#		if {[llength $genefile] != 1} {error "could not identify genefile"}
#		putslog "Sort gene file ($genefile)"
#		cg select -s "chromosome begin end" $genefile $target.temp
#		file rename -force -- $target.temp $target
#	}
	# annotated vars file
	job cg_annotvar-$sample -checkcompressed 1 -optional 1 -vars {split sample} \
	-deps {svar-$sample.tsv (sgene-$sample.tsv) (bcolall/coverage-cg-cg-$sample.bcol) (bcolall/refScore-cg-cg-$sample.bcol)} \
	-targets {annotvar-$sample.tsv} \
	-skip [list var-cg-cg-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv] -code {
		putslog "Create annotated varfile $target"
		if {[file exists $dep2]} {
			cg cg2tsv -split $split -sorted 1 $dep1 $dep2 $target.temp
		} else {
			cg cg2tsv -split $split -sorted 1 $dep1 $target.temp
		}
		set todo {}
		if {[file exists bcolall/coverage-cg-cg-$sample.bcol]} {
			lappend todo bcolall/coverage-cg-cg-$sample.bcol
		}
		if {[file exists bcolall/refScore-cg-cg-$sample.bcol]} {
			lappend todo bcolall/refScore-cg-cg-$sample.bcol
		}
		if {[llength $todo]} {
			cg annotate $target.temp $target.temp2 {*}$todo
			file rename -force -- $target.temp2 $target
			file delete -force $target.temp $target.temp2.analysisinfo $target.temp.index $target.temp2.index
		} else {
			file rename -force -- $target.temp $target
		}
	}
#	# if not from CGI, we do not have svar and sgene, take from first var_* file found
#	job cg_annotvar_var-$sample -checkcompressed 1 -optional 1 {vars_*.tsv} {annotvar-$sample.tsv} {
#		gzmklink $dep $target
#	}
#	# if we also do not find a var_* file, take from first variant* file found
#	job cg_annotvar_variant-$sample -checkcompressed 1 -optional 1 {variant*.tsv} {annotvar-$sample.tsv} {
#		set file [gzfile $dep]
#		gzmklink $file $target
#	}
	# only if reg file exists, otherwise extract from svar (next)
	job cg_sreg-cg-cg-$sample -checkcompressed 1 -optional 1 -deps {ori/ASM/reg-*-ASM*.tsv} \
	-targets {sreg-cg-cg-$sample.tsv.zst} -code {
		set regfile [gzfile $dep]
		putslog "Sort region file ($regfile)"
		cg select -s "chromosome begin end" $regfile $target.temp.zst
		file rename -force -- $target.temp.zst $target
	}
	job cg_regfromsvar-$sample -checkcompressed 1 -optional 1 -deps {svar-$sample.tsv} \
	-targets {sreg-cg-cg-$sample.tsv.zst} -code {
		set svarfile $dep
		putslog "Extract $target from $svarfile"
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $svarfile $target.temp
		exec cg regjoin $target.temp {*}[compresspipe .zst 9] > $target.temp2
		file rename -force -- $target.temp2 $target
		file delete $target.temp
	}
	job cg_reg_refcons-$sample -checkcompressed 1 -optional 1 -deps {svar-$sample.tsv} \
	-targets {reg_refcons-$sample.tsv} -code {
		putslog "Find refcons regions for $dep"
		cg refconsregions $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job cg_reg_nocall-$sample -checkcompressed 1 -optional 1 -deps {svar-$sample.tsv} \
	-targets {reg_nocall-$sample.tsv} -code {
		putslog "Find partial no-call regions for dep"
		if {[catch {
			nocallregions $dep $target.temp
		}]} {
			puts stderr "Could not make $job (old version files ?)"
		} else {
			file rename -force -- $target.temp $target
		}
	}
	job cg_cpSV-$sample -checkcompressed 1 -optional 1 -deps {ori/ASM/SV ^ori/ASM/SV/(.*)$} -targets {SV/\_} -code {
		putslog "Copying SV"
		set targetdir [file dir $target]
		set tempdir [filetemp $targetdir]
		file delete $tempdir
		file copy $dep $tempdir
		file delete -force $targetdir
		file rename -force -- $tempdir $targetdir
	}
	job cg_cgsv-$sample -checkcompressed 1 -optional 1 -deps {SV/allJunctionsBeta-*.tsv} -targets {cgsv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		file delete $tempfile
		convcgsv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cgsv_alpha-$sample -checkcompressed 1 -optional 1 -deps {SV/annotatedJunctionsAlpha-*.tsv} -targets {cgsv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		file delete $tempfile
		convcgsv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cpCNV-$sample -checkcompressed 1 -optional 1 -deps {ori/ASM/CNV ^ori/ASM/CNV/(.*)$} -targets {CNV/\_} -code {
		putslog "Copying CNV"
		set targetdir [file dir $target]
		set tempdir [filetemp $targetdir]
		file delete $tempdir
		file copy $dep $tempdir
		file delete -force $targetdir
		file rename -force -- $tempdir $targetdir
	}
	job cg_cgcnv -checkcompressed 1 -optional 1 -deps {CNV/cnvSegmentsBeta-*.tsv} -targets {cgcnv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		convcgcnv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cgcnv_diploid-$sample -checkcompressed 1 -optional 1 -deps {CNV/cnvSegmentsDiploidBeta-*.tsv} -targets {cgcnv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		convcgcnv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cgcnv_alpha-$sample -checkcompressed 1 -optional 1 {CNV/cnvSegmentsAlpha-*.tsv} {cgcnv-$sample.tsv} {
		set tempfile [filetemp $target]
		convcgcnv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	# multiarch
	job reg_cluster-$sample -checkcompressed 1 -optional 1 -deps {
		annotvar-$sample.tsv
	} -targets {
		reg_cluster-$sample.tsv.zst
	} -code {
		cg clusterregions < $dep > $target.temp
		zst $target.temp
		file rename -force -- $target.temp.zst $target
	}
	job reg_ns-$sample -checkcompressed 1 -optional 1 -deps {
		annotvar-$sample.tsv
	} -targets {
		reg_ns-$sample.tsv
	} -code {
		putslog "Find regions with N's for $dep"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job reg_lowscore-$sample -checkcompressed 1 -optional 1 -deps {
		annotvar-$sample.tsv
	} -targets {
		reg_lowscore-$sample.tsv
	} -code {
		set header [cg select -h $dep]
		if {[llength [list_common $header {totalScore1 totalScore2}]] == 2} {
			putslog "Find regions with lowscores for $dep"
			cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < $dep > $target.temp
			file rename -force -- $target.temp $target
		}
	}
	job cg_var-cg-cg-$sample -checkcompressed 1 -optional 1 -deps {
		annotvar-$sample.tsv (reg_refcons-$sample.tsv) (reg_cluster-$sample.tsv) (coverage-cg-$sample/bcol_coverage-$sample.tsv) (coverage-cg-$sample/bcol_refscore-$sample.tsv)
	} -targets {
		var-cg-cg-$sample.tsv.zst var-cg-cg-$sample.tsv.analysisinfo
	} -vars {
		sample
	} -code {
		set cgi_version ?
		set cgi_reference ?
		if {[file exists info.txt]} {
			set c [file_read info.txt]
			regexp {SOFTWARE_VERSION\t([^\n]+)} $c temp cgi_version
			regexp {GENOME_REFERENCE\t([^\n]+)} $c temp cgi_reference
			if {$cgi_reference eq "NCBI build 37"} {
				set reference hg19
			} elseif {$cgi_reference eq "NCBI build 36"} {
				set reference hg18
			} else {
				set reference $cgi_reference
			}
		}
		set tempfile [filetemp_ext $target 0]
		analysisinfo_write {} $dep sample cg-cg-$sample aligner cgi aligner_version $cgi_version varcaller cgi varcaller_version $cgi_version reference $reference
		cg annotate $dep $tempfile {*}[list_remove [lrange $deps 1 end] {}]
		file rename $tempfile $target
		file rename $tempfile.zsti $target.zsti
		file rename [analysisinfo_file $tempfile] [analysisinfo_file $target]
		file delete -force [gzroot $tempfile].index
		file delete -force [gzroot $dep].index
	}
	job reg_covered-$sample -checkcompressed 1 -optional 1 -deps {sreg-cg-cg-$sample.tsv} -targets {reg-$sample.covered} -code {
		putslog "Genomic coverage of sequenced regions"
		cg covered $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	file mkdir filtered
	file mkdir covered
	job cg_filteredrefcons-$sample -checkcompressed 1 -optional 1 -vars sample \
	-deps {sreg-cg-cg-$sample.tsv reg_refcons-$sample.tsv} \
	-targets {filtered/filteredrefcons-$sample.tsv covered/filteredrefcons-$sample.covered} \
	-code {
		putslog "Coverage of refcons region"
		set temp1 [filetemp $target1]
		set temp2 [filetemp $target2]
		cg regsubtract $dep1 $dep2 > $temp1
		cg covered $temp1 > $temp2
		zst -keep 0 -i 1 -o $target1.zst $temp1
		file rename -force -- $temp2 $target2
	}
	job cg_filteredns-$sample -checkcompressed 1 -optional 1 -deps {sreg-cg-cg-$sample.tsv reg_ns-$sample.tsv} \
	-targets {filtered/filteredns-$sample.tsv.zst} -code {
		putslog "Coverage of ns region"
		cg regsubtract $dep1 $dep2 > $target.temp
		zst -keep 0 -i 1 -o $target $target.temp
	}
	job cg_filteredns_covered-$sample -checkcompressed 1 -optional 1 -deps {filtered/filteredns-$sample.tsv} \
	-targets {covered/filteredns-$sample.covered} -code {
		putslog "Making $target"
		set temp [filetemp $target]
		cg covered $dep > $temp
		file rename -force -- $temp $target
	}
	job cg_filteredlowscore-$sample -checkcompressed 1 -optional 1 -deps {sreg-cg-cg-$sample.tsv reg_lowscore-$sample.tsv} \
	-targets {filtered/filteredlowscore-$sample.tsv.zst} -code {
		set temp [filetemp $target]
		cg regsubtract $dep1 $dep2 > $temp
		zst -keep 0 -i 1 -o $target $temp
	}
	job cg_filteredlowscore_covered-$sample -checkcompressed 1 -optional 1 -deps {filtered/filteredlowscore-$sample.tsv} \
	-targets {covered/filteredlowscore-$sample.covered} -code {
		set temp [filetemp $target]
		cg covered $dep > $temp
		file rename -force -- $temp $target
	}
	job cg_refcons_histo-$sample -checkcompressed 1 -optional 1 -deps {reg_refcons-$sample.tsv} -targets {histo-refcons-$sample.tsv} -code {
		putslog "Making $target"
		cg reghisto $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job cg_filteredcluster-$sample -checkcompressed 1 -optional 1 -deps {sreg-cg-cg-$sample.tsv reg_cluster-$sample.tsv} \
	-targets {filtered/filteredcluster-$sample.tsv.zst} -code {
		putslog "Coverage of clusters region"
		set temp [filetemp $target]
		cg regsubtract $dep1 $dep2 > $temp
		zst -keep 0 -i 1 -o $target $temp
	}
	job cg_filteredcluster_covered-$sample -checkcompressed 1 -optional 1 -deps {filtered/filteredcluster-$sample.tsv} \
	-targets {covered/filteredcluster-$sample.covered} -code {
		putslog "Making $target"
		cg covered $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job cg_process_cleanup-$sample -checkcompressed 1 -optional 1 \
		-deps {
			(svar-$sample.tsv) (annotvar-$sample.tsv) (annotvar-$sample.tsv.index) (sgene-$sample.tsv)
			var-cg-cg-$sample.tsv sreg-cg-cg-$sample.tsv
			reg_refcons-$sample.tsv reg_nocall-$sample.tsv
			reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv
		} -vars {
			sample
		} -rmtargets {
			svar-$sample.tsv annotvar-$sample.tsv sgene-$sample.tsv
		} -code {
			catch {file delete svar-$sample.tsv}
			catch {file delete annotvar-$sample.tsv}
			catch {file delete annotvar-$sample.tsv.index}
			catch {file delete sgene-$sample.tsv}
	}
	job cg_process_summary-$sample -checkcompressed 1 -deps {
		var-cg-cg-$sample.tsv
		sreg-cg-cg-$sample.tsv
		reg-$sample.covered
		bcolall/coverage-cg-cg-$sample.bcol
		(reg_refcons-$sample.tsv)
		(reg_nocall-$sample.tsv)
		(SV) (CNV)
		(cgsv-$sample.tsv)
		(cgcnv-$sample.tsv)
		(reg_cluster-$sample.tsv)
		(reg_ns-$sample.tsv)
		(reg_lowscore-$sample.tsv)
		(reg-$sample.covered)
		(covered/filteredrefcons-$sample.covered)
		(covered/filteredns-$sample.covered)
		(covered/filteredlowscore-$sample.covered)
		(covered/filteredcluster-$sample.covered)
		(histo-refcons-$sample.tsv)
	} -targets {summary-$sample.txt} -vars {sample} -code {
		set f [open $target.temp w]
		puts $f "finished\t\t[job_timestamp]"
		set c [split [exec cg select -f chromosome $dep | uniq -c] \n]
		puts $f ""
		puts $f "Variants:\nchromosome\tvalue"
		set total 0
		foreach line [lrange $c 1 end] {
			foreach {num chr} $line break
			puts $f "$chr\t$num"
			incr total $num
		}
		puts $f "total\t$total"
		puts $f ""
		puts $f "covered:\nchromosome\tvalue"
		set c [split [file_read $dep3] \n]
		set total 0
		foreach line [lrange $c 1 end] {
			if {![llength $line]} break
			foreach {chr num} $line break
			puts $f "$chr\t$num"
			incr total $num
		}
		puts $f "total\t$total"
		set header [cg select -h $dep]
		if {[inlist $header impact]} {
			puts $f ""
			puts $f "impact:\nimpact\tnumber_of_variants"
			set c [split [exec cg select -s impact -f impact $dep | uniq -c] \n]
			foreach line [lrange $c 1 end] {
				foreach {num impact} $line break
				puts $f "$impact\t$num"
			}
		}
		close $f
		file rename -force -- $target.temp $target
	}
	cd $keepdir
}

proc checkminreads {fastqdir minfastqreads numVar} {
	upvar $numVar num
	set files [bsort [jobglob \
		$fastqdir/*.fastq.gz $fastqdir/*.fastq $fastqdir/*.fq.gz $fastqdir/*.fq \
		$fastqdir/*.bam $fastqdir/*.cram $fastqdir/*.sam \
	]]
	if {![llength $files]} {
		set num 0
		return 0
	}
	set file [lindex $files 0]
	set f [open "| [convert_pipe $file .tsv]"]
	set header [tsv_open $f]
	set count $minfastqreads
	while {$count} {
		if {[gets $f line] == -1} {
			gzclose $f
			set num [expr {$minfastqreads-$count}]
			return 0
		}
		incr count -1
	}
	gzclose $f
	return 1
}

proc process_sample_reports_minfastqreads {sampledir sample reports num todoVar} {
	upvar $todoVar todo
	if {[inlist $reports fastqstats] || [inlist $reports all] || [inlist $reports basic]} {
		file mkdir $sampledir/reports
		file_write $sampledir/reports/report_fastq_fw-$sample.tsv [join [list \
			[join {sample source parameter value} \t] \
			[join [list $sample fastq-stats fw_numreads $num] \t] \
		] \n]\n
		file_write $sampledir/reports/report_fastq_rev-$sample.tsv [join [list \
			[join {sample source parameter value} \t] \
			[join [list $sample fastq-stats rev_numreads $num] \t] \
		] \n]\n
		lappend todo(reports) $sampledir/reports
	}
}

proc process_sample_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg process_sample {*}$args]
	set keepargs $args
	set dbdir {}
	set minfastqreads 1
	set clip 1
	set aligners bwa
	set ali_keepcomments {}
	set varcallers {gatk sam}
	set svcallers {}
	set methcallers {}
	set counters {}
	set reftranscripts {}
	set isocallers {}
	set organelles {}
	set realign 1
	set cleanup 1
	set paired 1
	set samBQ 0
	set adapterfile {}
	set removeskew {}
	set dt {}
	set targetfile {}
	set reports basic
	# if not set (stays empty, will be set later depending on type of experiment)
	set removeduplicates {}
	set amplicons {}
	set threads 2
	set distrreg 0
	set keepsams 0
	set datatype {}
	set aliformat bam
	set maxfastqdistr {}
	set hap_bam 0
	set depth_histo_max 1000
	set fastqdir {}
	set singlecell {}
	set sc_whitelist {}
	set sc_umisize 10
	set sc_barcodesize 16
	set sc_adaptorseq CTACACGACGCTCTTCCGATCT
	set sc_filters {}
	set sc_celltypers {}
	set sc_expectedcells {}
	set cellmarkerfile {}
	set tissue {}
	cg_options process_sample args {
		-preset {
			if {$value ne ""} {
				if {![command_exists preset_$value]} {
					error "preset $value does not exist, must be one of: [presets]"
				}
				foreach {var val} [preset_$value] {
					set $var $val
				}
				set preset $value
			}
		}
		-oridir {
			set oridir $value
		}
		-fastqdir {
			set fastqdir $value
		}
		-dbdir - -refdir {
			set dbdir $value
		}
		-minfastqreads {
			set minfastqreads $value
		}
		-clip {
			set clip $value
		}
		-p - -paired {
			set paired $value
		}
		-a - -aligner - -aligners {
			set aligners $value
		}
		-ali_keepcomments {
			set ali_keepcomments $value
		}
		-singlecell {
			if {$value ni {ontr10x {}}} {error "Unknown value $value for -singlecell, must be one of: ontr10x (or empty)"}
			set singlecell $value
		}
		-sc_whitelist {
			set sc_whitelist $value
		}
		-sc_umisize {
			set sc_umisize $value
		}
		-sc_barcodesize {
			set sc_barcodesize $value
		}
		-sc_adaptorseq {
			set sc_adaptorseq $value
		}
		-sc_filters {
			set sc_filters $value
		}
		-sc_celltypers {
			set sc_celltypers $value
		}
		-sc_expectedcells {
			set sc_expectedcells $value
		}
		-cellmarkerfile {
			if {$value ne "" && ![file exists $value]} {error "cellmarkerfile $value does not exists"}
			set cellmarkerfile [file_absolute $value]
		}
		-tissue {
			set tissue $value
		}
		-realign - -realign {
			set realign $value
		}
		-removeduplicates {
			set removeduplicates $value
		}
		-amplicons {
			set amplicons [file_absolute $value]
			if {$value ne "" && ![jobfileexists $amplicons]} {error "amplicons file $amplicons does not exists"}
		}
		-v - -varcallers {
			set varcallers $value
		}
		-svcallers {
			set svcallers $value
		}
		-methcallers {
			set methcallers $value
		}
		-counters {
			set counters $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-isocallers {
			set isocallers $value
		}
		-organelles {
			set organelles $value
		}
		-s - -split {
			set split $value
		}
		-samBQ {
			set samBQ $value
		}
		-a - -adapterfile {
			set adapterfile [file_absolute $value]
			if {$value ne "" && ![jobfileexists $adapterfile]} {error "adapterfile $adapterfile does not exists"}
		}
		-removeskew {
			set removeskew $value
		}
		-dt - -downsampling_type {
			if {$value ni "{} NONE ALL_READS BY_SAMPLE"} {error "-dt must be one of: NONE ALL_READS BY_SAMPLE"}
			set dt $value
		}
		-targetfile {
			set targetfile [file_absolute $value]
			if {$value ne "" && ![jobfileexists $targetfile]} {error "target file $targetfile does not exists"}
		}
		-r - -reports {
			set reports $value
		}
		-todoVar {
			upvar $value todo
		}
		-threads {
			set threads $value
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-c - -cleanup {
			set cleanup $value
		}
		-m - -maxopenfiles {
			maxopenfiles [expr {$value - 4}]
		}
		-keepsams {
			set keepsams $value
		}
		-datatype {
			set datatype $value
		}
		-aliformat {
			set aliformat [string tolower $value]
		}
		-maxfastqdistr {
			set maxfastqdistr $value
		}
		-hap_bam {
			set hap_bam $value
		}
		-depth_histo_max {
			set depth_histo_max $value
		}
		-extraopts {
			foreach {k v} $value {
				set ::cgextraopts($k) $v
			}
		}
		-*-* {
			set ::specialopt($key) $value
		}
	} {} 1 2
	if {[llength $args] == 1} {
		foreach {sampledir} $args break
	} elseif {[llength $args] == 2} {
		foreach {oridir sampledir} $args break
	}
	if {$ali_keepcomments eq "" && "remora" in $methcallers} {
		set ali_keepcomments 1
	}
	# If ubam dir is present, prefer this
	# we will handle ubams "as fastqs" mostly (keep in var fastqdir, etc.)
	if {$fastqdir eq ""} {
		if {[file exists $sampledir/ubam]} {
			set fastqdir $sampledir/ubam
		} else {
			set fastqdir $sampledir/fastq
		}
	}
	set reports [reports_expand $reports]
	set dbdir [file_absolute $dbdir]
	set sampledir [file_absolute $sampledir]
	adapterfile [adapterfile $adapterfile]
	if {![info exists todo]} {
		set todo(var) {}
		set todo(reports) {}
	}
	set sample [file tail $sampledir]
	if {[regexp -- - $sample]} {
		error "samplename $sample contains a -, which is not allowed; please rename the sample dir $sampledir"
	}
	if {$reftranscripts eq ""} {
		set reftranscripts [ref_tsvtranscripts $dbdir]
	}
	set organelles [getorganelles $dbdir $organelles]
	#
	putslog "Making $sampledir"
	catch {file mkdir $sampledir}
	set ref [file tail $dbdir]
	#
	# ampliconsfile
	set temp [ampliconsfile $sampledir $ref]
	if {$temp ne ""} {
		if {$amplicons ne ""} {puts stderr "Not overwriting existing ampliconsfile $temp"}
		set amplicons $temp
		if {$datatype eq ""} {set datatype amplicons}
	} elseif {$amplicons ne ""} {
		set temp [lindex [split [file root [gzroot [file tail $amplicons]]] -] end]
		mklink $amplicons $sampledir/reg_${ref}_amplicons-$temp.tsv
		set amplicons $sampledir/reg_${ref}_amplicons-$temp.tsv
		if {$datatype eq ""} {set datatype amplicons}
	}
	if {![catch {file link $amplicons} link]} {
		set ampliconsname [file tail $link]
	} else {
		set ampliconsname [file tail $amplicons]
	}
	if {$amplicons ne ""} {
		if {$removeduplicates eq ""} {set removeduplicates 0}
		if {$removeskew eq ""} {set removeskew 0}
		if {$dt eq ""} {set dt NONE}
		list_addnew dbfiles $amplicons
	} else {
		if {$removeduplicates eq ""} {set removeduplicates 1}
		if {$removeskew eq ""} {set removeskew 1}
	}
	# targetfile
	set temp [targetfile $sampledir $ref]
	if {$temp ne "" && [file exists $temp]} {
		if {$amplicons ne ""} {puts stderr "Not overwriting existing targetfile $temp"}
		set targetfile $temp
		if {$datatype eq ""} {set datatype exome}
	} elseif {$targetfile ne ""} {
		set tail [gzroot [file tail $targetfile]]
		if {[string match reg_${ref}targets*.tsv $tail] 
			|| [string match reg_targets*.tsv $tail] 
			|| [string match reg_*_targets*.tsv $tail]} {
			set filename [file tail $targetfile]
		} else {
			set temp [lindex [split [file root $tail] -] end]
			set filename reg_${ref}_targets-$temp.tsv[gzext $targetfile]
		}
		mklink $targetfile $sampledir/$filename 1
		set targetfile $sampledir/$filename
		if {$datatype eq ""} {set datatype exome}
	} elseif {$temp eq "" && $targetfile eq "" && $amplicons ne ""} {
		set targetfile $sampledir/reg_${ref}_targets.tsv.zst
		job reports_amplicons2targetfile -deps {$amplicons} -targets {$targetfile} -vars {sample dbdir ref} -code {
			cg regcollapse $dep | cg zst > $target
		}
	}
	# check projectinfo
	projectinfo $sampledir dbdir {split 1}
	set dbdir [dbdir $dbdir]
	if {[info exists oridir]} {
		if {[file exists $oridir]} {
			set oridir [file_absolute $oridir]
			file delete $sampledir/ori
			mklink $oridir $sampledir/ori 1
		}
	}
	job_logfile $sampledir/process_sample_[file tail $sampledir] $sampledir $cmdline \
		{*}[versions dbdir fastq-stats samtools gnusort8 zst os]
	# check if ori is a cg dir, if so use process_sample_cgi_job
	# ----------------------------------------------------------
	if {![job_getinfo] && [jobglob -checkcompressed 1 $sampledir/ori/ASM/var-*-ASM*.tsv] ne ""} {
		# analysis info
		# -------------
		info_analysis_file $sampledir/info_analysis.tsv $sample \
			{dbdir reports} \
			{genomecomb dbdir gnusort8 tabix zst os} \
			command [list cg process_sample {*}$keepargs]
		process_sample_cgi_job $sampledir $split
		lappend todo(var) var-cg-cg-$sample.tsv
		return $todo(var)
	}
	# analysis info
	# -------------
	if {![job_getinfo]} {
		info_analysis_file $sampledir/info_analysis.tsv $sample \
			{dbdir aligners varcallers svcallers methcallers counters reftranscripts isocallers realign paired samBQ adapterfile reports} \
			{genomecomb dbdir fastq-stats samtools gnusort8 tabix zst os} \
			command [list cg process_sample {*}$keepargs]
	}
	# convert existing vcfs
	# ----------------------
	set files [jobglob -checkcompressed 1 $sampledir/ori/var-*.vcf]
	foreach file $files {
		set target $sampledir/[file root [gzroot [file tail $file]]].tsv
		if {![file exists $target]} {
			job vcf2tsv-$file -deps {$file} -targets {$target} -vars split -code {
				set tempfile [filetemp $target]
				cg vcf2tsv -split $split $dep $tempfile
				file rename -force -- $tempfile $target
			}
			lappend todo(var) $target
		}
	}
	# add existing var files to todo(var)
	# ------------------------------
	set files [jobglob -checkcompressed 1 $sampledir/var-*.tsv]
	foreach file $files {
		set target [file root [gzroot $file]].tsv
		lappend todo(var) [file tail $target]
	}
	# get bam prefix
	# --------------
#	set keeppwd [pwd]
#	cd $sampledir
	set_job_logdir $sampledir/log_jobs
	set refseq [glob $dbdir/genome_*.ifas]
	set resultbamprefix {}
	if {$amplicons ne ""} {append resultbamprefix c}
	if {$realign} {append resultbamprefix r}
	if {$removeduplicates ni {0 {}}} {append resultbamprefix d}
	# allways sort
	append resultbamprefix s
	# single cell fastq adaptation
	# ----------------------------
	if {$singlecell ne ""} {
		# singlecell only works from fastqdir
		# check if we have the minimum number of reads required (at least 1)
		# if not, write minimum nr of reports, and return (quit processing this sample)
		set minreads [max 1 $minfastqreads]
		if {![checkminreads $fastqdir $minreads num]} {
			process_sample_reports_minfastqreads $sampledir $sample $reports $num todo
			return {}
		}
		set fastqfiles [gzfiles $fastqdir/*.fq $fastqdir/*.fastq]
		# put bams in skips (don't actually run sc_barcodes if already exis)
		set skips {}
		set skipsresult {}
		# make skips for clipping (do not do any of preliminaries if end product is already there)
		# clipped fastqs are used for all aligners!
		foreach aligner $aligners {
			set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
			set bamfile $sampledir/map-${aligner}-$sample.$aliformat
			lappend skips $bamfile $sampledir/barcode2celbarcode.tsv
			lappend skipsresult $resultbamfile $sampledir/barcode2celbarcode.tsv
			# if bam exists and is older than any of the fastqfiles -> remove (so older fastq files are not skipped)
			if {[file exists $bamfile] && (![jobtargetexists $bamfile $fastqfiles] || [file mtime $bamfile] < [file mtime [file dir [lindex $fastqfiles 0]]])} {
				putslog "$bamfile older than one of fastqfiles (renaming to .old)"
				file rename -force $bamfile $bamfile.old
			}
			if {[file exists $resultbamfile] && ![jobtargetexists $resultbamfile $fastqfiles]} {
				putslog "$resultbamfile older than one of fastqfiles (renaming to .old)"
				file rename -force $resultbamfile $resultbamfile.old
			}
		}
		sc_barcodes_job -skip $skips -skip $skipsresult \
			-whitelist $sc_whitelist \
			-umisize $sc_umisize \
			-barcodesize $sc_barcodesize \
			-adaptorseq $sc_adaptorseq \
			$fastqdir $sampledir
		set fastqdir $sampledir/bcfastq
		if {$sc_filters eq ""} {
			set sc_filters default
		}
		if {$cellmarkerfile ne ""} {
			if {$sc_celltypers eq ""} {set sc_celltypers {scsorter sctype}}
			set f [gzopen $cellmarkerfile]
			set h [tsv_open $f]
			if {$tissue ne ""} {
				if {marker ni $h} {
					error "cellmarkerfile $cellmarkerfile has no field named marker (or gene) (case sensitive)"
				}
				if {celltype ni $h} {
					error "cellmarkerfile $cellmarkerfile has no field named celltype (case sensitive)"
				}
				if {tissue ni $h} {
					error "-tissue option is given, but cellmarkerfile $cellmarkerfile has no field named tissue"
				}
			}
			gzclose $f
		} elseif {$tissue ne ""} {
			if {$sc_celltypers eq ""} {set sc_celltypers {sctype}}
		}
	}
	# use generic (fastq/bam source)
	# ------------------------------
	# find fastq files in fastq dir
	set fastqfiles [bsort [jobglob \
		$fastqdir/*.fastq.gz $fastqdir/*.fastq $fastqdir/*.fq.gz $fastqdir/*.fq \
		$fastqdir/*.bam $fastqdir/*.cram $fastqdir/*.sam \
	]]
	if {![llength $fastqfiles]} {
		file mkdir $fastqdir
		# if there are none in the fastq dir, check ori dir
		set fastqfiles [bsort [jobglob $sampledir/ori/*.fastq.gz $sampledir/ori/*.fastq $sampledir/ori/*.fq.gz $sampledir/ori/*.fq]]
		if {[llength $fastqfiles]} {
			set targets {}
			foreach file $fastqfiles {
				lappend targets $fastqdir/[file tail $file]
			}
			job fastq_from_ori-$sample -deps $fastqfiles -targets $targets -code {
				foreach file $deps target $targets {
					mklink $file $target
				}
			}
		} else {
			# check if there are bam files in ori to extract fastq from
			set files [bsort [jobglob $sampledir/ori/*.bam $sampledir/ori/*.cram]]
			foreach file $files {
				set base $fastqdir/[file tail [file root $file]]
				set target $base-R1.fastq.gz
				set target2 $base-R2.fastq.gz
				job [job_relfile2name bam2fastq- $file] -deps {$file} -cores $threads \
				-targets {$target $target2} -vars {threads} -code {
					cg bam2fastq -threads $threads $dep $target.temp.gz $target2.temp.gz
					file rename -force -- $target.temp.gz $target
					file rename -force -- $target2.temp.gz $target2
				}
			}
		}
		set fastqfiles [bsort [jobglob $fastqdir/*.fastq.gz $fastqdir/*.fastq $fastqdir/*.fq.gz $fastqdir/*.fq]]
	}
	set cleanedbams {}
	if {![llength $fastqfiles]} {
		# check for existing bam files
		foreach aligner $aligners {
			set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
			if {[file exists $resultbamfile]} {
				lappend cleanedbams $resultbamfile
			}
		}
	}
	# put check here, because fastqs might be generated from bams, etc.
	if {$singlecell eq "" && $minfastqreads > 0 && ![llength $cleanedbams]} {
		# check if we have the minimum number of reads required (default 1)
		# if not, write minimum nr of reports, and return (quit processing this sample)
		if {![checkminreads $fastqdir $minfastqreads num]} {
			process_sample_reports_minfastqreads $sampledir $sample $reports $num todo
			return {}
		}
	}
	# create bam from fastq files (if found)
	if {[llength $fastqfiles]} {
		set processlist {}
		# check even number of fastqs for paired analysis
		if {[file ext [lindex $fastqfiles 0]] in ".bam .cram .sam"} {set ubams 1} else {set ubams 0}
		if {$paired && !$ubams} {
			if {[expr {[llength $fastqfiles]%2}] != 0} {
				error "paired analysis (default, use -paired 0 to turn off), but number of fastqs is uneven ([llength $fastqfiles]) for sample $sampledir"
			}
		}
		if {[isint $maxfastqdistr]} {
			set len [llength $fastqfiles]
			if {$paired && !$ubams} {
				set perbatch [expr {($len + $maxfastqdistr -1)/$maxfastqdistr}]
				if {[expr {$perbatch%2}]} {incr perbatch}
			} else {
				set perbatch [expr {($len + $maxfastqdistr -1)/$maxfastqdistr}]
			}
			set pos 0
			for {set pos 0} {$pos < $len} {incr pos $perbatch} {
				lappend processlist [lrange $fastqfiles $pos [expr {$pos + $perbatch - 1}]]
			}
		} else {
			if {$paired && !$ubams} {
				foreach {fastq1 fastq2} $fastqfiles {
					lappend processlist [list $fastq1 $fastq2]
				}
			} else {
				foreach {fastq} $fastqfiles {
					lappend processlist [list $fastq]
				}
			}
		}
		if {[llength $processlist]} {
			set skips {}
			set skipsresult {}
			# make skips for clipping (do not do any of preliminaries if end product is already there)
			# clipped fastqs are used for all aligners!
			foreach aligner $aligners {
				set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
				set bamfile $sampledir/map-${aligner}-$sample.$aliformat
				lappend skips $bamfile
				lappend skipsresult $resultbamfile
				# if bam exists and is older than any of the fastqfiles -> remove (so older fastq files are not skipped)
				if {[file exists $bamfile] && (![jobtargetexists $bamfile $fastqfiles] || [file mtime $bamfile] < [file mtime [file dir [lindex $fastqfiles 0]]])} {
					putslog "$bamfile older than one of fastqfiles (renaming to .old)"
					file rename -force $bamfile $bamfile.old
				}
				if {[file exists $resultbamfile] && ![jobtargetexists $resultbamfile $fastqfiles]} {
					putslog "$resultbamfile older than one of fastqfiles (renaming to .old)"
					file rename -force $resultbamfile $resultbamfile.old
				}
			}
			unset -nocomplain partsa
			foreach pfastqfiles $processlist {
				set files $pfastqfiles
				set cleanupfiles {}
				set cleanupdeps {}
				# pbase takes the "root" of the filename, and will be the basis for naming results
				set pbase [file_root [file tail [lindex $pfastqfiles 0]]]
				if {[llength $pfastqfiles] > 1} {
					# if we have 2 fastqfiles (paired), and the second has a different root, add this to pbase
					set last [file_root [file tail [lindex $pfastqfiles end]]]
					if {$last ne $pbase} {
						set pos 0
						string_foreach c1 $pbase c2 $last {
							if {$c1 ne $c2} break
							incr pos
						}
						incr pos -1
						set size [string length $pbase]
						set lastlen [string length $last]
						# -2 for cnnector (__), -30 for extensions and prefix
						set minpos [expr {$lastlen - (255-2-30-$size)}]
						if {$minpos > $pos} {set pos $minpos}
						append pbase __[string range $last $pos end]
					}
				}
				foreach aligner $aligners {
					set bamfile $sampledir/map-${aligner}-$sample.$aliformat
					set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
					set target $resultbamfile.temp/[file_root [file tail [lindex $files 0]]].sam.zst
					lappend cleanupdeps $resultbamfile.temp/$pbase.sam.zst
				}
				if {$clip} {
					if {$ubams} {error "clipping not supported for ubams"}
					set files [fastq_clipadapters_job -adapterfile $adapterfile -paired $paired \
						-skip $skips -skip $skipsresult -skip $cleanupdeps \
						-removeskew $removeskew \
						{*}$files]
					foreach file $files {
						lappend cleanupfiles $file [analysisinfo_file $file]
					}
					lappend cleanupfiles [file dir [lindex $files 0]]
				}
				foreach aligner $aligners {
					# alignment per fastq per aligner
					# do not do any of preliminaries if end product is already there
					set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
					set bamfile $sampledir/map-${aligner}-$sample.$aliformat
					file mkdir $resultbamfile.temp
					job_cleanup_add $resultbamfile.temp
					set target $resultbamfile.temp/$pbase.sam.zst
					lappend partsa($aligner) $target
					# map using ${aligner}
					set opts {}
					if {[regexp {^(.*)_([^_]+)$} $aligner tmp aliprog alipreset]} {
						lappend opts -preset $alipreset
					} else {
						set aliprog $aligner
					}
					map_job -method $aliprog {*}$opts -threads $threads \
						-skip $skips \
						-joinfastqs 1 \
						-skip $skipsresult \
						-paired $paired \
						-sort coordinate \
						-compressionlevel 1 \
						-ali_keepcomments $ali_keepcomments \
						$target $refseq $sample {*}$files
				}
				if {$cleanup} {
					# clean up no longer needed intermediate files
					cleanup_job -skip $skips -skip $skipsresult [job_relfile2name cleanupclipped- $target] $cleanupfiles $cleanupdeps
				}
			}
			set cleanupdeps {}
			foreach aligner $aligners {
				# mergesort sams from individual fastq files
				# and distribute again over if -distrreg i set (for distributed cleaning)
				# do not do any of preliminaries if end product is already there
				set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
				set bamfile $sampledir/map-${aligner}-$sample.$aliformat
				if {$distrreg in {0 ""}} {
					set tempbamfile $bamfile
				} else {
					file mkdir $resultbamfile.temp
					job_cleanup_add $resultbamfile.temp
					set tempbamfile $resultbamfile.temp/map-${aligner}-$sample.$aliformat
				}
				set compressionlevel [defcompressionlevel 5]
				setdefcompressionlevel 1
				set udistrreg [distrreg_use $distrreg chr chr]
				set bamfiles [sam_catmerge_job \
					-skip [list $bamfile [analysisinfo_file $bamfile]] \
					-skip [list $resultbamfile [analysisinfo_file $resultbamfile]] \
					-name mergesams-$aligner-$sample -refseq $refseq \
					-sort coordinate -mergesort 1 \
					-distrreg $udistrreg \
					-deletesams [string is false $keepsams] -threads $threads \
					$tempbamfile {*}$partsa($aligner)]
				setdefcompressionlevel $compressionlevel
				
				# clean bamfile (mark duplicates, realign)
				# bam is already sorted, just add the s (-sort 2)
				if {[llength $bamfiles] == 1} {
					set cleanbam [bam_clean_job -sort 2 -outputformat $aliformat -distrreg $distrreg \
						-removeduplicates $removeduplicates -clipamplicons $amplicons -realign $realign \
						-regionfile 5 -refseq $refseq -threads $threads \
						 $bamfile]
					lappend cleanedbams $cleanbam
				} else {
					# distributed cleaning
					set cleanbams {}
					set compressionlevel [defcompressionlevel 5]
					setdefcompressionlevel 1
					foreach bam $bamfiles {
						lappend cleanbams [bam_clean_job -sort 2 -outputformat sam.zst -distrreg $distrreg \
							-skip [list $resultbamfile [analysisinfo_file $resultbamfile]] \
							-removeduplicates $removeduplicates -clipamplicons $amplicons -realign $realign \
							-regionfile 5 -refseq $refseq -threads $threads \
							 $bam]
					}
					# join $cleanbams \n
					setdefcompressionlevel $compressionlevel
					sam_catmerge_job \
						-name merge2bam-$aligner-$sample -refseq $refseq \
						-sort nosort \
						-index 1 \
						-deletesams [string is false $keepsams] -threads $threads \
						$resultbamfile {*}$cleanbams
					lappend cleanedbams $resultbamfile
				}
			}
		} else {
			# check if result bam file exists (without fastq), and use that if so
			foreach aligner $aligners {
				set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
				if {[file exists $resultbamfile]} {
					lappend cleanedbams $resultbamfile
				} else {
					putslog "no fastqs for $sample and file $resultbamfile does not exist"
				}
			}
		}
	}
	if {$singlecell ne ""} {
		# clean up bcfastq
		set bcfastqs [jobglob -checkcompressed 1 $sampledir/bcfastq/*]
		if {[llength $bcfastqs]} {
			job delete_bcfastq-$sample -optional 1 \
			-deps $cleanedbams -vars {sampledir} -rmtargets $bcfastqs -code {
				shadow_delete $sampledir/bcfastq
			}
		}
	}
	# make sure the bams are indexed
	foreach bam $cleanedbams {
		bam_index_job $bam
	}
	# varcaller from bams
	foreach cleanedbam $cleanedbams {
		set bambase [file_rootname $cleanedbam]
		# make 5x coverage regfile from cleanedbam
		set cov5reg [bam2reg_job -mincoverage 5 -distrreg $distrreg -refseq $refseq $cleanedbam]
		# make 20x coverage regfile
		set cov20reg [bam2reg_job -mincoverage 20 -compress 1 -distrreg $distrreg -refseq $refseq $cleanedbam]
		if {$amplicons eq ""} {
			set regionfile $cov5reg
		} else {
			set regionfile $amplicons
		}
		foreach varcaller $varcallers {
			switch $varcaller {
				gatk {set extraopts [list -dt $dt]}
				sam {set extraopts [list -dt $dt]}
				longshot {set extraopts [list -hap_bam $hap_bam]}
				default {set extraopts {}}
			}
			if {![auto_load var_${varcaller}_job]} {
				error "varcaller $varcaller not supported"
			}
			lappend cleanupdeps {*}[var_job -method ${varcaller} -distrreg $distrreg -datatype $datatype -regionfile $regionfile -split $split -threads $threads {*}$extraopts -cleanup $cleanup $cleanedbam $refseq]
			lappend todo(var) var-$varcaller-$bambase.tsv
		}
		foreach svcaller $svcallers {
			switch $svcaller {
				default {set extraopts {}}
			}
			if {[regexp {^(.*)_([^_]+)$} $svcaller tmp svcallerprog preset]} {
				lappend extraopts -preset $preset
			} else {
				set svcallerprog $svcaller
			}
			if {![auto_load sv_${svcallerprog}_job]} {
				error "svcaller $svcaller not supported"
			}
			lappend cleanupdeps {*}[sv_job -method ${svcallerprog} -distrreg $distrreg -regionfile $regionfile -split $split -threads $threads {*}$extraopts -cleanup $cleanup -refseq $refseq $cleanedbam]
			lappend todo(sv) [lindex [jobglob -checkcompressed 1 [file dir $cleanedbam]/sv-$svcaller-$bambase.tsv] 0]
		}
		foreach methcaller $methcallers {
			set extraopts {}
			if {[regexp {^(.*)_([^_]+)$} $methcaller tmp methcallerprog preset]} {
				lappend extraopts -preset $preset
			} else {
				set methcallerprog $methcaller
			}
			switch $methcallerprog {
				default {}
			}
			if {![auto_load meth_${methcallerprog}_job]} {
				error "methcaller $methcallerprog not supported"
			}
			set fast5dir [file dir $cleanedbam]/fast5
			set fastqdir [file dir $cleanedbam]/fastq
			lappend cleanupdeps {*}[meth_${methcallerprog}_job \
				-distrreg $distrreg -threads $threads {*}$extraopts -refseq $refseq \
				$fast5dir $fastqdir $cleanedbam]
			lappend todo(meth) [lindex [jobglob -checkcompressed 1 [file dir $cleanedbam]/meth-$methcaller-$bambase.tsv] 0]
		}
		foreach counter $counters {
			if {![auto_load count_${counter}_job]} {
				error "counter $counter not supported"
			}
			lappend cleanupdeps {*}[count_${counter}_job -threads $threads \
				-refseq $refseq $cleanedbam]
		}
		foreach isocaller $isocallers {
			set preset {}
			foreach {isocaller preset} [split $isocaller _] break
			if {![auto_load iso_${isocaller}_job]} {
				error "isocaller $isocaller not supported"
			}
			set options {}
			if {$preset ne ""} {lappend options -preset $preset}
			iso_${isocaller}_job \
				-reftranscripts $reftranscripts \
				{*}$options \
				-organelles $organelles \
				-cleanup $cleanup \
				-distrreg $distrreg -threads $threads \
				-refseq $refseq \
				$cleanedbam
		}
	}
	set sampledir [file_absolute $sampledir]
	set scgenefiles [jobgzfiles $sampledir/sc_gene_counts_raw-*.tsv]
	foreach scgenefile $scgenefiles {
		set scisoformfile [file dir $scgenefile]/[regsub ^sc_gene_counts_raw- [file tail $scgenefile] sc_isoform_counts_raw-]
		foreach sc_filter $sc_filters {
			if {![auto_load sc_filter_${sc_filter}_job]} {
				error "sc_filter $sc_filter not supported"
			}
			sc_filter_${sc_filter}_job \
				-reftranscripts $reftranscripts \
				-organelles $organelles \
				$scgenefile $scisoformfile $sc_expectedcells
		}
	}
	set scgenefiles [jobgzfiles $sampledir/sc_gene_counts_filtered-*.tsv]
	foreach scgenefile $scgenefiles {
		set scisoformfile [file dir $scgenefile]/[regsub ^sc_gene_counts_filtered- [file tail $scgenefile] sc_isoform_counts_filtered-]
		foreach sc_celltyper $sc_celltypers {
			if {![auto_load sc_celltyper_${sc_celltyper}_job]} {
				error "sc_celltyper $sc_celltyper not supported"
			}
			sc_celltyper_${sc_celltyper}_job -cellmarkerfile $cellmarkerfile -tissue $tissue $scgenefile $scisoformfile
		}
	}
	#calculate reports
	if {$singlecell eq "ontr10x"} {
		scywalker_report_job $sampledir $refseq
	}
	if {[llength $reports]} {
		process_reports_job -paired $paired -depth_histo_max $depth_histo_max -threads $threads $sampledir $dbdir $reports
		lappend todo(reports) $sampledir/reports
	}
	return $todo(var)
}

proc cg_process_sample {args} {
	set args [job_init {*}$args]
	process_sample_job {*}$args
	job_wait
}
