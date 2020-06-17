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
	if {![file exists $file.${field}_index]} {
		tsv_index $field $file
	}
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
	job_logdir $workdir/log_jobs
	set chromosomes {}
	set files [jobglob $workdir/ori/ASM/REF/coverage*.tsv]
	if {[llength $files]} {
		foreach file $files {
			lappend chromosomes [chr_clip [lindex [split [file tail $file] -] 1]]
		}
		set chromosomes [bsort [list_remdup $chromosomes]]
	} else {
		set files [jobglob $workdir/coverage/coverage*.tsv]
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
				job cg_coverage-cg-$sample-$outfield-$chr-$sample -deps $file \
				   -vars {sample chr field posfield} \
				   -skip [list $finaltarget $finaltarget.bin.zst $finaltarget.bin.zst.zsti] \
				   -targets {$target $target.bin} -code {
					# make coverage files
					set file $dep
					exec {*}[gzcat $file] $file | cg bcol make -n $chr -co 0 -p $posfield -t i $target $field >@ stdout 2>@ stderr
				}
			}
			job cg_coverage-cg-${sample}_merge-$outfield-$sample \
			    -vars {tomerge tomergebins finaltarget} \
			    -deps [list {*}$tomerge {*}$tomergebins] \
			    -rmtargets [list {*}$tomerge {*}$tomergebins] \
			    -targets {$finaltarget $finaltarget.bin.zst $finaltarget.bin.zst.zsti} -code {
				cg cat -c f {*}$tomerge > $finaltarget.temp
				exec cat {*}$tomergebins > $finaltarget.bin.temp
				cg_zst -keep 0 -i 1 -o $finaltarget.bin.zst $finaltarget.bin.temp
				file delete $finaltarget.bin.temp
				file rename $finaltarget.temp $finaltarget
				file delete {*}$deps
			}
		}
	}
	# variants
	job cg_svar-$sample -optional 1 -deps {ori/ASM/var-*-ASM*.tsv} -targets {svar-$sample.tsv} \
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
#	job cg_sgene-$sample -optional 1 -deps {ori/ASM/gene-*-ASM*.tsv} -targets {sgene-$sample.tsv} \
#	-skip [list var-cg-cg-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv] -code {
#		set genefile $dep
#		if {[llength $genefile] != 1} {error "could not identify genefile"}
#		putslog "Sort gene file ($genefile)"
#		cg select -s "chromosome begin end" $genefile $target.temp
#		file rename -force -- $target.temp $target
#	}
	# annotated vars file
	job cg_annotvar-$sample -optional 1 -vars {split sample} \
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
			cg annotate -analysisinfo 0 $target.temp $target.temp2 {*}$todo
			file rename -force -- $target.temp2 $target
			file delete -force $target.temp $target.temp.index $target.temp2.index
		} else {
			file rename -force -- $target.temp $target
		}
	}
#	# if not from CGI, we do not have svar and sgene, take from first var_* file found
#	job cg_annotvar_var-$sample -optional 1 {vars_*.tsv} {annotvar-$sample.tsv} {
#		gzmklink $dep $target
#	}
#	# if we also do not find a var_* file, take from first variant* file found
#	job cg_annotvar_variant-$sample -optional 1 {variant*.tsv} {annotvar-$sample.tsv} {
#		set file [gzfile $dep]
#		gzmklink $file $target
#	}
	# only if reg file exists, otherwise extract from svar (next)
	job cg_sreg-cg-cg-$sample -optional 1 -deps {ori/ASM/reg-*-ASM*.tsv} \
	-targets {sreg-cg-cg-$sample.tsv.zst} -code {
		set regfile [gzfile $dep]
		putslog "Sort region file ($regfile)"
		cg select -s "chromosome begin end" $regfile $target.temp.zst
		file rename -force -- $target.temp.zst $target
	}
	job cg_regfromsvar-$sample -optional 1 -deps {svar-$sample.tsv} \
	-targets {sreg-cg-cg-$sample.tsv.zst} -code {
		set svarfile $dep
		putslog "Extract $target from $svarfile"
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $svarfile $target.temp
		exec cg regjoin $target.temp {*}[compresspipe .zst 9] > $target.temp2
		file rename -force -- $target.temp2 $target
		file delete $target.temp
	}
	job cg_reg_refcons-$sample -optional 1 -deps {svar-$sample.tsv} \
	-targets {reg_refcons-$sample.tsv} -code {
		putslog "Find refcons regions for $dep"
		cg refconsregions $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job cg_reg_nocall-$sample -optional 1 -deps {svar-$sample.tsv} \
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
	job cg_cpSV-$sample -optional 1 -deps {ori/ASM/SV ^ori/ASM/SV/(.*)$} -targets {SV/\_} -code {
		putslog "Copying SV"
		set targetdir [file dir $target]
		set tempdir [filetemp $targetdir]
		file delete $tempdir
		file copy $dep $tempdir
		file delete -force $targetdir
		file rename -force -- $tempdir $targetdir
	}
	job cg_cgsv-$sample -optional 1 -deps {SV/allJunctionsBeta-*.tsv} -targets {cgsv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		file delete $tempfile
		convcgsv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cgsv_alpha-$sample -optional 1 -deps {SV/annotatedJunctionsAlpha-*.tsv} -targets {cgsv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		file delete $tempfile
		convcgsv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cpCNV-$sample -optional 1 -deps {ori/ASM/CNV ^ori/ASM/CNV/(.*)$} -targets {CNV/\_} -code {
		putslog "Copying CNV"
		set targetdir [file dir $target]
		set tempdir [filetemp $targetdir]
		file delete $tempdir
		file copy $dep $tempdir
		file delete -force $targetdir
		file rename -force -- $tempdir $targetdir
	}
	job cg_cgcnv -optional 1 -deps {CNV/cnvSegmentsBeta-*.tsv} -targets {cgcnv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		convcgcnv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cgcnv_diploid-$sample -optional 1 -deps {CNV/cnvSegmentsDiploidBeta-*.tsv} -targets {cgcnv-$sample.tsv} -code {
		set tempfile [filetemp $target]
		convcgcnv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	job cg_cgcnv_alpha-$sample -optional 1 {CNV/cnvSegmentsAlpha-*.tsv} {cgcnv-$sample.tsv} {
		set tempfile [filetemp $target]
		convcgcnv $dep $tempfile
		file rename -force -- $tempfile $target
	}
	# multiarch
	job reg_cluster-$sample -optional 1 -deps {annotvar-$sample.tsv} -targets {reg_cluster-$sample.tsv.zst} -code {
		cg clusterregions < $dep > $target.temp
		cg_zst $target.temp
		file rename -force -- $target.temp.zst $target
	}
	job reg_ns-$sample -optional 1 -deps {annotvar-$sample.tsv} -targets {reg_ns-$sample.tsv} -code {
		putslog "Find regions with N's for $dep"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job reg_lowscore-$sample -optional 1 -deps {annotvar-$sample.tsv} -targets {reg_lowscore-$sample.tsv} -code {
		set header [cg select -h $dep]
		if {[llength [list_common $header {totalScore1 totalScore2}]] == 2} {
			putslog "Find regions with lowscores for $dep"
			cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < $dep > $target.temp
			file rename -force -- $target.temp $target
		}
	}
	job cg_var-cg-cg-$sample -optional 1 \
	-deps {annotvar-$sample.tsv (reg_refcons-$sample.tsv) (reg_cluster-$sample.tsv) (coverage-cg-$sample/bcol_coverage-$sample.tsv) (coverage-cg-$sample/bcol_refscore-$sample.tsv)} \
	-targets {var-cg-cg-$sample.tsv.zst var-cg-cg-$sample.tsv.analysisinfo} -vars {sample} -code {
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
		cg annotate -analysisinfo 0 $dep $tempfile {*}[list_remove [lrange $deps 1 end] {}]
		analysisinfo_write $dep $target sample cg-cg-$sample aligner cgi aligner_version $cgi_version varcaller cgi varcaller_version $cgi_version reference $reference
		file rename $tempfile $target
		file rename $tempfile.zsti $target.zsti
		file delete -force [gzroot $tempfile].index
		file delete -force [gzroot $dep].index
	}
	job reg_covered-$sample -optional 1 -deps {sreg-cg-cg-$sample.tsv.zst} -targets {reg-$sample.covered} -code {
		putslog "Genomic coverage of sequenced regions"
		cg covered $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	file mkdir filtered
	file mkdir covered
	job cg_filteredrefcons-$sample -optional 1 -vars sample \
	-deps {sreg-cg-cg-$sample.tsv reg_refcons-$sample.tsv} \
	-targets {filtered/filteredrefcons-$sample.tsv covered/filteredrefcons-$sample.covered} \
	-code {
		putslog "Coverage of refcons region"
		set temp1 [filetemp $target1]
		set temp2 [filetemp $target2]
		cg regsubtract $dep1 $dep2 > $temp1
		cg covered $temp1 > $temp2
		cg_zst -keep 0 -i 1 -o $target1.zst $temp1
		file rename -force -- $temp2 $target2
	}
	job cg_filteredns-$sample -optional 1 -deps {sreg-cg-cg-$sample.tsv reg_ns-$sample.tsv} \
	-targets {filtered/filteredns-$sample.tsv.zst} -code {
		putslog "Coverage of ns region"
		cg regsubtract $dep1 $dep2 > $target.temp
		cg_zst -keep 0 -i 1 -o $target $target.temp
	}
	job cg_filteredns_covered-$sample -optional 1 -deps {filtered/filteredns-$sample.tsv} \
	-targets {covered/filteredns-$sample.covered} -code {
		putslog "Making $target"
		set temp [filetemp $target]
		cg covered $dep > $temp
		file rename -force -- $temp $target
	}
	job cg_filteredlowscore-$sample -optional 1 -deps {sreg-cg-cg-$sample.tsv reg_lowscore-$sample.tsv} \
	-targets {filtered/filteredlowscore-$sample.tsv.zst} -code {
		set temp [filetemp $target]
		cg regsubtract $dep1 $dep2 > $temp
		cg_zst -keep 0 -i 1 -o $target $temp
	}
	job cg_filteredlowscore_covered-$sample -optional 1 -deps {filtered/filteredlowscore-$sample.tsv} \
	-targets {covered/filteredlowscore-$sample.covered} -code {
		set temp [filetemp $target]
		cg covered $dep > $temp
		file rename -force -- $temp $target
	}
	job cg_refcons_histo-$sample -optional 1 -deps {reg_refcons-$sample.tsv} -targets {histo-refcons-$sample.tsv} -code {
		putslog "Making $target"
		cg reghisto $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job cg_filteredcluster-$sample -optional 1 -deps {sreg-cg-cg-$sample.tsv reg_cluster-$sample.tsv} \
	-targets {filtered/filteredcluster-$sample.tsv.zst} -code {
		putslog "Coverage of clusters region"
		set temp [filetemp $target]
		cg regsubtract $dep1 $dep2 > $temp
		cg_zst -keep 0 -i 1 -o $target $temp
	}
	job cg_filteredcluster_covered-$sample -optional 1 -deps {filtered/filteredcluster-$sample.tsv} \
	-targets {covered/filteredcluster-$sample.covered} -code {
		putslog "Making $target"
		cg covered $dep > $target.temp
		file rename -force -- $target.temp $target
	}
	job cg_process_cleanup-$sample -optional 1 \
		-deps {(svar-$sample.tsv) (annotvar-$sample.tsv) (annotvar-$sample.tsv.index) (sgene-$sample.tsv) var-cg-cg-$sample.tsv sreg-cg-cg-$sample.tsv reg_refcons-$sample.tsv reg_nocall-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} \
		-vars {sample} -rmtargets {svar-$sample.tsv annotvar-$sample.tsv sgene-$sample.tsv} -code {
			catch {file delete svar-$sample.tsv}
			catch {file delete annotvar-$sample.tsv}
			catch {file delete annotvar-$sample.tsv.index}
			catch {file delete sgene-$sample.tsv}
	}
	job cg_process_summary-$sample -deps {
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

proc process_sample_job {args} {
	upvar job_logdir job_logdir
	set cmdline [list cg process_sample {*}$args]
	set keepargs $args
	set dbdir {}
	set minfastqreads 0
	set clip 1
	set aligners bwa
	set varcallers {gatk sam}
	set svcallers {}
	set methcallers {}
	set realign 1
	set cleanup 1
	set paired 1
	set samBQ 0
	set adapterfile {}
	set removeskew {}
	set dt {}
	set targetfile {}
	set reports basic
	set removeduplicates {}
	set amplicons {}
	set threads 2
	set distrreg 0
	set keepsams 0
	set datatype {}
	set aliformat bam
	set maxfastqdistr {}
	set hap_bam 0
	cg_options process_sample args {
		-oridir {
			set oridir $value
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
	} {} 1 2
	if {[llength $args] == 1} {
		foreach {sampledir} $args break
	} elseif {[llength $args] == 2} {
		foreach {oridir sampledir} $args break
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
	#
	if {$minfastqreads > 0} {
		set files [bsort [jobglob $sampledir/fastq/*.fastq.gz $sampledir/fastq/*.fastq $sampledir/fastq/*.fq.gz $sampledir/fastq/*.fq]]
		if {![llength $files]} {return {}}
		set file [lindex $files 0]
		set f [gzopen $file]
		set count $minfastqreads
		while {$count} {
			for {set i 0} {$i < 4} {incr i} {
				if {[gets $f line] == -1} {
					gzclose $f
					if {[inlist $reports fastqstats] || [inlist $reports all] || [inlist $reports basic]} {
						set num [expr {$minfastqreads-$count}]
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
					return {}
				}
			}
			incr count -1
		}
		gzclose $f
	}
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
		mklink $amplicons $sampledir/reg_${ref}_amplicons-$temp.tsv 1
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
	} else {
		if {$removeduplicates eq ""} {set removeduplicates 1}
		if {$removeskew eq ""} {set removeskew 1}
		if {$dt eq ""} {set dt BY_SAMPLE}
	}
	# targetfile
	set temp [targetfile $sampledir $ref]
	if {$temp ne ""} {
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
		{*}[versions dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk gatk3 picard java gnusort8 zst os]
	# check if ori is a cg dir, if so use process_sample_cgi_job
	# ----------------------------------------------------------
	if {![job_getinfo] && [jobglob $sampledir/ori/ASM/var-*-ASM*.tsv] ne ""} {
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
			{dbdir aligners varcallers svcallers methcallers realign paired samBQ adapterfile reports} \
			{genomecomb dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk3 gatk gatkjava picard java gnusort8 tabix zst os} \
			command [list cg process_sample {*}$keepargs]
	}
	# convert existing vcfs
	# ----------------------
	set files [jobglob $sampledir/var-*.vcf]
	foreach file $files {
		set target [file root [gzroot $file]].tsv
		if {![file exists $target]} {
			job vcf2tsv-$file -deps {$file} -targets {$target} -vars split -code {
				cg vcf2tsv -split $split $dep $target.temp
				file rename -force -- $target.temp $target
			}
			lappend todo(var) $target
		}
	}
	# add existing var files to todo(var)
	# ------------------------------
	set files [jobglob $sampledir/var-*.tsv]
	foreach file $files {
		set target [file root [gzroot $file]].tsv
		lappend todo(var) [file tail $target]
	}
	# use generic (fastq/bam source)
	# ------------------------------
#	set keeppwd [pwd]
#	cd $sampledir
	job_logdir $sampledir/log_jobs
	set refseq [glob $dbdir/genome_*.ifas]
	set resultbamprefix {}
	if {$amplicons ne ""} {append resultbamprefix c}
	if {$realign} {append resultbamprefix r}
	if {$removeduplicates} {append resultbamprefix d}
	# allways sort
	append resultbamprefix s
	# find fastq files in fastq dir
	set fastqfiles [bsort [jobglob $sampledir/fastq/*.fastq.gz $sampledir/fastq/*.fastq $sampledir/fastq/*.fq.gz $sampledir/fastq/*.fq]]
	if {![llength $fastqfiles]} {
		file mkdir $sampledir/fastq
		# if there are none in the fastq dir, check ori dir
		set fastqfiles [bsort [jobglob $sampledir/ori/*.fastq.gz $sampledir/ori/*.fastq $sampledir/ori/*.fq.gz $sampledir/ori/*.fq]]
		if {[llength $fastqfiles]} {
			set targets {}
			foreach file $fastqfiles {
				lappend targets $sampledir/fastq/[file tail $file]
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
				set base $sampledir/fastq/[file tail [file root $file]]
				set target $base-R1.fastq.gz
				set target2 $base-R2.fastq.gz
				job bam2fastq-[file tail $file] -deps {$file} -cores $threads \
				-targets {$target $target2} -vars {threads} -code {
					cg bam2fastq -threads $threads $dep $target.temp.gz $target2.temp.gz
					file rename -force -- $target.temp.gz $target
					file rename -force -- $target2.temp.gz $target2
				}
			}
		}
		set fastqfiles [bsort [jobglob $sampledir/fastq/*.fastq.gz $sampledir/fastq/*.fastq $sampledir/fastq/*.fq.gz $sampledir/fastq/*.fq]]
	} 
	# create bam from fastq files (if found)
	set cleanedbams {}
	set processlist {}
	if {[isint $maxfastqdistr]} {
		set len [llength $fastqfiles]
		if {$paired} {
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
		if {$paired} {
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
		# clipped fastqs are used for all aligners for all aligners!
		foreach aligner $aligners {
			set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
			set bamfile $sampledir/map-${aligner}-$sample.bam
			lappend skips $bamfile [analysisinfo_file $bamfile]
			lappend skipsresult $resultbamfile [analysisinfo_file $resultbamfile]
		}
		unset -nocomplain partsa
		foreach pfastqfiles $processlist {
			set files $pfastqfiles
			set cleanupfiles {}
			set cleanupdeps {}
			set pbase [file_root [file tail [lindex $pfastqfiles 0]]]
			if {[llength $pfastqfiles] > 1} {
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
				set bamfile $sampledir/map-${aligner}-$sample.bam
				set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
				set target $resultbamfile.index/[file_root [file tail [lindex $files 0]]].sam.zst
				lappend cleanupdeps $resultbamfile.index/$pbase.sam.zst
			}
			if {$clip} {
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
				set bamfile $sampledir/map-${aligner}-$sample.bam
				file mkdir $resultbamfile.index
				set target $resultbamfile.index/$pbase.sam.zst
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
					$target $refseq $sample {*}$files
			}
			if {$cleanup} {
				# clean up no longer needed intermediate files
				cleanup_job cleanupclipped-[file tail $target] $cleanupfiles $cleanupdeps
			}
		}
		set cleanupdeps {}
		foreach aligner $aligners {
			# mergesort sams from individual fastq files
			# and distribute again over if -distrreg i set (for distributed cleaning)
			# do not do any of preliminaries if end product is already there
			set resultbamfile $sampledir/map-${resultbamprefix}${aligner}-$sample.$aliformat
			set bamfile $sampledir/map-${aligner}-$sample.bam
			if {$distrreg in {0 ""}} {
				set tempbamfile $bamfile
			} else {
				file mkdir $resultbamfile.index
				set tempbamfile $resultbamfile.index/map-${aligner}-$sample.bam
			}
			set compressionlevel [defcompressionlevel 5]
			setdefcompressionlevel 1
			set udistrreg [distrreg_use $distrreg chr chr]
			set bamfiles [sam_catmerge_job \
				-skip [list $bamfile $bamfile.analysisinfo] \
				-skip [list $resultbamfile $resultbamfile.analysisinfo] \
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
				set cleanbams {}
				set compressionlevel [defcompressionlevel 5]
				setdefcompressionlevel 1
				foreach bam $bamfiles {
					lappend cleanbams [bam_clean_job -sort 2 -outputformat sam.zst -distrreg $distrreg \
						-skip [list $resultbamfile $resultbamfile.analysisinfo] \
						-removeduplicates $removeduplicates -clipamplicons $amplicons -realign $realign \
						-regionfile 5 -refseq $refseq -threads $threads \
						 $bam]
				}
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
			switch {$varcaller} {
				gatk {set extraopts [list -dt $dt]}
				sam {set extraopts [list -dt $dt]}
				longshot {set extraopts [list -hap_bam $hap_bam}
				default {set extraopts {}}
			}
			if {![auto_load var_${varcaller}_job]} {
				error "varcaller $varcaller not supported"
			}
			lappend cleanupdeps {*}[var_job -method ${varcaller} -distrreg $distrreg -datatype $datatype -regionfile $regionfile -split $split -threads $threads {*}$extraopts -cleanup $cleanup $cleanedbam $refseq]
			lappend todo(var) var-$varcaller-$bambase.tsv
		}
		foreach svcaller $svcallers {
			switch {$svcaller} {
				default {set extraopts {}}
			}
			if {![auto_load sv_${svcaller}_job]} {
				error "svcaller $svcaller not supported"
			}
			lappend cleanupdeps {*}[sv_job -method ${svcaller} -distrreg $distrreg -regionfile $regionfile -split $split -threads $threads {*}$extraopts -cleanup $cleanup -refseq $refseq $cleanedbam]
			lappend todo(sv) sv-$svcaller-$bambase.tsv
		}
		foreach methcaller $methcallers {
			switch {$methcaller} {
				default {set extraopts {}}
			}
			if {![auto_load meth_${methcaller}_job]} {
				error "methcaller $methcaller not supported"
			}
			set fast5dir [file dir $cleanedbam]/fast5
			set fastqdir [file dir $cleanedbam]/fastq
			lappend cleanupdeps {*}[meth_${methcaller}_job \
				-distrreg $distrreg -threads $threads {*}$extraopts -refseq $refseq \
				$fast5dir $fastqdir $cleanedbam]
			lappend todo(meth) meth-$methcaller-$bambase.tsv
		}
	}
	#calculate reports
	if {[llength $reports]} {
		process_reports_job $sampledir $dbdir $reports
		lappend todo(reports) $sampledir/reports
	}
	return $todo(var)
}

proc cg_process_sample {args} {
	set args [job_init {*}$args]
	process_sample_job {*}$args
	job_wait
}
