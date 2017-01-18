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
		file rename -force $file.gz.temp $file.gz
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
		set chromosomes [ssort -natural [list_remdup $chromosomes]]
	} else {
		set files [jobglob $workdir/coverage/coverage*.tsv]
		foreach file $files {
			lappend chromosomes [chr_clip [lindex [split [file tail $file] -] end-1]]
		}
		set chromosomes [ssort -natural [list_remdup $chromosomes]]
	}
	# start from CGI data
	job cg_svar-$sample -optional 1 -deps {ori/ASM/var-*-ASM*.tsv} -targets {svar-$sample.tsv} \
	-skip {fannotvar-$sample.tsv sreg-$sample.tsv reg_refcons-$sample.tsv reg_nocall-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} \
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
		file rename -force $target.temp $target
	}
	job cg_sgene-$sample -optional 1 -deps {ori/ASM/gene-*-ASM*.tsv} -targets {sgene-$sample.tsv} \
	-skip {fannotvar-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} -code {
		set genefile $dep
		if {[llength $genefile] != 1} {error "could not identify genefile"}
		putslog "Sort gene file ($genefile)"
		cg select -s "chromosome begin end" $genefile $target.temp
		file rename -force $target.temp $target
	}
	# annotated vars file
	job cg_annotvar-$sample -optional 1 -vars {split} \
	-deps {svar-$sample.tsv (sgene-$sample.tsv)} -targets {annotvar-$sample.tsv} \
	-skip {fannotvar-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} -code {
		putslog "Create annotated varfile $target"
		if {[file exists $dep2]} {
			cg cg2tsv -split $split -sorted 1 $dep1 $dep2 $target.temp
		} else {
			cg cg2tsv -split $split -sorted 1 $dep1 $target.temp
		}
		file rename -force $target.temp $target
	}
	# if not from CGI, we do not have svar and sgene, take from first var_* file found
	job cg_annotvar_var-$sample -optional 1 {vars_*.tsv} {annotvar-$sample.tsv} {
		gzmklink $dep $target
	}
	# if we also do not find a var_* file, take from first variant* file found
	job cg_annotvar_variant-$sample -optional 1 {variant*.tsv} {annotvar-$sample.tsv} {
		set file [gzfile $dep]
		gzmklink $file $target
	}
	# only if reg file exists, otherwise extract from svar (next)
	job cg_sreg-$sample -optional 1 {ori/ASM/reg-*-ASM*.tsv} {sreg-$sample.tsv} {
		set regfile [gzfile $dep]
		putslog "Sort region file ($regfile)"
		cg select -s "chromosome begin end" $regfile $target.temp
		file rename -force $target.temp $target
	}
	job cg_regfromsvar-$sample -optional 1 {svar-$sample.tsv} {sreg-$sample.tsv} {
		set svarfile $dep
		putslog "Extract $target from $svarfile"
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $svarfile $target.temp
		cg regjoin $target.temp > $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
	job cg_reg_refcons-$sample -optional 1 {svar-$sample.tsv} {reg_refcons-$sample.tsv} {
		putslog "Find refcons regions for $dep"
		cg refconsregions $dep > $target.temp
		file rename -force $target.temp $target
	}
	job cg_reg_nocall-$sample -optional 1 {svar-$sample.tsv} {reg_nocall-$sample.tsv} {
		putslog "Find partial no-call regions for dep"
		if {[catch {
			nocallregions $dep $target.temp
		}]} {
			puts stderr "Could not make $job (old version files ?)"
		} else {
			file rename -force $target.temp $target
		}
	}
	set files [ssort -natural [glob -nocomplain ori/ASM/REF/coverage*-chr*-*]]
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
			exiterror "No position/offset field found in $file"
		}
		foreach field $fields {
			if {$field eq "uniqueSequenceCoverage"} {set outfield coverage} else {set outfield $field}
			set finaltarget coverage/$outfield-$sample.bcol
			set tomerge {}
			set tomergebins {}
			foreach file $files chr $chrs {
				set target coverage/$outfield-$chr-$sample.bcol
				lappend tomerge $target
				lappend tomergebins $target.bin
				job cg_coverage-$outfield-$chr-$sample -deps $file \
				   -vars {sample chr field posfield} \
				   -skip {$finaltarget $finaltarget.bin.lz4 $finaltarget.bin.lz4.lz4i} \
				   -targets {$target} -code {
					# make coverage files
					set file $dep
					file mkdir coverage
					exec {*}[gzcat $file] $file | cg bcol make -n $chr -co 0 -p $posfield -t i $target $field >@ stdout 2>@ stderr
				}
			}
			job cg_coverage_merge-$outfield-$sample \
			    -vars {tomerge tomergebins finaltarget} \
			    -deps [list {*}$tomerge {*}$tomergebins] \
			    -rmtargets [list {*}$tomerge {*}$tomergebins] \
			    -targets {$finaltarget $finaltarget.bin.lz4 $finaltarget.bin.lz4.lz4i} -code {
				cg cat -c f {*}$tomerge > $finaltarget.temp
				exec cat {*}$tomergebins > $finaltarget.bin.temp
				cg_lz4 -keep 0 -i 1 -o $finaltarget.bin.lz4 $finaltarget.bin.temp
				file delete $finaltarget.bin.temp
				file rename $finaltarget.temp $finaltarget
				file delete {*}$deps
			}
		}
	}
	job cg_cpSV-$sample -optional 1 {ori/ASM/SV ^ori/ASM/SV/(.*)$} {SV/\_} {
		putslog "Copying SV"
		set targetdir [file dir $target]
		file delete -force $targetdir
		file delete -force $targetdir.temp
		file copy $dep $targetdir.temp
		file rename -force $targetdir.temp $targetdir
	}
	job cg_cgsv-$sample -optional 1 {SV/allJunctionsBeta-*.tsv} {cgsv-$sample.tsv} {
		cg convcgsv $dep $target
	}
	job cg_cgsv_alpha-$sample -optional 1 {SV/annotatedJunctionsAlpha-*.tsv} {cgsv-$sample.tsv} {
		cg convcgsv $dep $target
	}
	job cg_cpCNV-$sample -optional 1 {ori/ASM/CNV ^ori/ASM/CNV/(.*)$} {CNV/\_} {
		putslog "Copying CNV"
		set targetdir [file dir $target]
		file delete -force $targetdir
		file delete -force $targetdir.temp
		file copy $dep $targetdir.temp
		file rename -force $targetdir.temp $targetdir
	}
	job cg_cgcnv -optional 1 {CNV/cnvSegmentsBeta-*.tsv} {cgcnv-$sample.tsv} {
		cg convcgcnv $dep $target
	}
	job cg_cgcnv_diploid-$sample -optional 1 {CNV/cnvSegmentsDiploidBeta-*.tsv} {cgcnv-$sample.tsv} {
		cg convcgcnv $dep $target
	}
	job cg_cgcnv_alpha-$sample -optional 1 {CNV/cnvSegmentsAlpha-*.tsv} {cgcnv-$sample.tsv} {
		cg convcgcnv [gzfile $dep] $target
	}
	# multiarch
	job reg_cluster-$sample -optional 1 {annotvar-$sample.tsv} {reg_cluster-$sample.tsv} {
		cg clusterregions < $dep > $target.temp
		file rename -force $target.temp $target
	}
	job reg_ns-$sample -optional 1 {annotvar-$sample.tsv} {reg_ns-$sample.tsv} {
		putslog "Find regions with N's for $dep"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < $dep > $target.temp
		file rename -force $target.temp $target
	}
	job reg_lowscore-$sample -optional 1 {annotvar-$sample.tsv} {reg_lowscore-$sample.tsv} {
		set header [cg select -h $dep]
		if {[llength [list_common $header {totalScore1 totalScore2}]] == 2} {
			putslog "Find regions with lowscores for $dep"
			cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < $dep > $target.temp
			file rename -force $target.temp $target
		}
	}
	job cg_fannotvar-$sample -optional 1 {annotvar-$sample.tsv (reg_refcons-$sample.tsv) (reg_cluster-$sample.tsv) (coverage/bcol_coverage-$sample.tsv) (coverage/bcol_refscore-$sample.tsv)} {fannotvar-$sample.tsv} {
		set temp [filetemp $target]
		cg annotate $dep $temp {*}[list_remove [lrange $deps 1 end] {}]
		cg_lz4 -keep 0 -i 1 -o $target.lz4 $temp
	}
	job cg_multitechlink_var-$sample -optional 1 {fannotvar-$sample.tsv} {var-cg-cg-$sample.tsv} {
		gzmklink $dep $target
	}
	job cg_multitechlink_sreg-$sample -optional 1 {sreg-$sample.tsv} {sreg-cg-cg-$sample.tsv} {
		gzmklink $dep $target
	}
	job cg_multitechlink_coverage-$sample -optional 1 {coverage} {coverage-cg-$sample} {
		gzmklink $dep $target
	}
	job reg_covered-$sample -optional 1 {sreg-$sample.tsv} {reg-$sample.covered} {
		putslog "Genomic coverage of sequenced regions"
		cg covered $dep > $target.temp
		file rename -force $target.temp $target
	}
	file mkdir filtered
	file mkdir covered
	job cg_filteredrefcons-$sample -optional 1 -vars sample {sreg-$sample.tsv reg_refcons-$sample.tsv} {filtered/filteredrefcons-$sample.tsv covered/filteredrefcons-$sample.covered} {
		putslog "Coverage of refcons region"
		set temp1 [filetemp $target1]
		set temp2 [filetemp $target2]
		cg regsubtract $dep1 $dep2 > $temp1
		cg covered $temp1 > $temp2
		cg_lz4 -keep 0 -i 1 -o $target1.lz4 $temp1
		file rename -force $temp2 $target2
	}
	job cg_filteredns-$sample -optional 1 {sreg-$sample.tsv reg_ns-$sample.tsv} {filtered/filteredns-$sample.tsv} {
		putslog "Coverage of ns region"
		cg regsubtract $dep1 $dep2 > $target.temp
		cg_lz4 -keep 0 -i 1 -o $target.lz4 $target.temp
	}
	job cg_filteredns_covered-$sample -optional 1 {filtered/filteredns-$sample.tsv} {covered/filteredns-$sample.covered} {
		putslog "Making $target"
		set temp [filetemp $target]
		cg covered $dep > $temp
		file rename -force $temp $target
	}
	job cg_filteredlowscore-$sample -optional 1 {sreg-$sample.tsv reg_lowscore-$sample.tsv} {filtered/filteredlowscore-$sample.tsv} {
		set temp [filetemp $target]
		cg regsubtract $dep1 $dep2 > $temp
		cg_lz4 -keep 0 -i 1 -o $target.lz4 $temp
	}
	job cg_filteredlowscore_covered-$sample -optional 1 {filtered/filteredlowscore-$sample.tsv} {covered/filteredlowscore-$sample.covered} {
		set temp [filetemp $target]
		cg covered $dep > $temp
		file rename -force $temp $target
	}
	job cg_refcons_histo-$sample -optional 1 {reg_refcons-$sample.tsv} {histo-refcons-$sample.tsv} {
		putslog "Making $target"
		cg reghisto $dep > $target.temp
		file rename -force $target.temp $target
	}
	job cg_filteredcluster-$sample -optional 1 {sreg-$sample.tsv reg_cluster-$sample.tsv} {filtered/filteredcluster-$sample.tsv} {
		putslog "Coverage of clusters region"
		set temp [filetemp $target]
		cg regsubtract $dep1 $dep2 > $temp
		cg_lz4 -keep 0 -i 1 -o $target.lz4 $temp
	}
	job cg_filteredcluster_covered-$sample -optional 1 {filtered/filteredcluster-$sample.tsv} {covered/filteredcluster-$sample.covered} {
		putslog "Making $target"
		cg covered $dep > $target.temp
		file rename -force $target.temp $target
	}
	job cg_process_cleanup-$sample -optional 1 -deps {(svar-$sample.tsv) (annotvar-$sample.tsv) (sgene-$sample.tsv) fannotvar-$sample.tsv sreg-$sample.tsv reg_refcons-$sample.tsv reg_nocall-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} \
		-vars {sample} -rmtargets {svar-$sample.tsv annotvar-$sample.tsv sgene-$sample.tsv} -code {
			catch {file delete svar-$sample.tsv}
			catch {file delete annotvar-$sample.tsv}
			catch {file delete sgene-$sample.tsv}
	}
	job cg_process_summary-$sample -deps {
		fannotvar-$sample.tsv
		sreg-$sample.tsv
		reg-$sample.covered
		coverage/coverage-$sample.bcol
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
		file rename -force $target.temp $target
	}
	cd $keepdir
}

proc process_sample_job {args} {
	set dbdir {}
	set realign 1
	set paired 1
	set adapterfile {}
	set reports all
	cg_options process_sample args {
		-oridir {
			set oridir $value
		}
		-realign {
			set realign $value
		}
		-dbdir - -refdir {
			set dbdir $value
		}
		-split {
			set split $value
		}
		-paired {
			set paired $value
		}
		-adapterfile {
			set adapterfile $value
		}
		-reports {
			set reports $value
		}
		-todoVar {
			upvar $value todo
		}
		-reportstodoVar {
			upvar $value reportstodo
		}
	} {} 1 2
	if {[llength $args] == 1} {
		foreach {destdir} $args break
	} elseif {[llength $args] == 2} {
		foreach {oridir destdir} $args break
	}
	set dbdir [file_absolute $dbdir]
	set dir [file_absolute $destdir]
	set sample [file tail $dir]
	putslog "Making $dir"
	catch {file mkdir $dir}
	# check projectinfo
	projectinfo $dir dbdir {split 1}
	if {[info exists oridir]} {
		if {[file exists $oridir]} {
			set oridir [file_absolute $oridir]
			file delete $dir/ori
			mklink $oridir $dir/ori
		}
	}
	# check if ori is a cg dir, if so use process_sample_cgi_job
	# ----------------------------------------------------------
	if {[jobglob $dir/ori/ASM/var-*-ASM*.tsv] ne ""} {
		process_sample_cgi_job $dir $split
		lappend todo cg-cg-$sample
		return
	}
	# not cgi, use generic (fastq/bam source)
	# ---------------------------------------
	set keeppwd [pwd]
	cd $dir
	job_logdir $dir/log_jobs
	set refseq [glob $dbdir/genome_*.ifas]
	set resultbamprefix rds
	# convert existing vcfs
	set files [jobglob var-*.vcf]
	foreach file $files {
		set target [file root [gzroot $file]].tsv
		if {![file exists $target]} {
			job vcf2tsv-$file -deps $file -targets $target -vars split -code {
				cg vcf2tsv -split $split $dep $target.temp
				file rename -force $target.temp $target
			}
			lappend todo [string range $target 4 end-4]
		}
	}
	# add existing var files to todo
	set files [jobglob var-*.tsv]
	foreach file $files {
		set target [file root [gzroot $file]].tsv
		lappend todo [string range $target 4 end-4]
	}
	# find fastq files in fastq dir
	set files [ssort -natural [jobglob fastq/*.fastq.gz fastq/*.fastq fastq/*.fq.gz fastq/*.fq]]
	if {![llength $files]} {
		file mkdir fastq
		# if there are none in the fastq dir, check ori dir
		set files [ssort -natural [jobglob ori/*.fastq.gz ori/*.fastq ori/*.fq.gz ori/*.fq]]
		if {[llength $files]} {
			set targets [list_regsub ^ori $files fastq]
			job fastq_from_ori-$sample -deps $files -targets $targets -code {
				foreach file $deps target $targets {
					mklink $file $target
				}
			}
		} else {
			# check if there are bam files in ori to extract fastq from
			set files [ssort -natural [jobglob ori/*.bam]]
			foreach file $files {
				set base fastq/[file tail [file root $file]]
				set target $base-R1.fastq
				set target2 $base-R2.fastq
				job bam2fastq-[file tail $file] -deps {$file} \
				-targets {$target $target2} -code {
					cg bam2fastq $dep $target.temp $target2.temp
					exec gzip $target.temp
					exec gzip $target2.temp
					file rename -force $target.temp.gz $target.gz
					file rename -force $target2.temp.gz $target2.gz
				}
			}
		}
		set files [ssort -natural [jobglob fastq/*.fastq.gz fastq/*.fastq fastq/*.fq.gz fastq/*.fq]]
	}
	# process fastq files (if found)
	if {[llength $files]} {
		# do not do any of preliminaries if end product is already there
		set bamfile map-bwa-$sample.bam
		set resultbamfile map-${resultbamprefix}bwa-$sample.bam
		# quality and adapter clipping
		set files [fastq_clipadapters_job $files -adapterfile $adapterfile -paired $paired -skips [list -skip $bamfile -skip $resultbamfile]]
		#
		# map using bwa
		map_bwa_job $refseq $files $sample $paired -skips [list -skip $resultbamfile]
		job rmclipped-$sample -optional 1 -deps $files -rmtargets $files -code {
			file delete {*}$deps
		}
	}
	# extract regions with coverage >= 5 (for cleaning)
	set cov5reg [bam2reg_job map-bwa-$sample.bam 5]
	set cov5bed [gatkworkaround_tsv2bed_job $cov5reg $refseq]
	# clean bamfile (mark duplicates, realign)
	set cleanedbam [bam_clean_job map-bwa-$sample.bam $refseq $sample -removeduplicates 1 -realign $realign -bed $cov5bed]
	# make 5x coverage regfile from cleanedbam
	set cov5reg [bam2reg_job $cleanedbam 5]
	set cov5bed [gatkworkaround_tsv2bed_job $cov5reg $refseq]
	# make 20x coverage regfile
	bam2reg_job $cleanedbam 20 1
	# samtools variant calling on map-rdsbwa
	var_sam_job $cleanedbam $refseq -bed $cov5bed -split $split
	lappend todo sam-rdsbwa-$sample
	# gatk variant calling on map-rdsbwa
	var_gatk_job $cleanedbam $refseq -bed $cov5bed -split $split
	lappend todo gatk-rdsbwa-$sample
	#calculate reports
	if {[llength $reports]} {
		proces_reports_job $dir $dbdir $reports
		lappend reportstodo $dir/$sample/reports
	}
	cd $keeppwd
	return $todo
}

proc cg_process_sample {args} {
	set args [job_init {*}$args]
	process_sample_job {*}$args
	job_wait
}
