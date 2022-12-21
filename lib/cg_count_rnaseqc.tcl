#proc rnaseqc_getref {} {
#	set tempdir [tempdir]
#	wgetfile https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/TOPMed_RNAseq_v2/gene_model/collapse_annotation.py $tempdir/collapse_annotation.py
#	exec gunzip /complgen/refseq/hg38/extra/gene_hg38_gencode.v39.gtf.gz
#	set gtf /complgen/refseq/hg38/extra/gene_hg38_gencode.v39.gtf
#	set collapsedgtf /complgen/refseq/hg38/extra/collapsedgene_hg38_gencode.v39.gtf
#	exec python3 $tempdir/collapse_annotation.py $gtf $collapsedgtf
#}

#proc rnaseqc_gtf_job {refseq {gtffile {}} {threads 4}} {
#	upvar job_logdir job_logdir
#	set refseq [refseq $refseq]
#	set gtffile [file_absolute $gtffile]
#	set rnaseqcrefseq $refseq.rnaseqc
#	if {[file exists $rnaseqcrefseq]} {return $rnaseqcrefseq}
#	if {[jobtargetexists [list $rnaseqcrefseq] $refseq]} return
#	job [job_relfile2name rnaseqc_2refseq- $refseq] -deps {
#		$refseq
#	} -targets {
#		$rnaseqcrefseq
#	} -vars {
#		refseq rnaseqcrefseq gtffile threads
#	} -code {
#		if {[gziscompressed $gtffile]} {
#			set tempfile [tempfile].gtf
#			exec cg zcat $gtffile > $tempfile
#			set gtffile $tempfile
#		}
#		file delete -force $rnaseqcrefseq.temp
#		file mkdir $rnaseqcrefseq.temp
#		set tail [file tail $refseq]
#		mklink $refseq $rnaseqcrefseq.temp/$tail
#		set extraopts {}
#		if {$gtffile ne ""} {
#			lappend extraopts --sjdbGTFfile $gtffile
#		}
#		set temp [catch_exec rnaseqc \
#			--runMode genomeGenerate \
#			--genomeFastaFiles $refseq \
#			--genomeDir $rnaseqcrefseq.temp \
#			--runThreadN $threads \
#			{*}$extraopts \
#		]
#		if {[regexp {loaded/built the index for 0 target sequence\(s\)} $temp]} {
#			error "could not properly index $dep: contains no sequences"
#		}
#		file rename -- $target.temp $target
#	}
#	return $rnaseqcrefseq
#}
#
#proc cg_rnaseqc_gtf args {
#	set args [job_init {*}$args]
#	set return [refseq_rnaseqc_job {*}$args]
#	job_wait
#	return $return
#}

proc count_mem_rnaseqc {mem threads preset} {
	return 1G
}

proc gct2tsv {gct tsv datafield datatype datadescr {idfield geneid} {namefield genename}} {
	set f [gzopen $gct]
	set o [wgzopen $tsv]
	set version [gets $f]
	set dimensions [gets $f]
	foreach {rows columns metarows metacolumns} $dimensions break
	set header [gets $f]
	if {[isint $metarows]} {
		for {set i 0} {$i < $metarows} {incr i} {
			gets $f
		}
	}
	set header [split $header \t]
	set newheader [list $idfield $namefield]
	foreach sample [lrange $header 2 end] {
		lappend newheader $datafield-$sample
	}
	puts $o [deindent [subst {
		#filetype tsv/countfile
		#fileversion    0.99
		#fields field	number	type	description
		#fields $idfield	1	String	id field
		#fields $namefield	1	String	name field
		#fields $datafield	1	$datatype	$datadescr
	}]]
	puts $o [join $newheader \t]
	fcopy $f $o
}

proc count_rnaseqc_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg count_rnaseqc {*}$args]
	set extraopts {}
	set stranded 1
	set keepargs $args
	set threads 2
	set mem 1G
	set preset {}
	set gtffile {}
	set resultfile {}
	cg_options count_rnaseqc args {
		-preset {
			set preset $value
		}
		-stranded {
			set stranded $value
		}
		-threads - -t {
			# not actually supported
			set threads $value
		}
		-refseq {
			set refseq $value
		}
		-gtffile {
			set gtffile $value
		}
		-mem {
			set mem $value
		}
	} {bamfile resultfile} 1 2 {
		count reads for genes in rna-seq experiment
	}
	set extraopts {}
	if {$preset eq "ht"} {
		lappend extraopts -q 10
	}
	foreach {key value} [specialopts -rnaseqc] {
		switch $key {
			default {
				if {$value eq ""} {
					lappend extraopts $key
				} else {
					lappend extraopts $key=$value
				}
			}
		}
	}
	set bamfile [file_absolute $bamfile]
	set bamtail [file tail $bamfile]
	set refseq [refseq $refseq]
	if {$gtffile eq ""} {
		set gtffile [lindex [jobgzfiles [bsort [file dir $refseq]/extra/collapsedgene_*gencode*.gtf]] end]
	}
	if {$gtffile eq ""} {
		set gtffile [lindex [jobgzfiles [bsort [file dir $refseq]/extra/collapsedgene*.gtf]] end]
	}
	if {$gtffile eq ""} {
		error "no gtffile found in refdir (looking for [file dir $refseq]/extra/collapsedgene_*gencode*.gtf)"
	}
	if {$resultfile eq ""} {
		set root [file_rootname $bamtail]
		set resultfile [file dir $bamfile]/gene_counts-rnaseqc-$root.tsv
	} else {
		set resultfile [file_absolute $resultfile]
	}
	set resultdir [file dir $resultfile]
	set root [file_rootname $resultfile]
	set rnaseqcdir [file dir $resultfile]/$root
	job_logfile $resultdir/count_rnaseqc_$root $resultdir $cmdline \
		{*}[versions rnaseqc os]
	#
	job rnaseqc-$root -deps {
		$bamfile
		$gtffile
		$refseq
	} -targets {
		$resultfile
		$resultdir/exon_counts-$root.tsv
		$resultdir/tpm-$root.tsv
	} -vars {
		bamfile resultfile extraopts resultdir root rnaseqcdir refseq gtffile
	} -code {
		analysisinfo_write $bamfile $resultfile counter rnaseqc counter_version [version rnaseqc] reference [file2refname $refseq]
		putslog "making $resultfile"
		catch_exec rnaseqc $gtffile $bamfile $rnaseqcdir \
			--sample=$root \
			--fasta=$refseq \
			--unpaired \
			--coverage \
			--detection-threshold=3 \
			{*}$extraopts
		gct2tsv $rnaseqcdir/$root.exon_reads.gct $rnaseqcdir/$root.exon_reads.tsv counts Float "exon counts per sample" exon genename
		gct2tsv $rnaseqcdir/$root.gene_reads.gct $rnaseqcdir/$root.gene_reads.tsv counts Integer "gene counts per sample"
		gct2tsv $rnaseqcdir/$root.gene_tpm.gct $rnaseqcdir/$root.gene_tpm.tsv tpm Float "transcripts per million"
		mklink $rnaseqcdir/$root.exon_reads.tsv $resultdir/exon_counts-$root.tsv
		mklink $rnaseqcdir/$root.gene_tpm.tsv $resultdir/tpm-$root.tsv
		mklink $rnaseqcdir/$root.gene_reads.tsv $resultfile
	}
	return {}
}

proc cg_count_rnaseqc {args} {
	set args [job_init {*}$args]
	count_rnaseqc_job {*}$args
	job_wait
}
