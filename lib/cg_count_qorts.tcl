proc count_mem_qorts {mem threads preset} {
	return 1G
}

# converts an ifas file (fasta with each sequence on one long line) to fasta with sequence lines split
# was put here because QoRTs fails when using a reference sequence with long sequence lines
# not currently using it (not giving reference)
proc cg_ifas2fas {ifasfile fasfile} {
	set f [open $ifasfile]
	set o [open $fasfile w]
	while 1 {
		if {[eof $f]} break
		set name [gets $f]
		set seq [gets $f]
		set len [string length $seq]
		puts "$name ($len)"
		puts $o $name
		for {set i 0} {$i < $len} {incr i 100} {
			puts $o [string range $seq $i [expr {$i+99}]]
		}
	}
	close $o
	close $f
}

proc version_qorts {} {
	set version ?
	set jar [findjar QoRTs]
	catch {exec java -jar $jar} c
	regexp {version: ([^\n]+)} $c temp version
	return $version
}

proc count_qorts_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg count_qorts_job {*}$args]
	set stranded {}
	set paired 1
	set keepargs $args
	set threads 2
	set mem 1G
	set preset {}
	set gtffile {}
	set resultfile {}
	set addfunctions {annotatedSpliceExonCounts,FPKM,calcDetailedGeneCounts,writeDocs}
	cg_options count_qorts args {
		-preset {
			set preset $value
		}
		-stranded {
			set stranded $value
		}
		-paired {
			set paired $value
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
		-addfunctions {
			set addfunctions $value
		}
		-mem {
			set mem $value
		}
	} {bamfile resultfile} 1 2 {
		count reads for genes in rna-seq experiment
	}
	set extraopts {}
	if {!$paired} {
		lappend extraopts --singleEnded
	}
	if {$stranded ne ""} {
		lappend extraopts --stranded
	}
	foreach {key value} [specialopts -qorts] {
		switch $key {
			default {
				if {$value eq ""} {
					lappend extraopts $key
				} else {
					lappend extraopts --$key $value
				}
			}
		}
	}
	set bamfile [file_absolute $bamfile]
	set bamtail [file tail $bamfile]
	set refseq [refseq $refseq]
	if {$gtffile eq ""} {
		set gtffile [lindex [jobgzfiles [bsort [file dir $refseq]/extra/gene*gencode*.gtf]] end]
	}
	if {$gtffile eq ""} {
		set gtffile [lindex [jobgzfiles [bsort [file dir $refseq]/extra/gene*.gtf]] end]
	}
	if {$gtffile eq ""} {
		set gtffile [lindex [jobgzfiles [bsort [file dir $refseq]/extra/gencode*.gtf]] end]
	}
	set flatgfffile [file dir $gtffile]/flat[file root [file tail $gtffile]].gff
	if {![file exists $flatgfffile]} {
		putslog "QoRTs: Making flat gfffile $flatgfffile needed"
		set jar [findjar QoRTs]
		catch_exec java -Xmx8G -XX:ParallelGCThreads=1 -jar $jar makeFlatGff \
			$gtffile $flatgfffile
	}
	if {$resultfile eq ""} {
		set root [file_rootname $bamtail]
		set resultfile [file dir $bamfile]/gene_counts-qorts-$root.tsv
	} else {
		set resultfile [file_absolute $resultfile]
	}
	set resultdir [file dir $resultfile]
	set root [file_rootname $resultfile]
	set qortsdir [file dir $resultfile]/$root
	set rqortsdir [file dir $resultfile]/reports-$root
	job_logfile $resultdir/count_qorts_$root $resultdir $cmdline \
		{*}[versions qorts os]
	#
	job qorts-$root -deps {
		$bamfile
		$gtffile
		$refseq
	} -targets {
		$qortsdir
		$resultfile
		$resultdir/exon_counts-$root.tsv
		$resultdir/junction_counts-$root.tsv
		$resultfile.analysisinfo
		$resultdir/exon_counts-$root.tsv.analysisinfo
		$resultdir/junction_counts-$root.tsv.analysisinfo
	} -vars {
		bamfile resultfile extraopts resultdir root qortsdir refseq
		gtffile flatgfffile addfunctions
	} -code {
		analysisinfo_write $bamfile $resultfile \
			counter qorts counter_version [version qorts] \
			reference [file2refname $refseq] \
			transcriptsfile [file tail $flatgfffile]
		analysisinfo_write $bamfile $resultdir/exon_counts-$root.tsv \
			counter qorts counter_version [version qorts] \
			reference [file2refname $refseq] \
			transcriptsfile [file tail $flatgfffile]
		analysisinfo_write $bamfile $resultdir/junction_counts-$root.tsv \
			counter qorts counter_version [version qorts] \
			reference [file2refname $refseq] \
			transcriptsfile [file tail $flatgfffile]
		analysisinfo_write $bamfile $resultdir/gene_fpkm-$root.tsv \
			counter qorts counter_version [version qorts] \
			reference [file2refname $refseq] \
			transcriptsfile [file tail $flatgfffile]
		if {[catch {
			set jar [findjar QoRTs-STABLE]
		}]} {
			set jar [findjar QoRTs]
		}
		# catch {exec java -jar $jar} c
		if {[file extension $bamfile] eq ".cram"} {
			# QoRTs gives an error on cram, make temp bamfile to solve
			set tempfile [tempdir]/[file root [file tail $bamfile]].bam
			exec samtools view -b -1 -h $bamfile > $tempfile
			exec samtools index $tempfile
			set bamfile $tempfile
		}
		file delete -force $qortsdir.temp
		catch_exec java -Xmx8G -XX:ParallelGCThreads=1 -jar $jar QC \
			--generatePlots \
			--addFunctions $addfunctions \
			--flatgff $flatgfffile \
			{*}$extraopts \
			$bamfile $gtffile $qortsdir.temp >@ stdout 2>@ stderr
		file delete -force $qortsdir
		file rename $qortsdir.temp $qortsdir
		#
		file_write $resultdir/gene_counts-$root.tsv.temp [deindent {
			#filetype tsv/countfile
			#fileversion    0.99
			#fields field	number	type	description
			#fields geneid	1	String	id field
			#fields counts	1	Integer	gene counts per sample
			#fields count_cds	1	Integer	gene counts per sample
			#fields count_utr	1	Integer	gene counts per sample
			#fields count_ambig_gene	1	Integer	gene counts per sample
		}]\n
		cg select -f "geneid=\$GENEID counts-$root=\$COUNT count_cds-$root=\$COUNT_CDS count_utr-$root=\$COUNT_UTR count_ambig_gene-$root=\$COUNT_AMBIG_GENE" \
			$qortsdir/QC.geneCounts.txt.gz >> $resultdir/gene_counts-$root.tsv.temp
		result_rename $resultdir/gene_counts-$root.tsv.temp $resultdir/gene_counts-$root.tsv
		#
		file_write $resultdir/gene_fpkm-$root.tsv.temp [deindent {
			#filetype tsv/countfile
			#fileversion    0.99
			#fields field	number	type	description
			#fields geneid	1	String	id field
			#fields fpkm	1	Integer	fpkm adjusted count per sample
		}]\n
		cg select -f "geneid=\$GENEID fpkm-$root=\$FPKM" \
			$qortsdir/QC.FPKM.txt.gz >> $resultdir/gene_fpkm-$root.tsv.temp
		result_rename $resultdir/gene_fpkm-$root.tsv.temp $resultdir/gene_fpkm-$root.tsv
		# exon_counts
		file_write $resultdir/exon_counts-$root.tsv.temp [deindent {
			#filetype tsv/countfile
			#fileversion    0.99
			#fields field	number	type	description
			#fields exonid	1	String	id field
			#fields count	Integer	exon counts per sample
			exonid	count
		}]\n
		set f [gzopen $qortsdir/QC.exonCounts.formatted.for.DEXSeq.txt.gz]
		set o [open $resultdir/exon_counts-$root.tsv.temp a]
		fcopy $f $o
		close $o
		close $f
		result_rename $resultdir/exon_counts-$root.tsv.temp $resultdir/exon_counts-$root.tsv
		# junction_counts
		file_write $resultdir/junction_counts-$root.tsv.temp [deindent {
			#filetype tsv/countfile
			#fileversion    0.99
			#fields field	number	type	description
			#fields spliceName	1	String	id field
			#fields chromosome      1       String  Chromosome/Contig
			#fields strand   1       String strand (+/-)
			#fields begin   1       Integer Begin of feature (0 based - half open)
			#fields end     1       Integer End of feature (0 based - half open)
			#fields count	Integer	junction counts per sample
		}]\n
		cg select -f "spliceName chromosome=\$chrom strand begin=\$start end counts-$root=\$CT" \
			$qortsdir/QC.spliceJunctionCounts.knownSplices.txt.gz >> $resultdir/junction_counts-$root.tsv.temp
		result_rename $resultdir/junction_counts-$root.tsv.temp $resultdir/junction_counts-$root.tsv
	}
	job qorts-report-$root -deps {
		$qortsdir
	} -targets {
		$rqortsdir
	} -vars {
		bamfile resultfile resultdir root qortsdir rqortsdir refseq
		gtffile flatgfffile addfunctions
	} -code {
		file delete -force $rqortsdir.temp
		mkdir $rqortsdir.temp
		file_write $rqortsdir.temp/decoder.txt "unique.ID\n$root\n" 
		set cmd [string_change {
			library(QoRTs)
			res <- read.qc.results.data("@base.dir@", decoder.files = "@decoder@", calc.DESeq2 = FALSE, calc.edgeR = FALSE)
			get.summary.table(res, outfile = paste0("@output.template@","summaryTable.txt"))
			makeMultiPlot.all(res, outfile.dir = "@output.template@")
			makeMultiPlot.all(res, outfile.dir = "@output.template@", plot.device.name = "pdf");
		} [list \
			@base.dir@ [file dir $qortsdir]/ \
			@decoder@ $rqortsdir.temp/decoder.txt \
			@output.dir@ $rqortsdir.temp \
			@output.template@ $rqortsdir.temp/qorts-$root- \
		]]
		set cmd [string trim $cmd]
		regsub -all {\n[\t ]*} $cmd " ; " cmd
		exec [findR] -e $cmd >@ stdout 2>@ stderr
		file delete -force $rqortsdir
		file rename $rqortsdir.temp $rqortsdir
	}
	return {}
}

proc cg_count_qorts {args} {
	set args [job_init {*}$args]
	count_qorts_job {*}$args
	job_wait
}

