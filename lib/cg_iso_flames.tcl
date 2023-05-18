proc version_flames {} {
	set bin [follow_links [exec which flames]]
	set flamesdir [file dir $bin]
	lindex [split [file tail $flamesdir] -] end
}

proc iso_flames_header_genecounts {root} {
	set header [deindent {
		#filetype	tsv/genecountsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	gene	1	String	gene id (usually gene_id from GTF)
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	category	1	String	gene category (known or novel)
		#fields	nrtranscripts	1	Integer	number of transcripts found
		#fields	counts	1	Integer	Number of reads mapping to gene
		#fields	type	1	String	Type of element
	}]
	set oheader {gene chromosome begin end strand type category nrtranscripts}
	lappend oheader counts-$root
	append header \n[join $oheader \t]
	return $header
}

proc cg_iso_flames_genecounts {isoformcounts genecounts {root {}}} {
	set tempfile [tempfile]
	if {$root eq ""} {
		set root [file_rootname $genecounts]
	}
	cg select -overwrite 1 \
		-g {geneid} \
		-gc {distinct(chromosome),min(begin),max(end),distinct(strand),distinct(structural_category),ucount(transcript),sum(counts-*)} \
		$isoformcounts $tempfile
	catch {close $f} ; catch {close $o}
	set f [open $tempfile]
	set inheader [tsv_open $f]
	set o [open $genecounts w]
	puts $o [iso_flames_header_genecounts $root]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set line [linsert $line 5 gene]
		if {[regexp intergenic [lindex $line 6]]} {set category novel} else {set category known}
		lset line 6 $category
		puts $o [join $line \t]
	}
	close $o
	close $f
}

proc iso_flames_header_isoformcounts {root} {
	set header [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	type	1	String	type of element
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	score	1	Float	Score
		#fields	geneid	1	String	gene id (e.g. gene_id from GTF)
		#fields	cdsStartStat	1	String	Status of CDS start annotation (none, unknown, incomplete, or complete)
		#fields	cdsEndStat	1	String	Status of CDS end annotation (none, unknown, incomplete, or complete)
		#fields	exonFrames	E	Integer	Exon frame offsets {0,1,2}
		#fields	category	1	String	one of the isoform categories (known, novel, intergenic)
		#fields	counts	1	Integer	Number of reads mapping to isoform
	}]
	append header \n[join [list \
		transcript chromosome strand begin end cdsStart cdsEnd exonCount exonStarts exonEnds score geneid cdsStartStat cdsEndStat exonFrames \
		type category counts-$root
	] \t]
	return $header
}

proc iso_flames_job {args} {
# putslog [list flames_job {*}$args]
	upvar job_logdir job_logdir
	global appdir
	set cmdline [clean_cmdline cg flames {*}$args]
	set refseq {}
	set skips {}
	set genes {}
	set extraopts {}
	set preset {}
	set referencecsv {}
	set hasumi {}
	set sqanti 1
	set threads 8
	set resultfile {}
	set reftranscripts {}
	cg_options flames args {
		-refseq {
			set refseq $value
		}
		-threads {
			set threads $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-distrreg {
			# not used
		}
		-preset {
			if {$preset ni "sc {} normal"} {error "isoflames only allows -preset sc or '' (empty)"}
			set preset $value
		}
		-referencecsv {
			set referencecsv $value
		}
		-hasumi {
			set hasumi $value
		}
		-extraopts {
			set extraopts $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {bam resultfile} 1 2
	set bam [file_absolute $bam]
	if {[file isdir $bam]} {
		# parameter given is a dir, suppose it is the fastqdir directly
		set fastqdir $bam
	} else {
		# flames needs the fastqs, find them for the given bam file
		set fastqdir [file dir $bam]/fastq
		if {![file exists $fastqdir]} {
			error "could not find fastq dir for $bam ($fastqdir)"
		}
	}
	if {$preset in {{} normal}} {
		set analysisname flames
	} else {
		set analysisname flames_$preset
	}
	if {$resultfile eq ""} {
		set sample [file tail [file dir $fastqdir]]
		set root $analysisname-fastqs-$sample
		set resultfile [file dir $fastqdir]/isoform_counts-$root.tsv
	} else {
		set resultfile [file_absolute $resultfile]
		set root [file_rootname $resultfile]
		set sample [lindex [split $root -] end]
	}
	set resultfile [file_absolute $resultfile]
	set reftranscripts [file_absolute $reftranscripts]
	set resulttail [file tail $resultfile]
	set destdir [file dir $resultfile]
	if {$hasumi eq ""} {
		if {$preset eq "sc"} {set hasumi true} else {set hasumi false}
	}
	if {[true $hasumi]} {
		set hasumi true
	} else {
		set hasumi false
	}
	#
	set refseq [refseq $refseq]
	if {$reftranscripts eq ""} {
		set reftranscripts [ref_gtftranscripts $refseq]
	} else {
		set reftranscripts [file_absolute $reftranscripts]
	}
	# 
	job_logfile $destdir/flames_[file tail $destdir] $destdir $cmdline \
		{*}[versions flames dbdir zstd os]
	# analysis per sample
	# set flamesdir $destdir/$root
	set fastqfiles [bsort [jobgzfiles $fastqdir/*.fq $fastqdir/*.fastq]]
	if {[llength $fastqfiles] == 0} {
		error "could not run flames: no fastq files found (tried in $fastqdir)"
	}
	job flames-[file tail $resultfile] {*}$skips \
	-cores $threads \
	-deps [list {*}$fastqfiles $refseq $reftranscripts $referencecsv] -targets {
		$resultfile
		$destdir/gene_counts-$root.tsv
	} -vars {
		fastqdir fastqfiles resultfile reftranscripts root destdir refseq threads extraopts sample
		preset referencecsv
	} -code {
		set extrainfo [list \
			analysis $root sample $sample \
			isocaller_reftranscripts [file tail $reftranscripts] \
			isocaller_distrreg 0 \
			isocaller flames isocaller_version [version flames] \
			isocaller_preset $preset \
		]
		analysisinfo_write [gzfile $fastqdir/*.fastq $fastqdir/*.fq] $resultfile \
			{*}$extrainfo
		analysisinfo_write [gzfile $fastqdir/*.fastq $fastqdir/*.fq] $destdir/gene_counts-$root.tsv \
			{*}$extrainfo
		set keeppwd [pwd]
		set flamesdir [tempdir]
		mkdir $flamesdir
		cd $flamesdir
		file_write config.json [deindent [subst {
			{
			    "comment":"this is the default config for SIRV spike-in data. use splice annotation on alignment.",
			    "pipeline_parameters":{
			        "do_genome_alignment":true,
			        "do_isoform_identification":true,
			        "do_read_realignment":true,
			        "do_transcript_quantification":true
			    },
			    "global_parameters":{
			        "generate_raw_isoform":true,
			        "has_UMI":$hasumi
			    },
			    "isoform_parameters":{
			        "MAX_DIST":10,
			        "MAX_TS_DIST":100,
			        "MAX_SPLICE_MATCH_DIST":10,
			        "min_fl_exon_len":40,
			        "Max_site_per_splice":3,
			        "Min_sup_cnt":10,
			        "Min_cnt_pct":0.01,
			        "Min_sup_pct":0.2,
			        "strand_specific":0,
			        "remove_incomp_reads":5
			    },
			    "alignment_parameters":{
			        "use_junctions":true,
			        "no_flank":true
			    },
			    "realign_parameters":{
			        "use_annotation":true
			    },
			    "transcript_counting":{
			        "min_tr_coverage":0.75,
			        "min_read_coverage":0.75
			    }
			}
		}]]
		mklink $refseq [file tail $refseq]
		mklink $reftranscripts [file tail $reftranscripts]
		mkdir fastq
		set mergedfastqfile fastq/$root.fastq
		if {$preset eq "sc"} {
			exec flames match_cell_barcode \
				$fastqdir \
				barcodestatistics-$root.tsv \
				$mergedfastqfile.gz \
				$referencecsv \
				2
		} else {
			if {[llength $fastqfiles] > 1} {
				exec cg zcat {*}$fastqfiles | gzip --fast > $mergedfastqfile.gz.temp
				file rename -force $mergedfastqfile.gz.temp $mergedfastqfile.gz
				set mergedfastqfile $mergedfastqfile.gz
			} else {
				set fastqfile [lindex $fastqfiles 0]
				mklink -absolute 0 $fastqfile $mergedfastqfile[gzext $fastqfile]
				set mergedfastqfile $mergedfastqfile[gzext $fastqfile]
			}
		}
		set f [gzopen $mergedfastqfile]
		for {set i 0} {$i < 20} {incr i} {gets $f}
		set short [eof $f]
		if {$short} {
			# only 4 reads aligned -> skip running flames
			file_write isoform_annotated.filtered.gff3 {}
			file_write transcript_assembly.fa {}
			set o [wgzopen transcript_count.csv.gz w] ; gzclose $o
			file_write flames-not_enough_reads ""
			cd $destdir
			file_write isoform_counts-${root}.tsv [iso_flames_header_isoformcounts $root]\n
			set inheader [list gene chromosomeb begin end strand type category nrtranscripts sum_counts-$root]
			file_write gene_counts-${root}.tsv [iso_flames_header_genecounts $root]\n
			return
		} elseif {$preset eq "sc"} {
			exec flames bulk_long_pipeline.py {*}$extraopts \
				--genomefa [file tail $refseq] \
				--gff3 [file tail $reftranscripts] \
				--outdir tempflames \
				--reference_csv $referencecsv \
				--config_file config.json \
				--fq_dir fastq \
				>@ stdout 2>@ stderr
			foreach file [glob tempflames/*] {
				file rename -force $file .
			}
			file delete tempflames
		} else {
			exec flames bulk_long_pipeline.py {*}$extraopts \
				--genomefa [file tail $refseq] \
				--gff3 [file tail $reftranscripts] \
				--outdir tempflames \
				--config_file config.json \
				--fq_dir fastq \
				>@ stdout 2>@ stderr
			foreach file [glob tempflames/*] {
				file rename -force $file .
			}
			file delete tempflames
		}
		#
		cd $destdir
		cg gzip -o transcripts-$root.isoforms.gff3.gz $flamesdir/isoform_annotated.filtered.gff3
		cg gzip -o transcripts-$root.isoforms.fa.gz $flamesdir/transcript_assembly.fa
		# make isoform_counts
		# read counts
		set f [gzopen $flamesdir/transcript_count.csv.gz]
		set header [split [gets $f] ,]
		unset -nocomplain a
		while {[gets $f line] != -1} {
			set line [split $line ,]
			set a([lrange $line 0 1]) [lindex $line 2]
		}
		close $f
		# make genepred and add counts
		exec flames gff3ToGenePred $flamesdir/isoform_annotated.filtered.gff3 $flamesdir/isoform_annotated.filtered.genepred
		set o [open $flamesdir/isoform_counts-$root.tsv w]
		puts $o [iso_flames_header_isoformcounts $root]
		set f [open $flamesdir/isoform_annotated.filtered.genepred]
		unset -nocomplain genea ; unset -nocomplain transcriptsa
		foreach {gene transcript count} [lrange [exec cg gtf2tsv $reftranscripts | cg select -g {gene * transcript *}] 3 end] {
			set genea($gene) 1
			set transcripta($transcript) 1
		}
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set transcript [lindex $line 0]
			regsub ^transcript: $transcript {} transcript
			lset line 0 $transcript
			set gene [lindex $line 11]
			regsub ^gene: $gene {} gene
			lset line 11 $gene
			lset line 8 [string trimright [lindex $line 8] ,]
			lset line 9 [string trimright [lindex $line 9] ,]
			if {[info exists a([list $transcript $gene])]} {
				set count $a([list $transcript $gene])
			} else {
				set count 0
			}
			if {![info exists genea($gene)]} {
				set category intergenic
			} elseif {![info exists transcripta($transcript)]} {
				set category novel
			} else {
				set category known
			}
			puts $o [join $line \t]\ttranscript\t$category\t$count
		}
		close $f
		close $o
		cg select -s - $flamesdir/isoform_counts-$root.tsv $flamesdir/isoform_counts-$root.stsv
		file rename -force $flamesdir/isoform_counts-$root.stsv isoform_counts-$root.tsv
		cg_iso_flames_genecounts isoform_counts-$root.tsv gene_counts-$root.tsv.temp $root
		file rename -force gene_counts-$root.tsv.temp gene_counts-$root.tsv
		file delete -force $flamesdir
	}
}

proc cg_iso_flames {args} {
	set args [job_init {*}$args]
	iso_flames_job {*}$args
	job_wait
}
