proc version_flames {} {
	set bin [follow_links [exec which flames]]
	set flamesdir [file dir $bin]
	lindex [split [file tail $flamesdir] -] end
}

proc cg_iso_flames_genecounts {isoformcounts genecounts} {
	set tempfile [tempfile]
	cg select -overwrite 1 \
		-g {geneid} \
		-gc {distinct(chromosome),min(begin),max(end),distinct(strand),ucount(transcript),sum(counts-*)} \
		$isoformcounts $tempfile
	catch {close $f} ; catch {close $o}
	set f [open $tempfile]
	set header [tsv_open $f]
	set oheader {gene chromosome begin end strand nrtranscripts}
	foreach field [lrange $header 6 end] {
		regsub sum_ $field {} field
		lappend oheader $field
	}
	set o [open $genecounts w]
	puts $o [deindent {
		#filetype	tsv/genecountsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	gene	1	String	gene id (usually gene_id from GTF)
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	nrtranscripts	1	Integer	number of transcripts found
		#fields	counts	1	Integer	Number of reads mapping to gene
	}]
	puts $o [join $oheader \t]
	while {[gets $f line] != -1} {
		puts $o $line
	}
	close $o
	close $f
}

proc iso_flames_job {args} {
	# putslog [list flames_job {*}$args]
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg flames {*}$args]
	global appdir
	set refseq {}
	set skips {}
	set genes {}
	set extraopts {}
	set sqanti 1
	set compar joint
	set threads 8
	set resultfile {}
	set reftranscripts {}
	cg_options flames args {
		-refseq {
			set refseq $value
		}
		-compar {
			set compar $value
		}
		-threads {
			set threads $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-extraopts {
			set extraopts $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {fastqdir resultfile} 1 2
	set fastqdir [file_absolute $fastqdir]
	if {![file isdir $fastqdir]} {
		set bam $fastqdir
		set fastqdir [file dir $bam]/fastq
		if {![file exists $fastqdir]} {
			error "could not find fastq dir for $bam ($fastqdir)"
		}
	}
	if {$resultfile eq ""} {
		set root [file_rootname [file dir $fastqdir]]
		set resultfile [file dir $fastqdir]/isoform_counts-flames-$root.tsv
	} else {
		set resultfile [file_absolute $resultfile]
		set root [file_rootname $resultfile]
	}
	set resultfile [file_absolute $resultfile]
	set reftranscripts [file_absolute $reftranscripts]
	set resulttail [file tail $resultfile]
	set destdir [file dir $resultfile]
	#
	set refseq [refseq $refseq]
	if {$reftranscripts eq ""} {
		set reftranscripts [ref_gtftranscripts $refseq]
	}
	# 
	job_logfile $destdir/flames_[file tail $destdir] $destdir $cmdline \
		{*}[versions flames dbdir zstd os]
	# analysis per sample
	set flamesdir $destdir/flames-$root
	set fastqfiles [bsort [jobgzfiles $fastqdir/*.fq $fastqdir/*.fastq]]
	mkdir $flamesdir
	job flames-[file tail $resultfile] {*}$skips -skip flames-$root/counts_matrix-flames-$root.tsv \
	-cores $threads \
	-deps [list {*}$fastqfiles $refseq $reftranscripts] -targets {
		$resultfile
		$destdir/gene_counts-flames-$root.tsv
	} -vars {
		fastqdir fastqfiles resultfile reftranscripts root destdir flamesdir refseq reftranscripts threads extraopts
	} -code {
		analysisinfo_write [gzfile $fastqdir/*.fastq $fastqdir/*.fq] $resultfile flames [version flames]
		set keeppwd [pwd]
		cd $flamesdir
		file_write config.json [deindent {
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
			        "has_UMI":false
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
		}]
		mklink $refseq [file tail $refseq]
		mklink $reftranscripts [file tail $reftranscripts]
		mkdir fastq
		set mergedfastqfile fastq/$root.fastq
		if {[llength $fastqfiles] > 1} {
			exec cg zcat {*}$fastqfiles | gzip --fast > $mergedfastqfile.gz.temp
			file rename -force $mergedfastqfile.gz.temp $mergedfastqfile.gz
		} else {
			set fastqfile [lindex $fastqfiles 0]
			mklink -absolute 0 $fastqfile $mergedfastqfile[gzext $fastqfile]
		}
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
		#
		cd $destdir
		mklink -absolute 0 flames-$root/isoform_annotated.filtered.gff3 transcripts-flames-$root.isoforms.gff3
		mklink -absolute 0 flames-$root/transcript_assembly.fa transcripts-flames-$root.isoforms.fa
		# make isoform_counts
		# read counts
		set f [gzopen flames-$root/transcript_count.csv.gz]
		set header [split [gets $f] ,]
		unset -nocomplain a
		while {[gets $f line] != -1} {
			set line [split $line ,]
			set a([lrange $line 0 1]) [lindex $line 2]
		}
		close $f
		# make genepred and add counts
		exec flames gff3ToGenePred flames-$root/isoform_annotated.filtered.gff3 flames-$root/isoform_annotated.filtered.genepred
		set o [open flames-$root/isoform_counts-flames-$root.tsv w]
		puts $o [deindent {
			#filetype	tsv/transcriptsfile
			#fileversion	0.99
			#fields	table
			#fields	field	number	type	description
			#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
			#fields	chromosome	1	String	Chromosome name
			#fields	strand	1	String	+ or - for strand
			#fields	begin	1	Integer	Transcription start position
			#fields	end	1	Integer	Transcription end position
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
			#fields	counts	1	Integer	Number of reads mapping to isoform
		}]
		puts $o [join [list \
			transcript chromosome strand begin end cdsStart cdsEnd exonCount exonStarts exonEnds score geneid cdsStartStat cdsEndStat exonFrames counts-$root
		] \t]
		set f [open flames-$root/isoform_annotated.filtered.genepred]
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set transcript [lindex $line 0]
			regsub ^transcript: $transcript {} transcript
			lset line 0 $transcript
			set gene [lindex $line 11]
			regsub ^gene: $gene {} gene
			lset line 11 $gene
			if {[info exists a([list $transcript $gene])]} {
				set count $a([list $transcript $gene])
			} else {
				set count 0
			}
			puts $o [join $line \t]\t$count
		}
		close $f
		close $o
		file rename -force flames-$root/isoform_counts-flames-$root.tsv isoform_counts-flames-$root.tsv
		cg_iso_flames_genecounts isoform_counts-flames-$root.tsv gene_counts-flames-$root.tsv.temp
		file rename -force gene_counts-flames-$root.tsv.temp gene_counts-flames-$root.tsv
	}
}

proc cg_iso_flames {args} {
	set args [job_init {*}$args]
	iso_flames_job {*}$args
	job_wait
}
