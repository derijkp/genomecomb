proc cg_tsv2sam {args} {
	set tsvfile -
	set samfile -
	set namecol name
	set outformat sam
	set refseq {}
	set threads 1
	cg_options tsv2sam args {
		-outformat {
			if {$value ni "sam bam cram"} {error "tsv2sam error: unknown format $value for -outformat"}
			set outformat $value
		}
		-refseq {
			set refseq [refseq $value]
		}
		-threads {set threads $value}
	} {tsvfile samfile} 0
	if {[file extension $samfile] eq ".bam"} {
		set outformat bam
	} elseif {[file extension $samfile] eq ".sam"} {
		set outformat sam
	} elseif {[file extension $samfile] eq ".cram"} {
		set outformat cram
	}
	set tsvcompressed [gziscompressed $tsvfile]
	if {$samfile eq "-"} {
		set pipe {}
	} elseif {$tsvcompressed} {
		set pipe [list {*}[gzcat $tsvfile] $tsvfile]
	}
	lappend pipe tsv2sam
	if {$outformat eq "bam"} {
		lappend pipe \| samtools view --threads $threads --no-PG -h -b --no-PG
	} elseif {$outformat eq "cram"} {
		if {$refseq eq ""} {error "tsv2sam error: outformat cram requires a reference sequence (user -refseq option)"}
		lappend pipe \| samtools view --threads $threads --no-PG -h -C -T $refseq --no-PG
	}
	if {$samfile eq "-"} {
		lappend pipe >@ stdout
	} else {
		lappend pipe > $samfile
	}
	if {$tsvfile eq "-"} {
		lappend pipe <@ stdin
	} elseif {!$tsvcompressed} {
		lappend pipe < $tsvfile
	}
	catch_exec {*}$pipe
}

proc cg_tsv2bam {args} {
	cg_tsv2sam -outformat bam {*}$args
}

