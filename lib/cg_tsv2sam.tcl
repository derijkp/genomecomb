proc cg_tsv2sam {args} {
	set tsvfile -
	set samfile -
	set namecol name
	set outformat sam
	cg_options tsv2sam args {
		-outformat {set outformat $value}
	} {tsvfile samfile} 0
	if {[file extension $samfile] eq ".bam"} {
		set outformat bam
	} elseif {[file extension $samfile] eq ".sam"} {
		set outformat sam
	}
	set tsvcompressed [gziscompressed $tsvfile]
	if {$samfile eq "-"} {
		set pipe {}
	} elseif {$tsvcompressed} {
		set pipe [list {*}[gzcat $tsvfile] $tsvfile]
	}
	lappend pipe tsv2sam
	if {$outformat eq "bam"} {
		lappend pipe \| samtools view -h -b
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

