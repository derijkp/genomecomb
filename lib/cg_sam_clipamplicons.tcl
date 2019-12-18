proc cg_sam_clipamplicons {args} {
	set refseq {}
	set sourcefile -
	set resultfile -
	set inputformat -
	set outputformat -
	cg_options sam_clipamplicons args {
		-refseq {
			set refseq $value
		}
		-inputformat - -if {
			set inputformat $value
		}
		-outputformat - -of {
			set outputformat $value
		}
	} {ampliconsfile sourcefile resultfile} 1 3 {
		realign around indels using gatk
	}
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile sam {bam cram sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile sam {bam cram sam}]}
	analysisinfo_write $sourcefile $resultfile clipamplicons genomecomb clipamplicons_version [version genomecomb] clipampliconsfile [file tail $ampliconsfile]
	set compressionlevel [defcompressionlevel 5]
	set f [open $ampliconsfile]
	set header [tsv_open $f]
	close $f
	set poss [tsv_basicfields $header 3]
	lappend poss [lsearch $header outer_begin] [lsearch $header outer_end]
	if {[inlist $poss -1]} {
		error "error in amplicons file: missing fields: [list_sub {chromosome begin end outer_begin outer_end} [list_find $poss -1]]"
	}
	set pipe {}
	set optsio {}
	if {$inputformat in "bam cram"} {
		if {$sourcefile eq "-"} {
			lappend optsio <@ stdin
			lappend pipe samtools view -h
		} else {
			lappend pipe samtools view -h $sourcefile
		}
		if {$inputformat eq "cram"} {
			lappend pipe -T [refseq $refseq]
		}
	} else {
		if {$sourcefile eq "-"} {
			lappend optsio <@ stdin
		} elseif {[gziscompressed $sourcefile]} {
			lappend pipe {*}[gzcat $file] $sourcefile
		} else {
			lappend optsio < $sourcefile
		}
	}
	if {[llength $pipe]} {lappend pipe |}
	lappend pipe sam_clipamplicons $ampliconsfile {*}$poss
	if {$outputformat eq "bam"} {
		lappend pipe | samtools view -hb --output-fmt-option level=$compressionlevel
	} elseif {$outputformat eq "cram"} {
		lappend pipe | samtools view -h -C -T [refseq $refseq]
	}
	if {[llength $optsio]} {lappend pipe {*}$optsio}
	if {$resultfile ne "-"} {
		lappend pipe > $resultfile
	} else {
		lappend pipe >@ stdout
	}
	catch_exec {*}$pipe
}
