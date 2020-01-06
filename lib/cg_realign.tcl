proc methods_realign {args} {
	set cmd realign
	set supportedmethods [list gatk]
	set start [expr {[string length cg_${cmd}] + 1}]
	foreach name [array names ::auto_index cg_${cmd}_*] {
		list_addnew supportedmethods [string range $name $start end]
	}
	if {[llength $args]} {
		set value [lindex $args 0]
		if {$value eq "1"} {set value [lindex $supportedmethods 0]}
		if {$value ni $supportedmethods} {error "$cmd: unsupported -method $value"}
		return $value
	} else {
		return $supportedmethods
	}
}

proc cg_realign {args} {
	set method gatk
	set regionfile {}
	set threads 2
	set refseq {}
	set bamfile -
	set sourcefile -
	set resultfile -
	set inputformat -
	set outputformat -
	cg_options realign_gatk args {
		-method {
			set method [methods_realign $value]
		}
		-regionfile {
			set regionfile $value
		}
		-refseq {
			set refseq $value
		}
		-inputformat - -if {
			set inputformat $value
		}
		-outputformat - -of {
			set outputformat $value
		}
		-threads - -t {
			set threads $value
		}
	} {sourcefile resultfile refseq} 0 3 {
		realign around indels using gatk
	}
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam cram sam}]}
	set inputformat [gzroot $inputformat]
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam cram sam}]}
	set refseq [refseq $refseq]
	cg_realign_${method} -regionfile $regionfile -refseq $refseq \
		-inputformat $inputformat -outputformat $outputformat -threads $threads \
		$sourcefile $resultfile
}
