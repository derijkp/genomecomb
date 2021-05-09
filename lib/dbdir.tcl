proc dbdir {args} {
	global dbdir
	if {[llength $args]} {
		set temp [lindex $args 0]
		if {$temp ne ""} {
			set dbdir $temp
		}
	}
	if {![info exists dbdir]} {
		if {[info exists ::env(GENOMECOMB_DBDIR)]} {
			set dbdir $::env(GENOMECOMB_DBDIR)
		} else {
			error "dbdir not specified, use options e.g. (-dbdir) or environment variable GENOMECOMB_DBDIR to set"
		}
	}
	if {[file isfile $dbdir] && [file extension $dbdir] eq ".ifas"} {
		set dbdir [file dir $dbdir]
	} elseif {![file isdir $dbdir]} {
		error "dbdir $dbdir does not exist, use options e.g. (-dbdir) or environment variable GENOMECOMB_DBDIR to change"
	}
	set dbdir [file_absolute $dbdir]
	# set REF_PATH for cram files
	set ::env(REFSEQ) [lindex [gzfiles $dbdir/genome_*.ifas] 0]
	set ::env(REF_PATH) $::env(REFSEQ).forcram
	return $dbdir
}

proc dbdir_ref dbdir {
	file tail $dbdir
}

proc refseq {{refseq {}} {dbdir {}}} {
	if {$refseq eq ""} {
		if {$dbdir ne ""} {
			set pattern $dbdir
		} elseif {[catch {
			set pattern [dbdir]
		} msg]} {
			error "refseq not defined, specify using the -refseq option. Alternatively you can also specify dbdir use options e.g. (-dbdir if supported) or environment variable GENOMECOMB_DBDIR to set"
		}
	} else {
		set pattern $refseq
	}
	if {[file isdir $pattern]} {
		set refseq [lindex [gzfiles $pattern/genome_*.ifas] 0]
	} else {
		set refseq [lindex [gzfiles $pattern $pattern/genome_*.ifas] 0]
	}
	if {![file exists $refseq]} {
		error "refseq not found ($pattern)"
	}
	# set REF_PATH for cram files
	set ::env(REFSEQ) $refseq
	set ::env(REF_PATH) $::env(REFSEQ).forcram
	if {![info exists ::dbdir]} {set ::dbdir [file dir $refseq]}
	return [file_absolute $refseq]
}

proc ref_chrsize {refseq chr} {
	global genomecomb_chrsizea
	if {![info exists genomecomb_chrsizea]} {
		list_foreach {tchr size} [split [string trim [cg select -sh /dev/null -hp {chromosome size} -f {chromosome size} $refseq.fai]] \n] {
			set genomecomb_chrsizea([chr_clip $tchr]) [list $tchr $size]
		}
	}
	set cchr [chr_clip $chr]
	if {[info exists genomecomb_chrsizea($cchr)]} {
		return [lindex $genomecomb_chrsizea($cchr) end]
	} else {
		return 536870912
	}
}

proc refcram {{dbdir {}}} {
	if {[info exists ::env(REF_PATH)]} {return $::env(REF_PATH)}
	set dbdir [dbdir $dbdir]
	set refseq [lindex [gzfiles $dbdir/genome_*.ifas] 0]
	set ::env(REF_PATH) $refseq.forcram
}
