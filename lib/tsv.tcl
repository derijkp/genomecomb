#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_empty {file} {
	set f [open $file]
	tsv_open $f
	set read [gets $f line]
	close $f
	if {$read == -1} {return 1} else {return 0}
}

proc tsv_open {f {commentVar {}} {lineVar {}}} {
	if {$commentVar ne ""} {
		upvar $commentVar comment
	}
	if {$lineVar ne ""} {
		upvar $lineVar line
	}
	set comment {}
	set keep 0
	set buffering [fconfigure $f -buffering]
	fconfigure $f -buffering line
	set split 1
	if {![info exists line]} {
		if {[gets $f line] == -1} {
			return {}
		}
	}
	if {[regexp {^@HD[\t]VN} $line]} {
		# sam file
		while {![eof $f]} {
			set fchar [string index $line 0]
			if {$fchar ne "@"} {
				break
			}
			lappend comment \#$line
			set pos [tell $f]
			set line [gets $f]
		}
		set comment [join $comment \n]\n
		seek $f $pos start
		set len [llength [split $line \t]]
		set header {qname flag rname pos mapq cigar rnext pnext tlen seq qual}
		set size [expr {$len-11}]
		for {set i 1} {$i <= $size} {incr i} {
			lappend header opt$i
		}
		return $header
	}
	set fchar [string index $line 0]
	set fchar2 [string index $line 1]
	if {[regexp {##fileformat=VCF} $line]} {
		set vcf 1
	} else {
		set vcf 0
	}
	set fpos 0
	while {![eof $f]} {
		while {![string length $line]} {
			lappend comment \#
			set fpos [tell $f]
			if {[gets $f line] == -1} break
		}
		set fchar [string index $line 0]
		if {$fchar eq ">"} {
			break
		} elseif {$fchar ne "#"} {
			break
		} elseif {$vcf} {
			set fchar2 [string index $line 1]
			if {$fchar2 ne "#"} {
				set split 0
				break
			}
		}
		lappend comment $line
		set fpos [tell $f]
		set line [gets $f]
	}
	set ::keepfpos($f) $fpos
	fconfigure $f -buffering $buffering
	set fchar [string index $line 0]
	if {[inlist {# >} $fchar]} {
		set comment [join $comment \n]\n
		if {$vcf} {append comment \#}
		if {!$split} {
			return [string range $line 1 end]
		} else {
			return [split [string range $line 1 end] \t]
		}
	} else {
		if {[llength $comment]} {
			set comment [join $comment \n]\n
		}
		return [split $line \t]
	}
}

proc tsv_next {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [getnotempty $f]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	if {![eof $f]} {
		return $fpos
	} else {
		return $x
	}
}

proc tsv_nextline {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [getnotempty $f]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	return $line
}

proc tsv_makeregular {source resultfile {extra *}} {
	set fields [cg select -h $source]
	set poss [tsv_basicfields $fields 14 0]
	set oldfields [list_sub $fields $poss]
	set makefields {}
	foreach newfield {chromosome begin end type ref alt strand exonStarts exonEnds cdsStart cdsEnd transcript gene_name geneid} \
		oldfield $oldfields pos $poss {
		if {$pos == -1} continue
		if {$newfield eq $oldfield} {
			lappend makefields $oldfield
		} else {
			lappend makefields $newfield=\$$oldfield
		}
	}
	lappend makefields {*}$extra
	cg select -overwrite 1 -f $makefields $source $resultfile.temp[gzext $resultfile]
	file rename -force $resultfile.temp[gzext $resultfile] $resultfile
}

proc tsv_basicfields {header {num 6} {giveerror 1}} {
	set poss [list_cor $header {chromosome begin end type ref alt strand exonStarts exonEnds cdsStart cdsEnd transcript gene_name geneid}]
	set nfposs [list_find $poss -1]
	foreach nfpos $nfposs {
		switch $nfpos {
			0 {
				foreach name {chrom chr chr1 genoName tName contig} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			1 {
				foreach name {start end1 chromStart genoStart tStart txStart pos offset} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			2 {
				foreach name {start2 chromEnd genoEnd tEnd txEnd} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			4 {
				set v [lsearch $header reference]
			}
			5 {
				set v [lsearch $header alternative]
			}
			11 {
				foreach name {transcript_id transcriptid name} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			12 {
				foreach name {gene geneid gene_id name2} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			13 {
				foreach name {gene_id gene_name gene name2} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			default {
				continue
			}
		}
		lset poss $nfpos $v
	}
	incr num -1
	set poss [lrange $poss 0 $num]
	set pos [lsearch $poss -1]
	if {$giveerror ne "0" && ($pos != -1)} {
		set notfound [list_sub {chromosome begin end type ref alt} [list_find $poss -1]]
		if {$giveerror eq "1"} {
			error "header error: fields (or alternatives) not found: $notfound"
		} else {
			error "header error in $giveerror: fields (or alternatives) not found: $notfound"
		}
	}
	return $poss
}

array set altfieldsa {
	chromosome {chrom chr chr1 genoName tName contig}
	begin {start end1 chromStart genoStart tStart txStart pos offset}
	end {start2 chromEnd genoEnd tEnd txEnd}
	ref {reference}
	alt {alternative}
	transcript {transcript_id transcriptid name}
	gene_name {gene geneid gene_id name2}
	gene {gene_name geneid gene_id name2}
	geneid {gene_id gene_name gene name2}
}

proc findfields {header fields} {
	global altfieldsa
	set result {}
	foreach field $fields {
		if {$field in $header} {
			lappend result $field
		} elseif {[info exists altfieldsa($field)]} {
			set found {}
			foreach name $altfieldsa($field) {
				set v [lsearch $header $name]
				if {$v != -1} {
					set found $name
					break
				}
			}
			lappend result $found
			
		} else {
			lappend result {}
		}
	}
	return $result
}

proc findfieldspos {header fields} {
	global altfieldsa
	set result {}
	foreach field $fields {
		set pos [lsearch $header $field]
		if {$pos != -1} {
			lappend result $pos
		} elseif {[info exists altfieldsa($field)]} {
			set found -1
			foreach name $altfieldsa($field) {
				set v [lsearch $header $name]
				if {$v != -1} {
					set found $v
					break
				}
			}
			lappend result $found
		} else {
			lappend result -1
		}
	}
	return $result
}

proc tsv_comment2var {comment varName} {
	upvar $varName a
	unset -nocomplain a
	set keepcomments {}
	foreach line [split $comment \n] {
		if {$line eq ""} continue
		set line [split $line \t]
		if {[lrange $line 0 8] eq "{#CHROM} POS ID REF ALT QUAL FILTER INFO FORMAT"} continue
		if {[llength $line] == 1} {
			lappend a() [lindex $line 0]
			continue
		}
		set key [string range [lindex $line 0] 1 end]
		set value [lrange $line 1 end]
		lappend a($key) $value
	}
}

proc tsv_var2comment {varName} {
	upvar $varName a
	set keys [list_remove [array names a] {}]
	set deffields {
		filetype fileversion split info refseq numsamples samplename
		fields contig
	}
	set keys [list_concat [list_common $deffields $keys] [list_lremove $keys $deffields]]
	set result {}
	foreach key $keys {
		foreach line $a($key) {
			lappend result \#$key\t[join $line \t]
		}
	}
	return [join $result \n]
}

proc tsv_count {tsvfile} {
	set countfile [indexdir_file $tsvfile vars.tsv.count ok]
	if {!$ok} {
		set varsfile [tsv_varsfile $tsvfile]
		if {$ok} {
			set count [lindex [cg select -g - $varsfile] end]
			file_write $countfile $count
			return $count
		} else {
			set count [lindex [cg select -g - [gzfile $tsvfile]] end]
			file_write $countfile $count
			return $count
		}
	}
	return [file_read $countfile]
}

# creates or updates (if needed) and returns a index file containing only the vars columns of $tsvfile
# if varsfile is given, the index file will be created as $varsfile instead of at the usual location
proc tsv_varsfile {tsvfile {varsfile {}}} {
	set tsvfile [gzfile $tsvfile]
	if {$varsfile eq ""} {
		set varsfile [indexdir_file $tsvfile vars.tsv ok]
	} else {
		set ok 0
	}
	if {!$ok} {
		set f [gzopen [gzfile $tsvfile]]
		set header [tsv_open $f comment]
		file_write $varsfile.temp $comment
		gzclose $f
		exec varsfile $tsvfile >> $varsfile.temp
		file rename -force -- $varsfile.temp $varsfile
	}
	return $varsfile
}

# if needed, create or update vars.tsv file in index to use in annotation (to avoid using the larger orifile)
proc tsv_varsfile_job {orifile {skips {}} {usefile {}}} {
	upvar job_logdir job_logdir
	if {$usefile eq ""} {
		set usefile [indexdir_file $orifile vars.tsv]
	}
	if {$usefile eq $orifile} {
		error "internal error: Trying to create smaller varfile as cache, but target is same as original"
	}
	job [job_relfile2name annot-createusefile- $usefile] {*}$skips -deps {
		$orifile
	} -targets {
		$usefile
	} -vars {
		ok orifile usefile dbfiles
	} -code {
		# usefile: smaller file with only variants used for actual annotation; 
		# if orifile is small, a link to it is made.
		# If it contains too many extra columns a cut down version is made
		set f [gzopen $orifile]
		set header [tsv_open $f]
		catch {gzclose $f}
		if {[gziscompressed $orifile] || [file dir $target] ne "$orifile.index" || ([llength $header] > 10 && [llength $dbfiles] >= 4)} {
			tsv_varsfile $orifile $usefile
		} else {
			mklink $orifile $target
		}
	}
	return $usefile
}

proc tsv_convert2var {file headerVar {commentVar {}}} {
	upvar $headerVar header
	if {$commentVar ne ""} {
		upvar $commentVar comment
	}
	set comment ""
	set f [gzopen $file]
	set line [gets $f]
	if {[regexp {##fileformat=VCF} $line]} {
		catch {gzclose $f}
		set tempfile [tempfile]
		exec vcf2tsv 1 [gztemp $file] | cg select -s - > $tempfile
		set f [open $tempfile]
		set header [tsv_open $f comment]
		close $f
	} else {
		set tempfile $file
		set header [tsv_open $f comment line]
		catch {gzclose $f}
	}
	return $tempfile
}

proc header {file} {
	set f [gzopen $file]
	set header [tsv_open $f]
	if {[catch {gzclose $f} msg]} {
		error "error reading header of file $file: $msg"
	}
	return $header
}
