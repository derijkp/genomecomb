#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc vcfheader2table {lines {sheader {}} {prekeys {}} {extra {}}} {
	set pos 0
	unset -nocomplain a
	foreach prekey $prekeys field $sheader {
		set a($prekey) $pos
		incr pos
	}
	unset -nocomplain donea
	set data {}
	foreach line $lines {
		set line [string trim $line <>]
		if {$extra ne ""} {append line ,$extra}
		set resultline [list_fill [llength $sheader] {}]
		set pos 0
		while 1 {
			set npos [string first = $line $pos]
			if {$npos == -1} break
			set key [string range $line $pos [expr {$npos-1}]]
			set pos [expr {$npos+1}]
			set c [string index $line $pos]
			if {$c == "\""} {
				incr pos
				set tpos $pos
				while 1 {
					set npos [string first \" $line $tpos]
					if {[string index $line [expr {$npos+1}]] != "\""} break
					set tpos [expr {$npos+2}]
				}
				incr npos
			} else {
				set npos [string first , $line $pos]
				if {$npos == -1} {set npos [expr {[string length $line]}]}
			}
			set value [string range $line $pos [expr {$npos-1}]]
			set pos [expr {$npos+1}]
			set value [string trim $value \"]
			if {$key eq "ID"} {
				if {[info exists donea($value)]} {
					set resultline {}
					break
				} else {
					set donea($value) 1
				}
			}
			if {![info exists a($key)]} {
				set a($key) [llength $sheader]
				lappend sheader $key
				lappend prekeys $key
				set newdata {}
				foreach l $data {
					lappend l {}
					lappend newdata $l
				}
				set data $newdata
				lappend resultline $value
			} else {
				lset resultline $a($key) $value
			}
		}
		if {![llength $resultline]} continue
		lappend data $resultline
	}
	list $sheader $data $prekeys
}

proc vcf2tsvheader_line {sheaderlen table line} {
	if {[llength $line] < $sheaderlen} {
		lappend line {*}[list_fill [expr {$sheaderlen - [llength $line]}] {}]
	}
	return "#$table\t[join $line \t]"
}

proc vcf2tsvheader {vcfheader header split meta typelist {nheaderVar {}}} {
	if {$nheaderVar ne ""} {upvar $nheaderVar nheader}
	array set conv_formata {
		AD alleledepth
		GT genotype
		DP coverage
		FT gfilter
		GL loglikelihood
		GQ genoqual
		HQ haploqual
		AN totalallelecount
		PS phaseset

		AC allelecount
		AF frequency
		AA Ancestralallele
		DB dbsnp
		H2 Hapmap2
	}
	array set number_conv [lrange $typelist 1 end]
	# read data in array a
	unset -nocomplain a
	unset -nocomplain donea
	foreach line $vcfheader {
		if {[regexp {##([^=]+)=(.*)$} $line temp key value]} {
			lappend a($key)	$value
		} else {
			lappend a(comments)	$value
		}
	}
	set result {}
	# start printing
	lappend result "#filetype\ttsv/varfile"
	lappend result "#fileversion\t[version fileformat]"
	lappend result "#split\t$split"
	lappend result "#info\ttsv converted from vcf"
	foreach {key value} $meta {
		lappend result "#$key\t$value"
	}
	# samples
	set samples [lrange $header 9 end]
	if {[llength $samples] == 1} {
		lappend result "#numsamples\t1"
		lappend result "#samplename\t[lindex $samples 0]"
	} else {
		lappend result "#numsamples\t[llength $samples]"
	}
	#
	set nheader {chromosome begin end type ref alt name quality filter}
	set formatfields {GT}
	set headerfields {alleleSeq1 alleleSeq2 zyg phased genotypes}
	# FORMAT fields
	set lines [get a(FORMAT) ""]
	foreach {fheader fdata prekeys} [vcfheader2table $lines {field number type description source} {ID Number Type Description source} source=format] break
	unset -nocomplain a(FORMAT)
	# INFO fields
	set lines [get a(INFO) ""]
	foreach {sheader idata} [vcfheader2table $lines $fheader $prekeys source=info] break
	unset -nocomplain a(INFO)
	# default fields
	set sheaderlen [llength $sheader]
	lappend result "#fields\ttable"
	lappend result "#fields\t[join $sheader \t]"
	lappend result [vcf2tsvheader_line $sheaderlen fields {chromosome 1 String Chromosome/Contig var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {begin 1 Integer {Begin of feature (0 based - half open)} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {end 1 Integer {End of feature (0 based - half open)} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {type 1 String {Type of feature (snp,del,ins,...)} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {ref 1 String {Reference sequence, can be a number for large features} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {alt 1 String {Alternative sequence, can be a number for large features} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {name 1 String {name of feature} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {quality 1 Float {Quality score of feature} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {filter 1 String {Filter value} var}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {alleleSeq1 1 String {allele present on first chromosome/haplotype} geno}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {alleleSeq2 1 String {allele present on second chromosome/haplotype} geno}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {sequenced 1 String {sequenced status: v = variant, r = reference (i.e. not this variant), u = unsequenced} geno}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {zyg 1 String {Zygosity status: m = homozygous, t = heterozygous, r = reference, o = other variant, v = variant but genotype unspecified, c = compound (i.e. genotype has this variant and other variant), u = unsequenced} geno}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {phased 1 Integer {Phased status: 0 if not phased, other integer if phased} geno}]
	lappend result [vcf2tsvheader_line $sheaderlen fields {genotypes H Integer {Genotypes} geno}]
	# print FORMAT fields
	foreach line $fdata {
		foreach {field number} $line break
		if {[inlist {GT} $field]} continue
		if {[info exists number_conv($field)]} {
			set number $number_conv($field)
			lset line 1 $number
		}
		if {[info exists conv_formata($field)]} {
			set field $conv_formata($field)
			lset line 0 $field
		}
		if {$number eq "R"} {
			set temp $line
			lset temp 0 ${field}_ref
			lset temp 1 1
			lset temp 3 "reference only value of: [lindex $temp 3]"
			lset line 3 "alleles only values of: [lindex $line 3]"
			lappend result [vcf2tsvheader_line $sheaderlen fields $temp]
			lappend headerfields ${field}_ref
			lset line 1 A
		}
		lappend result [vcf2tsvheader_line $sheaderlen fields $line]
		set donea($field) 1
		lappend headerfields $field
	}
	# sample fields
	if {[llength $samples] == 1} {
		lappend nheader {*}$headerfields
	} else {
		foreach sample $samples {
			foreach field $headerfields {
				lappend nheader $field-$sample
			}
		}
	}
	# print info fields
	foreach line $idata {
		foreach {field number} $line break
		if {[inlist {END SVLEN SVTYPE} $field]} continue
		if {[info exists number_conv($field)]} {
			set number $number_conv($field)
			lset line 1 $number
		}
		if {$field eq "DP"} {
			set field totalcoverage
		} elseif {[info exists conv_formata($field)]} {
			set field $conv_formata($field)
		}
		if {$number eq "R"} {
			set temp $line
			if {[info exists donea(${field}_ref)]} {
				lset temp 0 info_${field}_ref
			} else {
				lset temp 0 ${field}_ref
			}
			lset temp 1 1
			lset temp 3 "reference only value of: [lindex $temp 3]"
			lset line 3 "alleles only values of: [lindex $line 3]"
			lappend result [vcf2tsvheader_line $sheaderlen fields $temp]
			lappend nheader ${field}_ref
			lset line 1 A
		}
		if {[info exists donea($field)]} {set field info_$field}
		lset line 0 $field
		lappend result [vcf2tsvheader_line $sheaderlen fields $line]
		lappend nheader $field
	}
	# other data in vcf header
	set fields [array names a]
	foreach field $fields {
		set lines $a($field)
		if {[string index [lindex $lines 0] 0] ne "<"} {
			foreach line $a($field) {
				lappend result \#vcf_$field\t$line
			}
			continue
		}
		foreach {sheader data} [vcfheader2table $lines] break
		lappend result \#$field\ttable
		lappend result \#$field\t[join $sheader \t]
		foreach line $data {
			lappend result "\#$field\t[join $line \t]"
		}
	}
	return $result
}

proc cg_vcfheader2tsv {args} {
	array set changefieldsa {fileformat vcf_fileformat}
	set typelist ". AD R RPA R AC A AF A"
	set showheader 1
	set split 1
	set meta {}
	cg_options vcfheader2tsv args {
		-changefields {
			array set changefieldsa $value
		}
		-showheader {
			set showheader 0
		}
		-split {
			set split $value
		}
		-meta {
			set meta $value
		}
		-typelist {
			set typelist $value
		}
	} {infile outfile} 0 2
	if {[info exists infile]} {
		set f [gzopen $infile]
	} else {
		set f stdin
	}
	if {[info exists outfile]} {
		set o [open $outfile w]
	} else {
		set o stdout
	}
	#
	set vcfheader {}
	while 1 {
		if {[gets $f line] == -1} break
		if {[string range $line 0 1] ne "##"} break
		lappend vcfheader $line
	}
	if {[string index $line 0] eq "#"} {
		set header [string range $line 1 end]
	} else {
		error "line with fieldnames missing from vcf"
	}
	gzclose $f
	puts $o [join [vcf2tsvheader $vcfheader $header $split $meta $typelist nheader] \n]
	if {$showheader} {
		puts $o [join $nheader \t]
	}
}
