#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc vcfheader2table {lines {sheader {}} {prekeys {}}} {
	set pos 0
	unset -nocomplain a
	foreach prekey $prekeys field $sheader {
		set a($prekey) $pos
		incr pos
	}
	set data {}
	foreach line $lines {
		set line [string trim $line <>]
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
				
			} else {
				set npos [string first , $line $pos]
				if {$npos == -1} {set npos [expr {[string length $line]}]}
			}
			set value [string range $line $pos [expr {$npos-1}]]
			set pos [expr {$npos+1}]
			set value [string trim $value \"]
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
		lappend data $resultline
	}
	list $sheader $data $prekeys
}

proc putsvcf2tsvheader {o sheaderlen table line} {
	if {[llength $line] < $sheaderlen} {
		lappend line {*}[list_fill [expr {$sheaderlen - [llength $line]}] {}]
	}
	puts $o "#$table\t[join $line \t]"
}

proc cg_vcfheader2tsv {args} {
	array set changefieldsa {fileformat vcf_fileformat}
	set showheader 1
	set split 1
	cg_options vcf2tsv args {
		-changefields {
			array set changefieldsa $value
		}
		-showheader {
			set showheader 0
		}
		-split {
			set split $value
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

	array set conv_formata {
		AD alleledepth
		GT genotype
		DP coverage
		FT filter
		GL loglikelihood
		GQ genoqual
		HQ haploqual
		AN totalallelecount

		AC allelecount
		AF frequency
		AA Ancestralallele
		DB dbsnp
		H2 Hapmap2
	}
	# read data in array a
	unset -nocomplain a
	unset -nocomplain donea
	while 1 {
		if {[gets $f line] == -1} break
		if {[string range $line 0 1] ne "##"} break
		if {[regexp {##([^=]+)=(.*)$} $line temp key value]} {
			lappend a($key)	$value
		} else {
			lappend a(comments)	$value
		}
	}
	if {[string index $line 0] eq "#"} {
		set header [string range $line 1 end]
	}
	close $f
	# start printing
	puts $o "#filetype\ttsv/varfile"
	puts $o "#fileversion\t[version fileformat]"
	puts $o "#split\t$split"
	puts $o "#info\ttsv converted from vcf"
	# samples
	set samples [lrange $header 9 end]
	if {[llength $samples] == 1} {
		puts $o "#numsamples\t1"
		puts $o "#samplename\t[lindex $samples 0]"
	} else {
		puts $o "#numsamples\t[llength $samples]"
	}
	#
	set nheader {chromosome begin end type ref alt name quality filter}
	set formatfields {GT}
	set headerfields {alleleSeq1 alleleSeq2 zyg phased genotypes}
	# FORMAT fields
	set lines [get a(FORMAT) ""]
	foreach {fheader fdata prekeys} [vcfheader2table $lines {field number type description} {ID Number Type Description}] break
	unset -nocomplain a(FORMAT)
	# INFO fields
	set lines [get a(INFO) ""]
	foreach {sheader idata} [vcfheader2table $lines $fheader $prekeys] break
	unset -nocomplain a(INFO)
	# default fields
	set sheaderlen [llength $sheader]
	puts $o "#fields\ttable"
	puts $o "#fields\t[join $sheader \t]"
	putsvcf2tsvheader $o $sheaderlen fields {chromosome 1 String Chromosome/Contig}
	putsvcf2tsvheader $o $sheaderlen fields {begin 1 Integer {Begin of feature (0 based - half open)}}
	putsvcf2tsvheader $o $sheaderlen fields {end 1 Integer {End of feature (0 based - half open)}}
	putsvcf2tsvheader $o $sheaderlen fields {type 1 String {Type of feature (snp,del,ins,...)}}
	putsvcf2tsvheader $o $sheaderlen fields {ref 1 String {Reference sequence, can be a number for large features}}
	putsvcf2tsvheader $o $sheaderlen fields {alt 1 String {Alternative sequence, can be a number for large features}}
	putsvcf2tsvheader $o $sheaderlen fields {name 1 String {name of feature}}
	putsvcf2tsvheader $o $sheaderlen fields {quality 1 Float {Quality score of feature}}
	putsvcf2tsvheader $o $sheaderlen fields {filter 1 String {Filter value}}
	putsvcf2tsvheader $o $sheaderlen fields {alleleSeq1 1 String {allele present on first chromosome/haplotype}}
	putsvcf2tsvheader $o $sheaderlen fields {alleleSeq2 1 String {allele present on second chromosome/haplotype}}
	putsvcf2tsvheader $o $sheaderlen fields {sequenced 1 String {sequenced status: v = variant, r = reference (i.e. not this variant), u = unsequenced}}
	putsvcf2tsvheader $o $sheaderlen fields {zyg 1 String {Zygosity status: m = homozygous, t = heterozygous, r = reference, o = other variant, c = compound, i.e. genotype has this variant and other variant}}
	putsvcf2tsvheader $o $sheaderlen fields {phased 1 Integer {Phased status: 0 if not phased, other integer if phased (same as variants in phase)}}
	putsvcf2tsvheader $o $sheaderlen fields {genotypes H Integer {Genotypes}}
	# print FORMAT fields
	foreach line $fdata {
		foreach {field number} $line break
		if {[inlist {GT} $field]} continue
		if {[info exists conv_formata($field)]} {
			set field $conv_formata($field)
			lset line 0 $field
		}
		if {$number eq "R"} {
			set temp $line
			lset temp 0 ${field}_ref
			putsvcf2tsvheader $o $sheaderlen fields $temp
			lappend headerfields ${field}_ref
		}
		putsvcf2tsvheader $o $sheaderlen fields $line
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
		set field [lindex $line 0]
		if {[inlist {END SVLEN SVTYPE} $field]} continue
		if {$field eq "DP"} {
			set field totalcoverage
		} elseif {[info exists conv_formata($field)]} {
			set field $conv_formata($field)
		}
		if {[info exists donea($field)]} {set field info_$field}
		lset line 0 $field
		putsvcf2tsvheader $o $sheaderlen fields $line
		lappend nheader $field
	}
	# other data in vcf header
	set fields [array names a]
	foreach field $fields {
		set lines $a($field)
		if {[string index [lindex $lines 0] 0] ne "<"} {
			foreach line $a($field) {
				puts $o \#vcf_$field\t$line
			}
			continue
		}
		foreach {sheader data} [vcfheader2table $lines] break
		puts $o \#$field\ttable
		puts $o \#$field\t[join $sheader \t]
		foreach line $data {
			puts $o "\#$field\t[join $line \t]"
		}
	}
	if {$showheader} {
		puts $o [join $nheader \t]
	}
}
