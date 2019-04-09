#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv2vcf_fielddata {field faVar} {
	upvar $faVar fa
	global tsv2vcf_prefa
	if {[info exists fa($field)]} {
		return $fa($field)
	}
	if {![info exists tsv2vcf_prefa]} {
		foreach line [split [deindent {
			chromosome	1	String	Chromosome/Contig	var
			begin	1	Integer	Begin of feature (0 based - half open)	var
			end	1	Integer	End of feature (0 based - half open)	var
			type	1	String	Type of feature (snp,del,ins,...)	var
			ref	1	String	Reference sequence, can be a number for large features	var
			alt	1	String	Alternative sequence, can be a number for large features	var
			name	1	String	name of feature	var
			quality	1	Float	Quality score of feature	var
			filter	1	String	Filter value	var
			alleleSeq1	1	String	allele present on first chromosome/haplotype	geno
			alleleSeq2	1	String	allele present on second chromosome/haplotype	geno
			sequenced	1	String	sequenced status: v = variant, r = reference (i.e. not this variant), u = unsequenced	geno
			zyg	1	String	Zygosity status: m = homozygous, t = heterozygous, r = reference, o = other variant, c = compound, i.e. genotype has this variant and other variant	geno
			phased	1	Integer	Phased status: 0 if not phased, other integer if phased (same as variants in phase)	geno
			alleledepth_ref	1	Integer	reference only value of: Allelic depths for the ref and alt alleles in the order listed	format
			alleledepth	A	Integer	alleles only values of: Allelic depths for the ref and alt alleles in the order listed	format
			coverage	1	Integer	Read Depth (filtered)	format
			genotypes	H	Integer	Genotypes	geno
			genoqual	1	Integer	Genotype Quality	format
			phaseset	1	Integer	Phase Set	format
			TE	A	Integer	test for alt alleles in the order listed	format
			haploqual	2	Integer	Haplotype Quality	format
			NS	1	Integer	Number of Samples With Data	info
			totalcoverage	1	Integer	Total Depth	info
			frequency	A	Float	Allele Frequency, for each ALT allele, in the same order as listed	info
			Ancestralallele	1	String	Ancestral Allele	info
			dbsnp	0	Flag	dbSNP membership, build 129	info
			Hapmap2	0	Flag	HapMap2 membership	info
			PGT	1	String	Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another	format
			PID	1	String	Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group	format
			PL	G	Integer	Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification	format
			RGQ	1	Integer	Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)	format
			SB	4	Integer	Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.	format
			allelecount	A	Integer	Allele count in genotypes, for each ALT allele, in the same order as listed	info
			totalallelecount	1	Integer	Total number of alleles in called genotypes	info
			AS_BaseQRankSum	A	Float	allele specific Z-score from Wilcoxon rank sum test of each Alt Vs. Ref base qualities	info
			AS_FS	A	Float	allele specific phred-scaled p-value using Fisher's exact test to detect strand bias of each alt allele	info
			AS_InbreedingCoeff	A	Float	allele specific heterozygosity as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation; relate to inbreeding coefficient	info
			AS_MQ	A	Float	Allele-specific RMS Mapping Quality	info
			AS_MQRankSum	A	Float	Allele-specific Mapping Quality Rank Sum	info
			AS_QD	A	Float	Allele-specific Variant Confidence/Quality by Depth	info
			AS_RAW_BaseQRankSum	1	String	raw data for allele specific rank sum test of base qualities	info
			AS_RAW_MQ	1	String	Allele-specfic raw data for RMS Mapping Quality	info
			AS_RAW_MQRankSum	1	String	Allele-specfic raw data for Mapping Quality Rank Sum	info
			AS_RAW_ReadPosRankSum	1	String	allele specific raw data for rank sum test of read position bias	info
			AS_ReadPosRankSum	A	Float	allele specific Z-score from Wilcoxon rank sum test of each Alt vs. Ref read position bias	info
			AS_SB_TABLE	1	String	Allele-specific forward/reverse read counts for strand bias tests	info
			AS_SOR	A	Float	Allele specific strand Odds Ratio of 2x|Alts| contingency table to detect allele specific strand bias	info
			BaseQRankSum	1	Float	Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities	info
			ClippingRankSum	1	Float	Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases	info
			DS	0	Flag	Were any of the samples downsampled?	info
			ExcessHet	1	Float	Phred-scaled p-value for exact test of excess heterozygosity	info
			FS	1	Float	Phred-scaled p-value using Fisher's exact test to detect strand bias	info
			InbreedingCoeff	1	Float	Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation	info
			MLEAC	A	Integer	Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed	info
			MLEAF	A	Float	Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed	info
			MQ	1	Float	RMS Mapping Quality	info
			MQRankSum	1	Float	Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities	info
			NDA	1	Integer	Number of alternate alleles discovered (but not necessarily genotyped) at this site	info
			QD	1	Float	Variant Confidence/Quality by Depth	info
			RAW_MQ	1	Float	Raw data for RMS Mapping Quality	info
			ReadPosRankSum	1	Float	Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias	info
			SOR	1	Float	Symmetric Odds Ratio of 2x2 contingency table to detect strand bias	info
			SP	1	Integer	Phred-scaled strand bias P-value	format
			INDEL	0	Flag	Indicates that the variant is an INDEL.	info
			IDV	1	Integer	Maximum number of reads supporting an indel	info
			IMF	1	Float	Maximum fraction of reads supporting an indel	info
			VDB	1	Float	Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)	info	3
			SGB	1	Float	Segregation based metric.	info
			RPB	1	Float	Mann-Whitney U test of Read Position Bias (bigger is better)	info
			MQB	1	Float	Mann-Whitney U test of Mapping Quality Bias (bigger is better)	info
			MQSB	1	Float	Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)	info
			BQB	1	Float	Mann-Whitney U test of Base Quality Bias (bigger is better)	info
			MQ0F	1	Float	Fraction of MQ0 reads (smaller is better)	info
			AF1	1	Float	Max-likelihood estimate of the first ALT allele frequency (assuming HWE)	info
			AF2	1	Float	Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)	info
			AC1	1	Float	Max-likelihood estimate of the first ALT allele count (no HWE assumption)	info
			DP4	4	Integer	Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases	info
			FQ	1	Float	Phred probability of all samples being the same	info
			PV4	4	Float	P-values for strand bias, baseQ bias, mapQ bias and tail distance bias	info
			G3	3	Float	ML estimate of genotype frequencies	info
			HWE	1	Float	Chi^2 based HWE test P-value based on G3	info
			GQX	1	Integer	Empirically calibrated genotype quality score for variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}	format
			DPF	1	Integer	Basecalls filtered from input prior to site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block.	format
			MIN_DP	1	Integer	Minimum filtered basecall depth used for site genotyping within a non-variant multi-site block	format
			ADF	.	Integer	Allelic depths on the forward strand	format
			ADR	.	Integer	Allelic depths on the reverse strand	format
			DPI	1	Integer	Read depth associated with indel, taken from the site preceding the indel	format
			gfilter	1	String	Sample filter, 'PASS' indicates that all filters have passed for this sample	format
			BLOCKAVG_min30p3a	0	Flag	Non-variant multi-site block. Non-variant blocks are defined independently for each sample. All sites in such a block are constrained to be non-variant, have the same filter value, and have sample values {GQX,DP,DPF} in range [x,y], y <= max(x+3,(x*1.3)).	info
			SNVHPOL	1	Integer	SNV contextual homopolymer length	info
			CIGAR	A	String	CIGAR alignment for each alternate indel allele	info
			RU	A	String	Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases	info
			REFREP	A	Integer	Number of times RU is repeated in reference	info
			IDREP	A	Integer	Number of times RU is repeated in indel allele	info
		}] \n] {
			set line [split $line \t]
			set tsv2vcf_prefa([lindex $line 0]) $line
		}
	}
	get tsv2vcf_prefa($field) [list $field . String	"$field, no further description found"]
}

proc tsv2vcf_getvar {line poss keyposaVar genomef} {
	upvar $keyposaVar keyposa
	foreach {chromosome begin end type ref alt} [list_sub $line $poss] break
	if {[info exists keyposa(name)]} {
		set id [lindex $line $keyposa(name)]
	} else {
		set id .
	}
	if {[info exists keyposa(quality)]} {
		set qual [lindex $line $keyposa(quality)]
	} else {
		set qual .
	}
	if {[info exists keyposa(filter)]} {
		set filter [lindex $line $keyposa(filter)]
	} else {
		set filter .
	}
	if {$type eq "snp"} {
		set pos $end
		set oref $ref
		set oalt $alt
		# ref and alt are ok as is in this case
	} elseif {$type eq "ins"} {
		if {[isint $alt]} {
			set oalt <$alt>
		} else {
			set oalt $alt
		}
		if {$begin != 0} {
			set pos $begin
			set oref [string toupper [genome_get $genomef $chromosome [expr {$begin - 1}] $end]]
			set oalt $oref$oalt
		} else {
			set pos 1
			set oref [string toupper [genome_get $genomef $chromosome 0 1]]
			set oalt $oalt$oref
		}
	} elseif {$type eq "del"} {
		set pos $begin
		if {$begin != 0} {
			set oref [string toupper [genome_get $genomef $chromosome [expr {$begin - 1}] $end]]
			set oalt [string index $oref 0]
		} else {
			set pos 1
			set oref [string toupper [genome_get $genomef $chromosome $begin [expr {$end + 1}]]]
			set oalt [string index $oref end]
		}
	} elseif {$type eq "sub"} {
		if {[isint $alt]} {
			set oalt <$alt>
		} else {
			set oalt $alt
		}
		set pos $begin
		set oref [string toupper [genome_get $genomef $chromosome [expr {$begin-1}] $end]]
		if {[isint $alt]} {
			set oalt [string index $oref 0]<$alt>
		} else {
			set oalt [string index $oref 0]$alt
		}
	}
	list $chromosome $pos $id $oref $oalt $qual $filter $ref $alt
}

proc tsv2vcf_outputheaderfield {o key aVar} {
	upvar $aVar a
	if {[info exists a($key)]} {
		set value $a($key)
		unset a($key)
	} elseif {[info exists a(vcf_$key)]} {
		set value $a(vcf_$key)
		unset a(vcf_$key)
	} else {
		return
	}
	if {[llength $value] == 1} {
		puts $o "##$key=[lindex $value 0]"
	} elseif {[lindex $value 0] eq "table"} {
		set tfields [lindex $value 1]
		foreach line [lrange $value 2 end] {
			set list {}
			foreach field $tfields v $line {
				if {[string first " " $v] != -1} {
					lappend list $field=\"$v\"
				} else {
					lappend list $field=$v
				}
			}
			puts $o "##$key=<[join $list ,]>"
		}
	} else {
		foreach v $value {
			puts $o "##$key=$v"
		}
	}
}

proc tsv2vcf_printlines {lines infofields infoposs infoflags infonumbers analyses formataVar} {
	#
	# position format fields
	upvar $formataVar formata
# putsvars lines infofields infoposs infoflags infonumbers analyses
# puts [list array set formata [array get formata]]
	# [llength $lines] == 1
	set resultline {}
	set vars [list_subindex $lines 0]
	set lines [list_subindex $lines 1]
	set refalts [list_subindex $vars 3 4]
	set alts [list_subindex $refalts 1]
	set rdup [list_remdup $refalts]
	set extra {}
	if {[llength $rdup] < [llength $refalts]} {
		set cor [list_cor $rdup $refalts]
		set poss [list_find $cor -1]
		foreach pos [list_reverse $poss] {
			set var [list_pop vars $pos]
			set line [list_pop lines $pos]
			lappend extra [tsv2vcf_printlines [list [list $var $line]] $infofields $infoposs $infoflags $infonumbers $analyses formata]
		}
	}
	foreach {chromosome chrpos id fref falt fqual ffilter fref falt} [lindex $vars 0] break
	#
	# ref and alt
	# -----------
	set types [list_subindex $vars 1]
	set refs [list_subindex $vars 3]
	set alts [list_subindex $vars 4]
	set linerefs [list_subindex $vars [expr {[llength [lindex $vars 0]] - 2}]]
	set linealts [list_subindex $vars [expr {[llength [lindex $vars 0]] - 1}]]
	set max 0
	set maxref {}
	foreach ref $refs {
		set len [string length $ref]
		if {$len > $max} {set max $len ; set maxref $ref}
	}
	set newalts {} ; set appends {}
	foreach ref $refs alt $alts {
		set len [string length $ref]
		if {$len < $max} {
			set append [string range $maxref end-[expr {$max-$len-1}] end]
			append alt $append
		} else {
			set append {}
		}
		lappend appends $append
		lappend newalts $alt
	}
	set vcfref $maxref
	set vcfalt [join $newalts ,]
	#
	# qual and filter
	# ---------------
	set quals [list_remdup [list_subindex $vars 5]]
	if {[llength $quals] > 1} {
		set qual [max [list_remove $quals {{} .}]]
	} else {
		set qual [lindex $quals 0]
	}
	set filters [list_remdup [list_subindex $vars 6]]
	if {[llength $filters] > 1} {
		if {[inlist $filters PASS]} {set filter PASS} else {set filter [lindex $filters 0]}
	} else {
		set filter [lindex $filters 0]
	}
	#
	# make info
	# ---------
	set info {}
	foreach field $infofields ipos $infoposs infoflag $infoflags infonumber $infonumbers {
		set value [list_subindex $lines $ipos]
		if {$infonumber in "A R" || ([isint $infonumber] && $infonumber > 1)} {
			if {[list_remove $value {}] eq ""} {
				set value {}
			} else {
				set value [join $value ,]
			}
		} else {
			set value [list_remove [list_remdup $value] {}]
		}
		if {$value eq ""} continue
		if {!$infoflag} {
			lappend info $field=[join $value ,]
		} elseif {$value == 1} {
			lappend info $field
		}
	}
	if {![llength $info]} {set info .}
	set info [join $info \;]

	#
	# make format data
	# ----------------
	set format {}
	set genolist [list]
	foreach analysis $analyses {
		set temp {}
		foreach field $formata(fields) formatnumber $formata(numbers) fpos $formata(fields,$analysis)  {
			if {[llength $fpos] == 1} {
				set value [list_subindex $lines $fpos]
			} else {
				# join  is always between a 1 (_ref) and A
				set v1 [lindex $lines 0 [lindex $fpos 0]]
				set v2 [list_subindex $lines [lindex $fpos 1]]
				set value [list $v1 {*}$v2]
			}
			if {$formatnumber in "A R" || ([isint $infonumber] && $infonumber > 1)} {
				if {[list_remove $value {}] eq ""} {
					set value {}
				} else {
					set value [join $value ,]
				}
			} else {
				set value [list_remove [list_remdup $value] {}]
			}
			lappend temp $value
		}
		lappend genolist $temp
	}
	set format GT
	set pos 0
	set useposs [list]
	foreach field $formata(fields) {
		if {[llength [list_remove [list_subindex $genolist $pos] {} ?]]} {
			append format :$field
			lappend useposs $pos
		}
		incr pos
	}
	set resultline $chromosome\t$chrpos\t$id\t$vcfref\t$vcfalt\t$qual\t$filter\t$info\t$format
	foreach analysis $analyses geno $genolist {
		set a1 . ; set a2 . ; set vcfphased /
		set altnr 1
		foreach line $lines ref $refs append $appends lref $linerefs lalt $linealts {
			foreach {alleleSeq1 alleleSeq2 genotypes sequenced zyg phased} [list_sub $line $formata(keyposs,$analysis)] break
			if {$alleleSeq1 eq $lref && $a1 eq "."} {
				set a1 0
			} elseif {$alleleSeq1 eq "$lalt"} {
				set a1 $altnr
			}
			if {$alleleSeq2 eq $lref && $a2 eq "."} {
				set a2 0
				if {$phased > 0} {set vcfphased |} else {set vcfphased /}
			} elseif {$alleleSeq2 eq "$lalt"} {
				set a2 $altnr
				if {$phased > 0} {set vcfphased |} else {set vcfphased /}
			}
			incr altnr
		}
		if {$a1 eq "."} {set a1 0} ; if {$a2 eq "."} {set a2 0}
		set gt $a1$vcfphased$a2
		set formatdata [list_sub $geno $useposs]
		set formatdata [join [list_change $formatdata {{} . ? .}] \:]
		regsub {(:\.)+$} $formatdata {} formatdata
		if {$formatdata ne ""} {
			append resultline \t$gt:$formatdata
		} else {
			append resultline \t$gt
		}
	}
	foreach line $extra {
		append resultline \n$line
	}
	return $resultline
}

proc cg_tsv2vcf {args} {
	set split 0
	set refseq {}
	set dbdir {}
	cg_options vcf2tsv args {
		-s - -split {
			set split $value
		}
		-dbdir {
			set dbdir $value
		}
		-refseq {
			set refseq $value
		}
		-sample {
			set sample $value
		}
	} {infile outfile} 0 2
	set refseq [refseq $refseq $dbdir]
	if {[info exists infile]} {
		set f [gzopen $infile]
	} else {
		set f stdin
	}
	if {[info exists outfile]} {
		set o [wgzopen $outfile]
	} else {
		set o stdout
	}

	array set conv_formata [list_reverse {
		AD alleledepth
		AD_ref alleledepth_ref
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
		tDP totaldepth
		PH phased
	}]
	set header [tsv_open $f comment]
	unset -nocomplain a
	set keepcomments {}
	foreach line [split $comment \n] {
		set line [split $line \t]
		if {[llength $line] == 1} {
			lappend keepcomments [lindex $line 0]
			continue
		}
		set key [string range [lindex $line 0] 1 end]
		if {[llength $line] > 2} {
			set value [lrange $line 1 end]
		} else {
			set value [lindex $line 1]
		}
		lappend a($key) $value
	}
	if {[info exists a(info)] && [regexp {original comments follow} $a(info)]} {
		set vcfinfo [vcf2tsvheader $keepcomments $header $split {} {}]
		foreach line $vcfinfo {
			set line [split $line \t]
			set key [string range [lindex $line 0] 1 end]
			if {[llength $line] > 2} {
				set value [lrange $line 1 end]
			} else {
				set value [lindex $line 1]
			}
			lappend a($key) $value
		}
		if {[info exists a(info)]} {set a(info) [list_remove $a(info) {tsv converted from vcf, original comments follow}]}
	}
	set poss [tsv_basicfields $header 6]
	set analyses [listanalyses $header {} analysisfields]
	if {![llength $analyses]} {
		if {[info exists sample]} {
			# keep manually set sample (using -sample)
		} elseif {[info exists a(samplename)]} {
			set sample $a(samplename)
		} elseif {[info exists infile]} {
			set sample [file_sample $infile]
		} else {
			set sample sample
		}
		set analyses [list $sample]
		set multicompar 0
	} else {
		set multicompar 1
	}
	if {[info exists a(fields)]} {
		# parse fields meta data
		# ----------------------
		unset -nocomplain fa
		set fheader [lindex $a(fields) 1]
		set fposs [list_cor $fheader {field number type description source}]
		if {-1 in [lrange $fposs 0 end-1]} {
			error "fields metadata table is missing fields: [list_sub {field number type description} [list_find $fposs -1]]"
		}
		foreach line [lrange $a(fields) 2 end] {
			set line [list_sub $line $fposs]
			set fa([lindex $line 0]) $line
		}
	}
	unset -nocomplain keyposa
	unset -nocomplain infoa
	unset -nocomplain formata
	set formata(fields) {}
	unset -nocomplain analysesa
	set formata(numbers) {}
	set infofields {}
	set infoposs {}
	set infoflags {}
	set infonumbers {}
	set pos -1
	foreach field $header {
		incr pos
		set hasanalysis [regexp {^([^-]+)-(.+)$} $field temp field analysis]
		set fielddata [tsv2vcf_fielddata $field fa]
		if {!$multicompar} {
			if {[lindex $fielddata end] in {format geno}} {
				set hasanalysis 1
				set analysis $sample
			}
		}
		if {$hasanalysis} {
			# sample specific field
			set analysesa($analysis) 1
			set fieldname [get conv_formata($field) $field]
			if {$field in {alleleSeq1 alleleSeq2 genotypes sequenced zyg phased}} {
				set keyposa($field,$analysis) $pos
				# unset -nocomplain fa($field)
				continue
			}
			if {[regexp {^(.*)_ref$} $field temp base]} {
				set join 0
				if {[info exists fa($base)] && [regexp {^alleles only values of:} [lindex $fa($base) 3]]} {
					set join 1
				}
				if {!$join && $field in "alleledepth_ref RPA_ref" && ([lsearch $header $base-$analysis] != -1 || [lsearch $header $base] != -1)} {
					set join 1
				}
				if {$join} {
					regsub {_ref$} $fieldname {} basename
					if {![info exists formata(header,$basename)]} {
						foreach {temp number type description} $fielddata break
						# unset fa($base) ; unset -nocomplain fa($field)
						set number R
						regsub {^reference only value of: } $description {} description
						lappend formata(fields) $basename
						lappend formata(numbers) R
						set formata(header,$basename) "##FORMAT=<ID=$basename,Number=$number,Type=$type,Description=\"$description\">"
					}
					set fieldname $basename
				}
			}
			if {![info exists formata(header,$fieldname)]} {
				lappend formata(fields) $fieldname
				foreach {temp number type description} $fielddata break
				lappend formata(numbers) $number
				set formata(header,$fieldname) "##FORMAT=<ID=$fieldname,Number=$number,Type=$type,Description=\"$description\">"
			}
			lappend formata(pos,$fieldname,$analysis) $pos
		} else {
			# general field
			set fieldname [get conv_formata($field) $field]
			if {$field in {chromosome begin end type ref alt name quality filter alleleSeq1 alleleSeq2 genotypes sequenced phased}} {
				set keyposa($field) $pos
				# unset -nocomplain fa($field)
				continue
			}
			if {$fieldname eq "totalcoverage"} {
				set fieldname DP
			}
			set fielddata [tsv2vcf_fielddata $field fa]
			foreach {temp number type description} $fielddata break
			lappend infofields $fieldname
			lappend infoposs $pos
			lappend infonumbers $number
			set infoa($fieldname) "##INFO=<ID=$fieldname,Number=$number,Type=$type,Description=\"$description\">"
			if {$type eq "Flag"} {lappend infoflags 1} else {lappend infoflags 0}
		}
	}
	# prepare for format output
	# -------------------------
	set analyses [array names analysesa]
	unset -nocomplain analysesa
	foreach analysis $analyses {
		set formata(keyposs,$analysis) [list]
		set formata(fields,$analysis) [list]
		foreach key {alleleSeq1 alleleSeq2 genotypes sequenced zyg phased} {
			lappend formata(keyposs,$analysis) [get keyposa($key,$analysis) -1]
		}
		foreach fieldname $formata(fields) {
			lappend formata(fields,$analysis) [get formata(pos,$fieldname,$analysis) -1]
		}
	}
	set format GT\;[join $formata(fields) \;]
	# output vcf header
	# -----------------
	puts $o "##fileformat=VCFv4.2"
	# tsv2vcf_outputheaderfield will print a line for the given key only if it exists
	tsv2vcf_outputheaderfield $o fileDate a
	tsv2vcf_outputheaderfield $o source a
	tsv2vcf_outputheaderfield $o reference a
	tsv2vcf_outputheaderfield $o phasing a
	foreach key {vcf_fileformat {} { ----} fields filetype fileversion numsamples} {
		unset -nocomplain a($key)
	}
	foreach field [lsort -dict [array names a vcf_*]] {
		puts $o "##[string range $field 4 end]=[lindex $a($field) 0]"
		unset a($field)
	}
	tsv2vcf_outputheaderfield $o FILTER a
	puts $o {##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">}
	foreach name $formata(fields) {
		puts $o $formata(header,$name)
	}
	foreach name $infofields {
		puts $o $infoa($name)
	}
	foreach key [array names a] {
		tsv2vcf_outputheaderfield $o $key a
	}
	# puts $o $keepcomments
	set vcfheader "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
	foreach analysis $analyses {
		append vcfheader "\t$analysis"
	}
	puts $o $vcfheader

	#
	# print variants
	# ==============
	catch {genome_close $genomef} ; set genomef [genome_open $refseq]
	set prevchr {}
	unset -nocomplain allelelista
	while 1 {
		if {[gets $f line] == -1} {
			set pos -1
			set print 1
		} else {
			set print $split
			set line [split $line \t]
			set var [tsv2vcf_getvar $line $poss keyposa $genomef]
			foreach {chr pos} $var break
			if {$chr ne $prevchr} {
				set print 1
				set prevchr $chr
			}
		}
		if {$split} {
			if {$pos == -1} break
			puts $o [tsv2vcf_printlines [list [list $var $line]] $infofields $infoposs $infoflags $infonumbers $analyses formata]
		} else {
			foreach p [ssort -natural [array names allelelista]] {
				if {$print || $p < $pos} {
					set lines $allelelista($p)
					puts $o [tsv2vcf_printlines $lines $infofields $infoposs $infoflags $infonumbers $analyses formata]
					unset allelelista($p)
				}
			}
			# break if file is finisehd
			if {$pos == -1} break
			# otherwise add to list to process, use pos +1 to make sure sub after snp gets processed with previous set
			lappend allelelista([expr {$pos+1}]) [list $var $line]
		}
	}
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}
