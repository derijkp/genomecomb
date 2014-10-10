proc code2typevar {typeaVar} {
	upvar $typeaVar typea
	foreach {code type number descr} {
		id Integer 1 "id"
		chromosome String 1 "chromosome name"
		begin Integer 1 "begin of feature (half open format)"
		end Integer 1 "end of feature (half open format)"
		type String 1 "type of feature (snp, ins, del, sub, ...)"
		ref String 1 "reference sequence"
		alt String . "alternative sequence(s)"
		quality Float 1 "Variant Quality"
		AC Integer A "Allele count in genotypes for each ALT allele, in the same order as listed"
		allelecount Integer A "AC Allele count in genotypes for each ALT allele, in the same order as listed"
		AC1 Float 1 "Max-likelihood estimate of the first ALT allele count (no HWE assumption)"
		AD Integer . "Allelic depths for the ref and alt alleles in the order listed"
		AF Float A "Allele Frequency, for each ALT allele, in the same order as listed"
		frequency Float A "AF Allele Frequency, for each ALT allele, in the same order as listed"
		AF1 Float 1 "Max-likelihood estimate of the first ALT allele frequency (assuming HWE)"
		AN Integer 1 "Total number of alleles in called genotypes"
		totalallelecount Integer 1 "AN Total number of alleles in called genotypes"
		Ancestralallele String 1 "Ancestral allele"
		BaseQRankSum Float 1 "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities"
		CGT String 1 "The most probable constrained genotype configuration in the trio"
		CLR Float 1 "Log ratio of genotype likelihoods with and without the constraint"
		DP Integer 1 "# read depth"
		coverage Integer 1 "DP (geno) Approximate read depth"
		totalcoverage Integer 1 "DP (info) Raw read depth (total)"
		DP4 Integer 4 "# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"
		DS Flag 0 "Were any of the samples downsampled?"
		DV Integer 1 "# high-quality non-reference bases"
		Dels Float 1 "Fraction of Reads Containing Spanning Deletions"
		FQ Float 1 "Phred probability of all samples being the same"
		FS Float 1 "Phred-scaled p-value using Fisher's exact test to detect strand bias"
		G3 Float 3 "ML estimate of genotype frequencies"
		GL Float 3 "Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)"
		loglikelihood Float 3 "GL Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)"
		GQ Float 1 "Genotype Quality"
		genoqual Float 1 "GQ Genotype Quality"
		haploqual Float 1 "HQ Haplotype Quality"
		GT String 1 "Genotype"
		genotype String 1 "GT Genotype"
		HWE Float 1 "Chi^2 based HWE test P-value based on G3"
		HaplotypeScore Float 1 "Consistency of the site with at most two segregating haplotypes"
		INDEL Flag 0 "Indicates that the variant is an INDEL."
		IS Float 2 "Maximum number of reads supporting an indel and fraction of indel reads"
		InbreedingCoeff Float 1 "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"
		MDV Integer 1 "Maximum number of high-quality nonRef reads in samples"
		MLEAC Float A "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed"
		MLEAF Float A "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"
		MQ Float 1 "Root-mean-square mapping quality of covering reads"
		MQ0 Float 1 "Total Mapping Quality Zero Reads"
		MQRankSum Float 1 "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"
		NDA Integer 1 "Number of alternate alleles discovered (but not necessarily genotyped) at this site"
		PC2 Float 2 "Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2."
		PCHI2 Float 1 "Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples."
		PL Float G "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"
		PR Integer 1 "# permutations yielding a smaller PCHI2."
		PV4 Float 4 "P-values for strand bias, baseQ bias, mapQ bias and tail distance bias"
		QBD Float 1 "Quality by Depth: QUAL/#reads"
		QCHI2 Float 1 "Phred scaled PCHI2."
		QD Float 1 "Variant Confidence/Quality by Depth"
		RPA Integer . "Number of times tandem repeat unit is repeated, for each allele (including reference)"
		RPB Float 1 "Read Position Bias"
		RU String 1 "Tandem repeat unit (bases)"
		ReadPosRankSum Float 1 "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"
		SP Float 1 "Phred-scaled strand bias P-value"
		STR Flag 0 "Variant is a short tandem repeat"
		UGT String 1 "The most probable unconstrained genotype configuration in the trio"
		VDB Float 1 "Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken."
		dbsnp String 1 "DB dbsnp xref"
		Hapmap2 String 1 "H2 hapmap xref"
		filter String . "filter"
	} {
		switch $type {
			String {set type text}
			Integer {set type int}
			Flag {set type char}
			default {set type [string tolower $type]}
		}
		if {$number eq "1"} {
			set typea($code) $type
		} else {
			set typea($code) ${type}_list
		}
	}
}

proc sqlite_tablefromtsv {table tempfile {replace {}}} {
	array set replacea $replace
	code2typevar typea
	set f [open $tempfile]
	set header [split [gets $f] \t]
	set line [split [gets $f] \t]
	close $f
	set sql ""
	foreach field $header v $line {
		if {[info exists replacea($field)]} {
			lappend sql $replacea($field)
			continue
		} elseif {$field eq "id"} {
			lappend sql "\"id\" integer primary key"
			continue
		} elseif {[info exists typea($field)]} {
			set type $typea($field)
		} elseif {[isint $v]} {
			set type integer
		} elseif {[isdouble $v]} {
			set type real
		} else {
			set type text
		}
		lappend sql "\"$field\" $type"
	}
	set sql "create table \"$table\" ([join $sql ,])"
}

proc cg_tosqlite {args} {
	global scriptname action
	if {[llength $args] != 3} {
		puts stderr "format is: $scriptname $action dbfile table tsvfile"
		exit 1
	}
	foreach {dbfile table tsvfile} $args break
	set sql [sqlite_tablefromtsv $table $tsvfile]
	catch {package require dbi}
	package require dbi_sqlite3
	dbi_sqlite3 db
	db create $dbfile
	db open $dbfile
	db exec $sql
	db import abort $table $tsvfile \t {}
	db close
}
