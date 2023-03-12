proc empty_vcf {file {samples SAMPLE} {extracomment {}} {extrainfo {}}} {
	set f [wgzopen $file w]
	if {$extrainfo ne ""} {
		regsub -all {\n\t*} [string trim $extrainfo] "\n\t\t" temp
		set extrainfo \n\t\t$temp
	}
	puts $f [deindent [string_change {
		##fileformat=VCFv4.0
		##fileDate=20090805
		##source=myImputationProgramV3.1
		##phasing=partial
		##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
		##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
		##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
		##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">@extrainfo@
		##FILTER=<ID=q10,Description="Quality below 10">
		##FILTER=<ID=s50,Description="Less than 50% of samples have data">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=TE,Number=A,Type=Integer,Description="test for alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
	} [list @extrainfo@ $extrainfo]]]
	if {$extracomment ne ""} {puts -nonewline $f [string trim $extracomment]\n}
	set header {CHROM POS ID REF ALT QUAL FILTER INFO FORMAT}
	lappend header {*}$samples
	puts $f \#[join $header \t]
	close $f
}

proc empty_gvcf {file {samples SAMPLE} {extracomment {}} {extrainfo {}}} {
	set f [wgzopen $file w]
	if {$extrainfo ne ""} {
		regsub -all {\n\t*} [string trim $extrainfo] "\n\t\t" temp
		set extrainfo \n\t\t$temp
	}
	puts $f [deindent [string_change {
		##fileformat=VCFv4.2
		##FILTER=<ID=PASS,Description="All filters passed">
		##FILTER=<ID=LowQual,Description="Low quality variant">
		##FILTER=<ID=RefCall,Description="Reference call">
		##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
		##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">@extrainfo@
		##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
		##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
		##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
		##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
		##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">
	} [list @extrainfo@ $extrainfo]]]
	if {$extracomment ne ""} {puts -nonewline $f [string trim $extracomment]\n}
	set header {CHROM POS ID REF ALT QUAL FILTER INFO FORMAT}
	lappend header {*}$samples
	puts $f \#[join $header \t]
	close $f
}

