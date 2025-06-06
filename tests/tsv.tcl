#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc makepastetest {num} {
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
	}
	set files tmp/vars1.tsv
	for {set i 1} {$i < $num} {incr i} {
		write_tab tmp/sample$i.tsv [subst {
			#
			# rc
			zyg-sample$i	value-sample$i
			m	$i
		}]
		lappend files tmp/sample$i.tsv
	}
	# expected
	set f [open tmp/expected.tsv w]
	puts $f "#	varcomment"
	puts -nonewline $f "chromosome	begin	end	type	ref	alt"
	for {set i 1} {$i < $num} {incr i} {
		puts -nonewline $f "	zyg-sample$i	value-sample$i"
	}
	puts -nonewline $f "\nchr1	4001	4002	snp	A	G,C"
	for {set i 1} {$i < $num} {incr i} {
		puts -nonewline $f "	m	$i"
	}
	puts $f ""
	close $f
	return $files
}

test tsv_open {tsv} {
	catch {close $f}
	set f [open data/reg1.tsv]
	set header [tsv_open $f keepheader]
	set line [gets $f]
	close $f
	list $header $keepheader $line
} {{chromosome test begin end} {} {1	t	10	20}}

test tsv_open {cg} {
	catch {close $f}
	set f [open data/cgtest.tsv]
	set header [tsv_open $f keepheader]
	set line [gets $f]
	close $f
	list $header $keepheader $line
} {{offset refScore uniqueSequenceCoverage weightSumSequenceCoverage gcCorrectedCoverage grossWeightSumSequenceCoverage} {#ASSEMBLY_ID	GS19240-1100-37-ASM
#CHROMOSOME	chr22
#FORMAT_VERSION	1.5
#GENERATED_AT	2010-Oct-29 08:59:22.556619
#GENERATED_BY	ExportReferenceSupport
#SAMPLE	GS00028-DNA_C01
#SOFTWARE_VERSION	1.10.0.17
#TYPE	REFMETRICS
#
} {16050000	0	0	0	0	0}}

test tsv_open {rtg} {
	catch {close $f}
	set f [open data/rtgsnps.tsv]
	set header [tsv_open $f keepheader]
	tsv_hcheader $f keepheader header
	set line [gets $f]
	close $f
	list $header $keepheader $line
} {{name position type reference prediction posterior coverage correction support_statistics} {#Version v2.0-RC build 27737 (2010-05-06), SNP output v2.0
#CL	snp -t /rtgshare/data/human/sdf/hg18 -o snps_GS00102-DNA-D06_m4_u4 --max-as-mated=4 --max-as-unmated=4 --max-ih=1 -m cg_errors --no-complex-calls /rtgshare3/users/richard/vib_20100324_GS00102-DNA-D06-hg18ref/map_GS000004945/mated.sam.gz /rtgshare3/users/richard/vib_20100324_GS00102-DNA-D06-hg18ref/map_GS000004945/unmated.sam.gz ...
} {chr1	231	e	C	A:C	1.5	109	4.982	A	17	0.655	C	91	4.004	T	1	0.323}}

test tsv_open {vcf} {
	catch {close $f}
	set f [open data/test.vcf]
	set header [tsv_open $f keepheader]
	set line [gets $f]
	close $f
	list $header $keepheader $line
} {{CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003} {##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=TE,Number=A,Type=Integer,Description="test for alt alleles in the order listed">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#} {20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.}}

test tsv_cat {one} {
	cg cat data/reg1b.tsv
} {chromosome	test	begin	end
1	t	10	20
1	t	50	60}

test tsv_cat {two same header -c 1} {
	cg cat -c 1 data/reg2.tsv data/reg4.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# ++++ data/reg2.tsv ++++
		# comments added
		#
		# ++++ data/reg4.tsv ++++
		## header
		##
		test	chromosome	begin	end	test2
		t	1	15	25	t2
		t	1	45	55	t2
		t	2	150	160	t2
		t	2	170	180	t2
		t	2	300	400	t2
		t	2	400	500	t2
		t	3	1000	1100	t2
		t	M	10	25	t2
		t	X	90	200	t2
		t	Y	1010	1900	t2
		t	chr1	15	25	t2
		t	chr1	45	55	t2
		t	chr2	150	160	t2
		t	chr2	170	180	t2
		t	chr2	300	400	t2
		t	chr2	400	500	t2
		t	chr3	1000	1100	t2
		t	chrM	10	25	t2
		t	chrX	90	200	t2
		t	chrY	1010	1900	t2
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two same header} {
	cg cat data/reg2.tsv data/reg4.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		test	chromosome	begin	end	test2
		t	1	15	25	t2
		t	1	45	55	t2
		t	2	150	160	t2
		t	2	170	180	t2
		t	2	300	400	t2
		t	2	400	500	t2
		t	3	1000	1100	t2
		t	M	10	25	t2
		t	X	90	200	t2
		t	Y	1010	1900	t2
		t	chr1	15	25	t2
		t	chr1	45	55	t2
		t	chr2	150	160	t2
		t	chr2	170	180	t2
		t	chr2	300	400	t2
		t	chr2	400	500	t2
		t	chr3	1000	1100	t2
		t	chrM	10	25	t2
		t	chrX	90	200	t2
		t	chrY	1010	1900	t2
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two same header -n file -c 0} {
	cg cat -c 0 -n file data/reg2.tsv data/reg4.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# comments added
		#
		## header
		##
		test	chromosome	begin	end	test2	file
		t	1	15	25	t2	reg2.tsv
		t	1	45	55	t2	reg2.tsv
		t	2	150	160	t2	reg2.tsv
		t	2	170	180	t2	reg2.tsv
		t	2	300	400	t2	reg2.tsv
		t	2	400	500	t2	reg2.tsv
		t	3	1000	1100	t2	reg2.tsv
		t	M	10	25	t2	reg2.tsv
		t	X	90	200	t2	reg2.tsv
		t	Y	1010	1900	t2	reg2.tsv
		t	chr1	15	25	t2	reg4.tsv
		t	chr1	45	55	t2	reg4.tsv
		t	chr2	150	160	t2	reg4.tsv
		t	chr2	170	180	t2	reg4.tsv
		t	chr2	300	400	t2	reg4.tsv
		t	chr2	400	500	t2	reg4.tsv
		t	chr3	1000	1100	t2	reg4.tsv
		t	chrM	10	25	t2	reg4.tsv
		t	chrX	90	200	t2	reg4.tsv
		t	chrY	1010	1900	t2	reg4.tsv
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two diff header} {
	exec cg cat data/reg1b.tsv data/reg2.tsv
} {headers do not match, use -f to force or -m to merge (at file data/reg2.tsv)} error

test tsv_cat {two diff header} {
	write_tab tmp/reg.tsv {
		begin end
		1 2
	}
	exec cg cat data/reg1b.tsv tmp/reg.tsv
} {headers do not match, use -f to force or -m to merge (at file tmp/reg.tsv)} error

test tsv_cat {two diff header -m -c 1} {
	cg cat -m -c 1 data/reg2.tsv data/reg1b.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# ++++ data/reg2.tsv ++++
		# ++ test chromosome begin end test2
		# comments added
		#
		# ++++ data/reg1b.tsv ++++
		# ++ chromosome test begin end
		test	chromosome	begin	end	test2
		t	1	15	25	t2
		t	1	45	55	t2
		t	2	150	160	t2
		t	2	170	180	t2
		t	2	300	400	t2
		t	2	400	500	t2
		t	3	1000	1100	t2
		t	M	10	25	t2
		t	X	90	200	t2
		t	Y	1010	1900	t2
		t	1	10	20	
		t	1	50	60	
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two diff header -m 1} {
	cg cat -m 1 -catfiles 1 data/reg2.tsv data/reg1b.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		#catfiles	data/reg2.tsv
		#catfiles	data/reg1b.tsv
		test	chromosome	begin	end	test2
		t	1	15	25	t2
		t	1	45	55	t2
		t	2	150	160	t2
		t	2	170	180	t2
		t	2	300	400	t2
		t	2	400	500	t2
		t	3	1000	1100	t2
		t	M	10	25	t2
		t	X	90	200	t2
		t	Y	1010	1900	t2
		t	1	10	20	
		t	1	50	60	
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two diff header -m, -n -c 0} {
	cg cat -n file -c 0 -m 1 data/reg2.tsv data/reg1b.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# comments added
		#
		test	chromosome	begin	end	test2	file
		t	1	15	25	t2	reg2.tsv
		t	1	45	55	t2	reg2.tsv
		t	2	150	160	t2	reg2.tsv
		t	2	170	180	t2	reg2.tsv
		t	2	300	400	t2	reg2.tsv
		t	2	400	500	t2	reg2.tsv
		t	3	1000	1100	t2	reg2.tsv
		t	M	10	25	t2	reg2.tsv
		t	X	90	200	t2	reg2.tsv
		t	Y	1010	1900	t2	reg2.tsv
		t	1	10	20		reg1b.tsv
		t	1	50	60		reg1b.tsv
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two diff header -m -c 1} {
	cg cat -m -c 1 data/reg1.tsv data/reg2.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# ++++ data/reg1.tsv ++++
		# ++ chromosome test begin end
		# ++++ data/reg2.tsv ++++
		# ++ test chromosome begin end test2
		# comments added
		#
		chromosome	test	begin	end	test2
		1	t	10	20	
		1	t	50	60	
		2	t	100	200	
		2	t	450	480	
		3	t	1000	1100	
		3	t	2000	2100	
		M	t	10	20	
		X	t	100	200	
		Y	t	1000	2000	
		1	t	15	25	t2
		1	t	45	55	t2
		2	t	150	160	t2
		2	t	170	180	t2
		2	t	300	400	t2
		2	t	400	500	t2
		3	t	1000	1100	t2
		M	t	10	25	t2
		X	t	90	200	t2
		Y	t	1010	1900	t2
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two diff header -f -c 1} {
	# ! This is not something you would want to do normally
	cg cat -f -c 1 data/reg2.tsv data/reg1b.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# ++++ data/reg2.tsv ++++
		# ++ test chromosome begin end test2
		# comments added
		#
		# ++++ data/reg1b.tsv ++++
		# ++ chromosome test begin end
		test	chromosome	begin	end	test2
		t	1	15	25	t2
		t	1	45	55	t2
		t	2	150	160	t2
		t	2	170	180	t2
		t	2	300	400	t2
		t	2	400	500	t2
		t	3	1000	1100	t2
		t	M	10	25	t2
		t	X	90	200	t2
		t	Y	1010	1900	t2
		1	t	10	20
		1	t	50	60
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {two diff header -fields} {
	write_tab tmp/reg.tsv {
		test    chromosome      begin   end     test2
		t       1       15      25      t2
		t       1       45      55      t2
	}
	exec cg cat -fields {chromosome begin end test} -c 1 data/reg1b.tsv tmp/reg.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		# ++++ data/reg1b.tsv ++++
		# ++++ tmp/reg.tsv ++++
		chromosome	begin	end	test
		1	10	20	t
		1	50	60	t
		1	15	25	t
		1	45	55	t
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {-n file} {
	write_tab tmp/reg1.tsv {
		chromosome	begin	end
		1	15	25
	}
	write_tab tmp/reg2.tsv {
		chromosome	begin	end
		1	45	55
	}
	exec cg cat -n file -c 0 tmp/reg1.tsv tmp/reg2.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		chromosome	begin	end	file
		1	15	25	reg1.tsv
		1	45	55	reg2.tsv
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {-m 1 -n file} {
	write_tab tmp/reg1.tsv {
		chromosome	test begin	end
		1	test	15	25
	}
	write_tab tmp/reg2.tsv {
		chromosome	begin	end
		1	45	55
	}
	exec cg cat -m 1 -n file -c 0 tmp/reg1.tsv tmp/reg2.tsv > tmp/result.tsv
	write_deindent tmp/expected.tsv {
		chromosome	test	begin	end	file
		1	test	15	25	reg1.tsv
		1		45	55	reg2.tsv
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {cg cat -m 1 -c m} {
	write_deindent tmp/var1.tsv {
		#filetype	tsv/varfile
		#fileversion	0.99
		#split	1
		#refseq	genome_hg19.ifas
		#numsamples	1
		#extra	1
		#info	test 1
		#info	test 2
		#plain comment
		#samplename	NA19240m
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	test	1	String	test	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#contig	table
		#contig	ID	length
		#contig	chr1	249250621
		chromosome	test	begin	end
		chr1	test	15	25
	}
	write_deindent tmp/var2.tsv {
		#filetype	tsv/varfile
		#fileversion	0.99
		#split	1
		#info	tsv converted from vcf
		#info	test 3
		#refseq	genome_hg19.ifas
		#numsamples	1
		#extra2	2
		#samplename	NA19240m
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#fields	type	1	String	Type of feature (snp,del,ins,...)	var
		#contig	table
		#contig	ID	length
		#contig	chr2	243199373
		chromosome	begin	end	type
		chr2	45	55	del
	}
	write_deindent tmp/expected.tsv {
		#filetype	tsv/varfile
		#fileversion	0.99
		#split	1
		#info	test 1
		#info	test 2
		#info	tsv converted from vcf
		#info	test 3
		#refseq	genome_hg19.ifas
		#numsamples	1
		#samplename	NA19240m
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	test	1	String	test	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#fields	type	1	String	Type of feature (snp,del,ins,...)	var
		#contig	table
		#contig	ID	length
		#contig	chr1	249250621
		#contig	chr2	243199373
		#extra	1
		#extra2	2
		chromosome	test	begin	end	type
		chr1	test	15	25	
		chr2		45	55	del
	}
	exec cg cat -m 1 -c m tmp/var1.tsv tmp/var2.tsv > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {cg cat -m 1 -c m} {
	write_deindent tmp/var1.tsv {
		chromosome	test	begin	end
		chr1	test	15	25
	}
	write_deindent tmp/var2.tsv {
		chromosome	begin	end	type
		chr2	45	55	del
	}
	write_deindent tmp/expected.tsv {
		chromosome	test	begin	end	type
		chr1	test	15	25	
		chr2		45	55	del
	}
	exec cg cat -m 1 -c m tmp/var1.tsv tmp/var2.tsv > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {cg cat -m 1 -c 1} {
	write_deindent tmp/var1.tsv {
		chromosome	test	begin	end
		chr1	test	15	25
	}
	write_deindent tmp/var2.tsv {
		chromosome	begin	end	type
		chr2	45	55	del
	}
	write_deindent tmp/expected.tsv {
		chromosome	test	begin	end	type
		chr1	test	15	25	
		chr2		45	55	del
	}
	exec cg cat -m 1 -c m tmp/var1.tsv tmp/var2.tsv > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {cg cat -m 1 -c m -sample testsample} {
	write_deindent tmp/var1.tsv {
		#filetype	tsv/varfile
		#fileversion	0.99
		#split	1
		#refseq	genome_hg19.ifas
		#numsamples	1
		#extra	1
		#info	test 1
		#info	test 2
		#plain comment
		#samplename	x1
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	test	1	String	test	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#contig	table
		#contig	ID	length
		#contig	chr1	249250621
		chromosome	test	begin	end
		chr1	test	15	25
	}
	write_deindent tmp/var2.tsv {
		#filetype	tsv/varfile
		#fileversion	0.99
		#split	1
		#info	tsv converted from vcf
		#info	test 3
		#refseq	genome_hg19.ifas
		#numsamples	1
		#extra2	2
		#samplename	x2
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#fields	type	1	String	Type of feature (snp,del,ins,...)	var
		#contig	table
		#contig	ID	length
		#contig	chr2	243199373
		chromosome	begin	end	type
		chr2	45	55	del
	}
	write_deindent tmp/expected.tsv {
		#filetype	tsv/varfile
		#fileversion	0.99
		#split	1
		#info	test 1
		#info	test 2
		#info	tsv converted from vcf
		#info	test 3
		#refseq	genome_hg19.ifas
		#numsamples	1
		#samplename	testsample
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	test	1	String	test	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#fields	type	1	String	Type of feature (snp,del,ins,...)	var
		#contig	table
		#contig	ID	length
		#contig	chr1	249250621
		#contig	chr2	243199373
		#extra	1
		#extra2	2
		chromosome	test	begin	end	type
		chr1	test	15	25	
		chr2		45	55	del
	}
	exec cg cat -m 1 -c m -sample testsample tmp/var1.tsv tmp/var2.tsv > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_cat {cg cat empty file} {
	write_deindent tmp/var1.tsv {
		chromosome	test	begin	end
		chr1	test	15	25
	}
	write_deindent tmp/var2.tsv {
	}
	exec cg cat tmp/var1.tsv tmp/var2.tsv > tmp/result.tsv
	exec diff tmp/result.tsv tmp/var1.tsv
} {}

test tsv_cat {cg cat compressed empty file} {
	write_deindent tmp/var1.tsv {
		chromosome	test	begin	end
		chr1	test	15	25
	}
	file_write tmp/var2.tsv {}
	cg zst tmp/var2.tsv
	exec cg cat tmp/var1.tsv tmp/var2.tsv.zst > tmp/result.tsv
	exec diff tmp/result.tsv tmp/var1.tsv
} {}

test check_sort {sort error 1 in vars} {
	exec cg checksort data/vars_sorterror1.tsv
} {error in file data/vars_sorterror1.tsv: file is not correctly sorted (sort correctly using "cg select -s -")
chr10:43198434-43198435:snp:G came before chr3:52847042-52847060:del:} error

test check_sort {sort error 2 in vars} {
	exec cg checksort data/vars_sorterror2.tsv
} {error in file data/vars_sorterror2.tsv: file is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:snp:G came before chr3:52847042-52847060:del:} error

test check_sort {sort error 3 in vars} {
	exec cg checksort data/vars_sorterror3.tsv
} {error in file data/vars_sorterror3.tsv: file is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847310:del: came before chr3:52847303-52847304:snp:G} error

test check_sort {sort error 4 in vars} {
	exec cg checksort data/vars_sorterror4.tsv
} {error in file data/vars_sorterror4.tsv: file is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:snp:G came before chr3:52847303-52847304:ins:G} error

test check_sort {sort error 5 in vars} {
	cg checksort data/vars_sorterror5.tsv
} {error in file data/vars_sorterror5.tsv: file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test check_sort {sort error alt} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	4001	4002	snp	A	G
		chr1	204434325	204434325	ins	{}	CT
		chr1	204434325	204434326	snp	A	G
		chr1	204434325	204434326	snp	A	C
	}
	cg checksort tmp/vars.tsv
} {error in file tmp/vars.tsv: file is not correctly sorted (sort correctly using "cg select -s -")
chr1:204434325-204434326:snp:G came before chr1:204434325-204434326:snp:C} error match

test check_sort {sort error alt sub} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	4001	4002	snp	A	G
		chr1	204434325	204434325	ins	{}	CT
		chr1	204434325	204434326	snp	A	G
		chr1	204434325	204434326	sub	A	CTA
	}
	cg checksort tmp/vars.tsv
} {}

test check_sort {sort error alt sub} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	4001	4002	snp	A	G
		chr1	204434325	204434325	ins	{}	CT
		chr1	204434325	204434326	sub	A	G
		chr1	204434325	204434326	sub	A	CTA
	}
	cg checksort tmp/vars.tsv
} {error in file tmp/vars.tsv: file is not correctly sorted (sort correctly using "cg select -s -")
chr1:204434325-204434326:sub:G came before chr1:204434325-204434326:sub:CTA} error match

test indexdir {basic same name different dirs} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	file mkdir tmp/tmp
	write_tab tmp/tmp/vars1.tsv {
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
	}
	set result {}
	set d1 [indexdir tmp/vars1.tsv]
	set d2 [indexdir tmp/tmp/vars1.tsv]
	set d12 [indexdir tmp/vars1.tsv]
	set d22 [indexdir tmp/tmp/vars1.tsv]
	list [expr {$d1 eq $d12}] [expr {$d2 eq $d22}]  [expr {$d1 eq $d2}]
} {1 1 0}

test indexdir {basic same name different dirs compressed} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	file mkdir tmp/tmp
	write_tab tmp/tmp/vars1.tsv {
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
	}
	cg razip tmp/vars1.tsv
	set result {}
	set d1 [indexdir tmp/vars1.tsv]
	set d2 [indexdir tmp/tmp/vars1.tsv]
	set d12 [indexdir tmp/vars1.tsv]
	set d22 [indexdir tmp/tmp/vars1.tsv]
	list [expr {$d1 eq $d12}] [expr {$d2 eq $d22}]  [expr {$d1 eq $d2}]
} {1 1 0}

test indexdir_cache {tsv_varsfile} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	set varsfile [tsv_varsfile tmp/vars1.tsv]
	cg select -f {chromosome begin end type ref alt} data/vars1.tsv tmp/test
	exec diff $varsfile tmp/test
} {}

test indexdir_cache {tsv_varsfile compressed} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	cg razip tmp/vars1.tsv
	set varsfile [tsv_varsfile tmp/vars1.tsv.rz]
	cg select -f {chromosome begin end type ref alt} data/vars1.tsv tmp/test
	exec diff $varsfile tmp/test
} {}

test indexdir_cache {tsv_varsfile compressed use plain filename} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	cg razip tmp/vars1.tsv
	set varsfile [tsv_varsfile tmp/vars1.tsv]
	cg select -f {chromosome begin end type ref alt} data/vars1.tsv tmp/test
	exec diff $varsfile tmp/test
} {}

test indexdir_cache {tsv_varsfile fix hang compressed but empty} {
	test_cleantmp
	file_write_gz tmp/vars1.tsv.zst {}
	set varsfile [tsv_varsfile tmp/vars1.tsv.zst]
	file_read $varsfile
} {chromosome	begin	end	type	ref	alt
}

test indexdir_cache {varsfile fix hang compressed but empty} {
	test_cleantmp
	file_write_gz tmp/vars1.tsv.zst {}
	exec varsfile tmp/vars1.tsv.zst
	file_write tmp/vars1.tsv.gz {}
	exec varsfile tmp/vars1.tsv.gz
} {chromosome	begin	end	type	ref	alt}

test indexdir_cache {tsv_count} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	tsv_count tmp/vars1.tsv
} 14

test indexdir_cache {tsv_count compressed} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	cg razip tmp/vars1.tsv
	tsv_count tmp/vars1.tsv.rz
} 14

test indexdir_cache {tsv_count compressed gzfile} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	cg razip tmp/vars1.tsv
	tsv_count tmp/vars1.tsv
} 14

test indexdir_cache {tsv_count invalid cache} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	file mkdir tmp/vars1.tsv.index
	file_write tmp/vars1.tsv.index/vars.tsv.count 10
	after 1000
	file mtime tmp/vars1.tsv [clock seconds]
	list [tsv_count tmp/vars1.tsv] [file_read tmp/vars1.tsv.index/vars.tsv.count]
} {14 14}

test indexdir_cache {tsv_count invalid cache that cannot be overwritten in sib indexdir} {
	test_cleantmp
	file copy data/vars1.tsv tmp/vars1.tsv
	file mkdir tmp/vars1.tsv.index
	file_write tmp/vars1.tsv.index/vars.tsv.count 10
	file mtime tmp/vars1.tsv [clock seconds]
	file mtime tmp/vars1.tsv.index/vars.tsv.count [expr {[file mtime tmp/vars1.tsv] - 1000}]
	file attributes tmp/vars1.tsv.index -permissions ugo-xw
	set countfile [indexdir_file tmp/vars1.tsv vars.tsv.count ok]
	set result {}
	lappend result [tsv_count tmp/vars1.tsv]
	file attributes tmp/vars1.tsv.index -permissions u+xw
	lappend result [file_read tmp/vars1.tsv.index/vars.tsv.count]
	set result  
} {14 10}

test tsv_split {split} {
	test_cleantmp
	exec cg split data/reg4.tsv tmp/split- .tsv
	list [bsort [glob tmp/*]] [file_read tmp/split-3.tsv]
} {{tmp/split-1.tsv tmp/split-2.tsv tmp/split-3.tsv tmp/split-M.tsv tmp/split-X.tsv tmp/split-Y.tsv} {## header
##
test	chromosome	begin	end	test2
t	chr3	1000	1100	t2
}}

test tsv_split {split compressed} {
	test_cleantmp
	exec cg split data/reg4.tsv tmp/split- .tsv.zst
	exec cg zcat tmp/split-3.tsv.zst > tmp/split-3.tsv
	list [bsort [glob tmp/*]] [file_read tmp/split-3.tsv]
} {{tmp/split-1.tsv.zst tmp/split-2.tsv.zst tmp/split-3.tsv tmp/split-3.tsv.zst tmp/split-M.tsv.zst tmp/split-X.tsv.zst tmp/split-Y.tsv.zst} {## header
##
test	chromosome	begin	end	test2
t	chr3	1000	1100	t2
}}

test tsv_split {split compressed sorted} {
	test_cleantmp
	exec cg split -s 1 data/reg4.tsv tmp/split- .tsv.zst
	exec cg zcat tmp/split-3.tsv.zst > tmp/split-3.tsv
	list [bsort [glob tmp/*]] [file_read tmp/split-3.tsv]
} {{tmp/split-1.tsv.zst tmp/split-2.tsv.zst tmp/split-3.tsv tmp/split-3.tsv.zst tmp/split-M.tsv.zst tmp/split-X.tsv.zst tmp/split-Y.tsv.zst} {## header
##
test	chromosome	begin	end	test2
t	chr3	1000	1100	t2
}}

test tsv_split {split compressed sorted partly done} {
	test_cleantmp
	exec cg split -s 1 data/reg4.tsv tmp/split- .tsv.zst
	file delete tmp/split-2.tsv.zst
	exec cg split -s 1 data/reg4.tsv tmp/split- .tsv.zst 2> /dev/null
	exec cg zcat tmp/split-3.tsv.zst > tmp/split-3.tsv
	list [bsort [glob tmp/*]] [file_read tmp/split-3.tsv]
} {{tmp/split-1.tsv.zst tmp/split-2.tsv.zst tmp/split-3.tsv tmp/split-3.tsv.zst tmp/split-M.tsv.zst tmp/split-X.tsv.zst tmp/split-Y.tsv.zst} {## header
##
test	chromosome	begin	end	test2
t	chr3	1000	1100	t2
}}

test tsv_split {split stdin} {
	test_cleantmp
	exec cat data/reg4.tsv | cg split -- - tmp/split- .tsv
	list [bsort [glob tmp/*]] [file_read tmp/split-3.tsv]
} {{tmp/split-1.tsv tmp/split-2.tsv tmp/split-3.tsv tmp/split-M.tsv tmp/split-X.tsv tmp/split-Y.tsv} {## header
##
test	chromosome	begin	end	test2
t	chr3	1000	1100	t2
}}

test tsv_split {split and cat} {
	test_cleantmp
	exec cg split data/reg4.tsv tmp/split- .tsv
	exec cg cat -s -c f {*}[glob tmp/split-*.tsv] > tmp/cat.tsv
	exec diff tmp/cat.tsv data/reg4.tsv
} {}

test tsv_histo {cg histo basic} {
	test_cleantmp
	exec cg histo coverage | cg select -s value < data/coverage.tsv
} {value	count
0	4
1	7
10	1
11	2
12	1
20	2
22	2
30	2
31	2
39	1}

test tsv_histo {cg histo -header} {
	test_cleantmp
	exec cg histo -header 0 2 | cg select -s value < data/coverage.tsv
} {value	count
0	4
1	7
10	1
11	2
12	1
20	2
22	2
30	2
31	2
39	1
coverage	1}

test tsv_paste {basic} {
	test_cleantmp
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
		chr1 4002 4003 snp A T
	}
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample2	value-sample2
		m	2
		t	x2
	}]
	exec cg paste tmp/vars1.tsv tmp/sample1.tsv
} {#	varcomment
chromosome	begin	end	type	ref	alt	zyg-sample2	value-sample2
chr1	4001	4002	snp	A	G,C	m	2
chr1	4002	4003	snp	A	T	t	x2}

test tsv_paste {basic -o} {
	test_cleantmp
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
		chr1 4002 4003 snp A T
	}
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample2	value-sample2
		m	2
		t	x2
	}]
	exec cg paste -o tmp/result.tsv tmp/vars1.tsv tmp/sample1.tsv
	file_read tmp/result.tsv
} {#	varcomment
chromosome	begin	end	type	ref	alt	zyg-sample2	value-sample2
chr1	4001	4002	snp	A	G,C	m	2
chr1	4002	4003	snp	A	T	t	x2
}

test tsv_paste {missing column} {
	test_cleantmp
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
		chr1 4002 4003 snp A T
	}
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample2	value-sample2
		m	2
		t
	}]
	exec cg paste tmp/vars1.tsv tmp/sample1.tsv
} {#	varcomment
chromosome	begin	end	type	ref	alt	zyg-sample2	value-sample2
chr1	4001	4002	snp	A	G,C	m	2
chr1	4002	4003	snp	A	T	t	}

test tsv_paste {extra column error} {
	test_cleantmp
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
		chr1 4002 4003 snp A T
	}
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample2	value-sample2
		m	2
		t	x2	2
	}]
	exec cg paste tmp/vars1.tsv tmp/sample1.tsv > /dev/null
} {*file tmp/sample1.tsv has more columns in a line than header*} error match

test tsv_paste {multiple files} {
	test_cleantmp
	set files [makepastetest 4]
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample1	value-sample1
		m	1
	}]
	exec cg paste {*}$files > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {multiple files -o} {
	test_cleantmp
	set files [makepastetest 4]
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample1	value-sample1
		m	1
	}]
	exec cg paste -o tmp/result.tsv {*}$files
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {4 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 4]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {7 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 7]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
	bsort [glob tmp/*]
} {tmp/expected.tsv tmp/result.tsv tmp/sample1.tsv tmp/sample2.tsv tmp/sample3.tsv tmp/sample4.tsv tmp/sample5.tsv tmp/sample6.tsv tmp/vars1.tsv}

test tsv_paste {8 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 8]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {9 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 9]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {10 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 10]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {11 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 11]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {50 files > maxopenfiles -o} {
	test_cleantmp
	set files [makepastetest 50]
	exec cg paste -m 5 -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {max files == maxopenfiles -o} {
	test_cleantmp
	set maxfiles [exec sh -c {ulimit -n}]
	# if limit is too big, we get too many files to handle
	# reduce the number (but we are no longer really testing the limit)
	if {$maxfiles > 5000} {set maxfiles 5000}
	# stdin, stdout and stderr are also file descriptors (and prog)
	incr maxfiles -10
	set files [makepastetest $maxfiles]
	exec cg paste -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {max files > maxopenfiles -o} {
	test_cleantmp
	set maxfiles [exec sh -c {ulimit -n}]
	# if limit is too big, we get too many files to handle
	# reduce the number (but we are no longer really testing the limit)
	if {$maxfiles > 5000} {set maxfiles 5000}
	incr maxfile 10
	# stdin, stdout and stderr are also file descriptors
	set files [makepastetest $maxfiles]
	exec cg paste -o tmp/result.tsv {*}$files >@ stdout 2>@ stderr
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv_paste {size mismatch} {
	test_cleantmp
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
		chr1 4002 4003 snp A T
	}
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample2	value-sample2
		m	2
	}]
	exec cg paste tmp/vars1.tsv tmp/sample1.tsv > /dev/null
} {*file tmp/sample1.tsv has less lines than other files in paste*} error match

test tsv_paste {size mismatch} {
	test_cleantmp
	write_tab tmp/vars1.tsv {
		# varcomment
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
	}
	write_tab tmp/sample1.tsv [subst {
		# comment
		zyg-sample2	value-sample2
		m	2
		t	x2
	}]
	exec cg paste tmp/vars1.tsv tmp/sample1.tsv
} {*file tmp/sample1.tsv has more lines than other files in paste*} error match

test mergesorted {basic} {
	test_cleantmp
	write_deindent tmp/vars1.tsv {
		# varcomment
		chromosome	name	begin	other
		chr1	t1	4000	a
		chr1	t2	4100	a
		chr2	t3	4050	a
	}
	write_deindent tmp/vars2.tsv {
		# varcomment
		chromosome	name	begin	other
		chr1	t4	3000	b
		chr1	t4	4050	b
		chr2	t5	4100	b
	}
	write_deindent tmp/vars3.tsv {
		# varcomment
		chromosome	name	begin	other
		chr1	t6	2000	c
	}
	exec mergesorted \# 1 {} {0 2} tmp/vars1.tsv tmp/vars2.tsv tmp/vars3.tsv
} {# varcomment
chromosome	name	begin	other
chr1	t6	2000	c
chr1	t4	3000	b
chr1	t1	4000	a
chr1	t4	4050	b
chr1	t2	4100	a
chr2	t3	4050	a
chr2	t5	4100	b}

test mergesorted {no header line, comments with @} {
	test_cleantmp
	write_deindent tmp/vars1.tsv {
		@ varcomment
		@ varcomment 1
		chr1	t1	4000	a
		chr1	t2	4100	a
		chr2	t3	4050	a
	}
	write_deindent tmp/vars2.tsv {
		@ varcomment
		@ varcomment 2
		chr1	t4	3000	b
		chr1	t4	4050	b
	}
	write_deindent tmp/vars3.tsv {
		@ varcomment
		@ varcomment 3
		chr1	t6	2000	c
		chr2	t5	4100	b
	}
	exec mergesorted @ 0 {} {0 2} tmp/vars1.tsv tmp/vars2.tsv tmp/vars3.tsv
} {@ varcomment
@ varcomment 1
chr1	t6	2000	c
chr1	t4	3000	b
chr1	t1	4000	a
chr1	t4	4050	b
chr1	t2	4100	a
chr2	t3	4050	a
chr2	t5	4100	b}

test mergesorted {header} {
	write_deindent tmp/vars1.tsv {
		# varcomment
		chromosome	name	begin	other
		chr1	t1	4000	a
	}
	write_deindent tmp/vars2.tsv {
		# varcomment
		chromosome	name	begin	other
		chr1	t4	3000	b
		chr1	t4	4050	b
		chr2	t5	4100	b
	}
	write_deindent tmp/vars3.tsv {
		# varcomment
		chromosome	name	begin	other
		chr1	t6	2000	c
	}
	exec mergesorted \# 1 {chr	n	b	o} {0 2} tmp/vars1.tsv tmp/vars2.tsv tmp/vars3.tsv
} {chr	n	b	o
chr1	t6	2000	c
chr1	t4	3000	b
chr1	t1	4000	a
chr1	t4	4050	b
chr2	t5	4100	b}

test mergesorted {basic} {
	write_deindent tmp/reg1.tsv {
		chromosome	begin	end	test
		1	10	20	a1
		1	90	100	b1
		2	10	20	c1
		10	40	50	d1
	}
	write_deindent tmp/reg2.tsv {
		chromosome	begin	end	test
		1	8	9	a2
		1	40	41	b2
		5	42	43	c2
		10	51	51	d2
	}
	cg mergesorted tmp/reg1.tsv tmp/reg2.tsv
} {chromosome	begin	end	test
1	8	9	a2
1	10	20	a1
1	40	41	b2
1	90	100	b1
2	10	20	c1
5	42	43	c2
10	40	50	d1
10	51	51	d2}

test mergesorted {short} {
	write_deindent tmp/reg1.tsv {
		chromosome	begin	end	test
		1	10	20	a1
		1	90	100	b1
	}
	write_deindent tmp/reg2.tsv {
		chromosome	begin	end	test
		1	8	9	a2
	}
	cg mergesorted tmp/reg1.tsv tmp/reg2.tsv
} {chromosome	begin	end	test
1	8	9	a2
1	10	20	a1
1	90	100	b1}

test mergesorted {compressed (gz_popen test)} {
	write_deindent tmp/reg1.tsv {
		chromosome	begin	end	test
		1	10	20	a1
		1	90	100	b1
	}
	write_deindent tmp/reg2.tsv {
		chromosome	begin	end	test
		1	8	9	a2
	}
	cg zst tmp/reg1.tsv
	cg mergesorted tmp/reg1.tsv.zst tmp/reg2.tsv
} {chromosome	begin	end	test
1	8	9	a2
1	10	20	a1
1	90	100	b1}

test mergesorted {error in gz_popen test} {
	write_deindent tmp/reg1.tsv {
		chromosome	begin	end	test
		1	10	20	a1
		1	90	100	b1
	}
	file rename tmp/reg1.tsv tmp/reg1.tsv.zst
	write_deindent tmp/reg2.tsv {
		chromosome	begin	end	test
		1	8	9	a2
	}
	cg mergesorted tmp/reg1.tsv.zst tmp/reg2.tsv
} {error reading header of file tmp/reg1.tsv.zst: error closing file "*/reg1.tsv.zst": zstd-mt: */reg1.tsv.zst: Malformed input} match error

test mergesorted {error file does not exist} {
	write_deindent tmp/reg2.tsv {
		chromosome	begin	end	test
		1	8	9	a2
	}
	cg mergesorted tmp/reg1.tsv.zst tmp/reg2.tsv
} {could not read "tmp/reg1.tsv.zst": no such file or directory} error

test renamesamples {basic renamesamples test} {
	write_deindent tmp/test.tsv {
		chromosome	begin	end	test-sample1	test-sample2 test-sample3
		1	8	9	a	b	c
	}
	cg renamesamples tmp/test.tsv sample1 samplea	sample3 samplec sample4 nosample
	file_read tmp/test.tsv
} {chromosome	begin	end	test-samplea	test-sample2 test-samplec
1	8	9	a	b	c
}

test fixtsv {basic fixtsv} {
	write_deindent tmp/test.tsv {
		chromosome	begin	end	test-sample1	test-sample2	begin	test-sample3
		1	8	9	a	b	8	c
		2	8	9	a	b
		3	8	9	a	b	8	c	d
	}
	write_deindent tmp/expected.tsv {
		chromosome	begin	end	test-sample1	test-sample2	test-sample3
		1	8	9	a	b	c
		2	8	9	a	b	
		3	8	9	a	b	c
	}
	catch {cg fixtsv tmp/test.tsv tmp/test2.tsv} msg
	exec diff tmp/test2.tsv tmp/expected.tsv
	set msg
} {header shows duplicate fields (begin), removing
line 1 is of wrong length: 5 iso 7	2 8 9 a b, fixing
line 2 is of wrong length: 8 iso 7	3 8 9 a b 8 c d, fixing}

test findfields {basic findfields} {
	set header {chromosome chromStart chromEnd type reference alternative name2 transcriptid}
	findfields $header {chromosome ref type begin end alt gene geneid transcript}
} {chromosome reference type chromStart chromEnd alternative name2 name2 transcriptid}

testsummarize
