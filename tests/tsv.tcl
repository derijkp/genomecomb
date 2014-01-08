#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test tsv_open {sft} {
	set f [open data/reg1.tsv]
	set header [tsv_open $f keepheader]
	close $f
	list $header $keepheader
} {{chromosome test begin end} {}}

#test tsv_open {cg} {
#	set f [open data/]
#	set header [tsv_open $f keepheader]
#	close $f
#	list $header $keepheader
#} {{chromosome test begin end} {}}

test tsv_open {rtg} {
	set f [open data/rtgsnps.tsv]
	set header [tsv_open $f keepheader]
	close $f
	tsv_hcheader $f keepheader header
	list $header $keepheader
} {{name position type reference prediction posterior coverage correction support_statistics} {#Version v2.0-RC build 27737 (2010-05-06), SNP output v2.0
#CL	snp -t /rtgshare/data/human/sdf/hg18 -o snps_GS00102-DNA-D06_m4_u4 --max-as-mated=4 --max-as-unmated=4 --max-ih=1 -m cg_errors --no-complex-calls /rtgshare3/users/richard/vib_20100324_GS00102-DNA-D06-hg18ref/map_GS000004945/mated.sam.gz /rtgshare3/users/richard/vib_20100324_GS00102-DNA-D06-hg18ref/map_GS000004945/unmated.sam.gz ...
}}

test tsv_open {vcf} {
	catch {close $f}
	set f [open data/test.vcf]
	set header [tsv_open $f keepheader]
	close $f
	list $header $keepheader
} {{CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003} {##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#}}

test tsv_cat {one} {
	cg cat data/reg1b.tsv
} {chromosome	test	begin	end
1	t	10	20
1	t	50	60}

test tsv_cat {two same header} {
	cg cat data/reg2.tsv data/reg4.tsv
} {# ++++ data/reg2.tsv ++++
# comments added
#
# ++++ data/reg4.tsv ++++
## vcf style header
##
#
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
t	chrY	1010	1900	t2}

test tsv_cat {two diff header} {
	exec cg cat data/reg1b.tsv data/reg2.tsv
} {headers do not match, use -f to force or -m to merge} error

test tsv_cat {two diff header -m} {
	cg cat -m data/reg2.tsv data/reg1b.tsv
} {# ++++ data/reg2.tsv ++++
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
t	1	50	60	}

test tsv_cat {two diff header -m} {
	cg cat -m data/reg1.tsv data/reg2.tsv
} {# ++++ data/reg1.tsv ++++
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
Y	t	1010	1900	t2}

test tsv_cat {two diff header -f} {
	cg cat -m data/reg2.tsv data/reg1b.tsv
} {# ++++ data/reg2.tsv ++++
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
t	1	50	60	}

set ::env(PATH) $keeppath

test check_sort {sort error 1 in vars} {
	exec cg checksort data/vars_sorterror1.sft
} {error in file data/vars_sorterror1.sft: file is not correctly sorted (sort correctly using "cg select -s -")
chr10:43198434-43198435:del: came before chr3:52847042-52847060:del:} error

test check_sort {sort error 2 in vars} {
	exec cg checksort data/vars_sorterror2.sft
} {error in file data/vars_sorterror2.sft: file is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:del: came before chr3:52847042-52847060:del:} error

test check_sort {sort error 3 in vars} {
	exec cg checksort data/vars_sorterror3.sft
} {error in file data/vars_sorterror3.sft: file is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847310:snp:G came before chr3:52847303-52847304:snp:G} error

test check_sort {sort error 4 in vars} {
	exec cg checksort data/vars_sorterror4.sft
} {error in file data/vars_sorterror4.sft: file is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:ins:G came before chr3:52847303-52847304:ins:G} error

test check_sort {sort error 5 in vars} {
	cg checksort data/vars_sorterror5.sft
} {error in file data/vars_sorterror5.sft: file is not correctly sorted (sort correctly using "cg select -s -")*} error match

testsummarize
