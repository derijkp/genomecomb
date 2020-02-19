#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc misc_testfiles {} {
	file mkdir tmp
	file mkdir tmp/test
	file_write tmp/test/test1 1
	file_write tmp/test/test2 2
	file mkdir tmp/test/subdir
	file_write tmp/test/subdir/test3 3
	file mkdir tmp/out
}

proc dirinfo {dir} {
	set result {}
	foreach file [bsort [glob -nocomplain $dir/*]] {
		if {[catch {file link $file} link]} {set link -}
		if {[file isdir $file]} {
			set data [list dir [dirinfo $file]]
		} else {
			if {[file exists $file]} {
				set data [file_read $file]
			} else {
				set data __file_does_not_exist__
			}
		}
		lappend result [list $file $link $data]
	}
	return $result
}

test cplinked {relative} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {over existing} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {relative, new dir} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out/test2
	set result [dirinfo tmp/out/test2]
	file delete -force tmp/out/test2
	set result
} {{tmp/out/test2/subdir - {dir {{tmp/out/test2/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test2/test1 ../../test/test1 1} {tmp/out/test2/test2 ../../test/test2 2}} 

test cplinked {relative, new dir 2} {
	misc_testfiles
	file delete -force tmp/out
	exec cg cplinked tmp/test tmp/out/test
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {relative extra level} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out/test2
	exec cg cplinked tmp/test tmp/out/test2/testl
	set result [dirinfo tmp/out/test2/testl]
	file delete -force tmp/out
	set result
} {{tmp/out/test2/testl/subdir - {dir {{tmp/out/test2/testl/subdir/test3 ../../../../test/subdir/test3 3}}}} {tmp/out/test2/testl/test1 ../../../test/test1 1} {tmp/out/test2/testl/test2 ../../../test/test2 2}} 

set dir [testdir cplinked {absolute}]
test cplinked {absolute} {
	misc_testfiles
	file delete -force /tmp/test2
	exec cg cplinked tmp/test /tmp/test2
	set result [dirinfo /tmp/test2]
	file delete -force /tmp/test2
	set result
} [list [list /tmp/test2/subdir - [list dir [list [list /tmp/test2/subdir/test3 $dir/tmp/test/subdir/test3 3]]]] [list /tmp/test2/test1 $dir/tmp/test/test1 1] [list /tmp/test2/test2 $dir/tmp/test/test2 2]]

test cplinked {backup existing} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out/test
	file_write tmp/out/test/test1 pre
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test1.old1 - pre} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {backup existing 2} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out/test
	file_write tmp/out/test/test1 pre
	file_write tmp/out/test/test1.old1 old
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test1.old1 - old} {tmp/out/test/test1.old2 - pre} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {file} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test/test1 tmp/out/testl
	set result [dirinfo tmp/out/]
	file delete -force tmp/out/test
	set result
} {{tmp/out/testl ../test/test1 1}} 

test cplinked {file replace link} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test/test1 tmp/out/testl
	exec cg cplinked tmp/test/test2 tmp/out/testl
	set result [dirinfo tmp/out/]
	file delete -force tmp/out/test
	set result
} {{tmp/out/testl ../test/test2 2}} 

test cplinked {file replace file} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	file_write tmp/out/testl pre
	exec cg cplinked tmp/test/test1 tmp/out/testl
	set result [dirinfo tmp/out/]
	file delete -force tmp/out/test
	set result
} {{tmp/out/testl ../test/test1 1} {tmp/out/testl.old1 - pre}} 

test cplinked {dangling link} {
	misc_testfiles
	file delete -force tmp/out
	file mkdir tmp/out
	file_write tmp/out/test1 pre
	mklink tmp/out/remove tmp/out/testlink
	file delete tmp/outresult
	exec cg cplinked tmp/out tmp/outresult
	set result [dirinfo tmp/outresult/]
	file delete -force tmp/out tmp/outresult
	set result
} {{tmp/outresult/test1 ../out/test1 pre} {tmp/outresult/testlink ../out/testlink __file_does_not_exist__}}

test mklink {nklink time} {
	file_write tmp/test.txt ""
	exec touch -d "2018-01-01 12:00" tmp/test.txt
	mklink tmp/test.txt tmp/link.txt
	expr {[job_file_mtime tmp/test.txt] == [job_file_mtime tmp/link.txt]}
} 1

test distr2chr {basic} {
	test_cleantmp
	exec distr2chr tmp/distrvars1- 0 < data/vars1.sft
	list [bsort [glob tmp/*]] [file_read tmp/distrvars1-chr1]
} {{tmp/distrvars1-chr1 tmp/distrvars1-chr2 tmp/distrvars1-chromosome} {chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
}} 

test distr2chr {many} {
	test_cleantmp
	exec distr2chr tmp/distrvars1- 0 < data/manychr.tsv
	exec cat {*}[glob tmp/distrvars1-*] > tmp/cat.tsv
	exec sort tmp/cat.tsv > tmp/scat.tsv
	exec sort data/manychr.tsv > tmp/check.tsv
	set result {}
	lappend result [exec diff tmp/scat.tsv tmp/check.tsv]
	set files [glob tmp/distrvars1-*]
	foreach file $files {
		set temp [lindex [split $file -] end]
		catch {exec grep -v $temp $file} r
		if {$r ne "child process exited abnormally"} {error "file $file error: $r"}
	}
	lappend result [llength $files]
	lappend result [file_read tmp/distrvars1-chr2]
	set result
} {{} 81 {chr2	2
chr2	2-2
}} 

test distr2chr {reglist} {
	test_cleantmp
	exec distr2chr tmp/distrvars1- 0 < data/vars1.sft
	list [bsort [glob tmp/*]] [file_read tmp/distrvars1-chr1]
} {{tmp/distrvars1-chr1 tmp/distrvars1-chr2 tmp/distrvars1-chromosome} {chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
}} 

test distrreg {basic} {
	test_cleantmp
	exec distrreg tmp/distrvars1- {} 1 \
		{chr1-1000-5000 chr1-5000-10000 chr1-10000-20000 chr2} \
		0 1 2 1 \# < data/vars1.sft
	list [bsort [glob tmp/*]] \
		[file_read tmp/distrvars1-chr1-1000-5000] \
		[file_read tmp/distrvars1-chr1-5000-10000] \
		[lindex [cg select -g all tmp/distrvars1-chr1-10000-20000] end] \
		[lindex [cg select -g all tmp/distrvars1-chr2] end]
} {{tmp/distrvars1-chr1-1000-5000 tmp/distrvars1-chr1-5000-10000 tmp/distrvars1-chr1-10000-20000 tmp/distrvars1-chr2} {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2
} {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
} 0 9} 

test distrreg {basic with compression} {
	test_cleantmp
	set temp [file_read data/vars1.sft]
	file_write tmp/vars1.sft \#test\ comment\n$temp
	exec distrreg tmp/distrvars1- .zst 1 \
		{chr1-1000-5000 chr1-5000-10000 chr1-10000-20000 chr2} \
		0 1 2 1 \# {zstd-mt -k -q -8 -b 0.5 -T 1 -c > } < tmp/vars1.sft
	list [bsort [glob tmp/*]] \
		[exec cg zcat tmp/distrvars1-chr1-1000-5000.zst] \
		[exec cg zcat tmp/distrvars1-chr1-5000-10000.zst] \
		[lindex [cg select -g all tmp/distrvars1-chr1-10000-20000.zst] end] \
		[lindex [cg select -g all tmp/distrvars1-chr2.zst] end]
} {{tmp/distrvars1-chr1-1000-5000.zst tmp/distrvars1-chr1-5000-10000.zst tmp/distrvars1-chr1-10000-20000.zst tmp/distrvars1-chr2.zst tmp/vars1.sft} {#test comment
chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2} {#test comment
chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3} 0 9} 

test distrreg {chr1_} {
	test_cleantmp
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	name
		chr1	2	3	c1-2
		chr1	80	81	c1-80
		chr1	180	181	c1-180
		chr1_random1	90	91	c11-90
		chr1_random2	99	100	c12-99
		chr2	1000	1001	c2-1000
	}]]\n
	exec distrreg tmp/distrvars1- {} 1 \
		{chr1-1-100 chr1-100-200 chr1_ chr2} \
		0 1 2 1 \# < tmp/test.tsv
	file delete tmp/test.tsv
	set files [bsort [glob tmp/*]]
	set result $files\n
	foreach file $files {
		append result \#$file\n[file_read $file]
	}
	set result
} {tmp/distrvars1-chr1-1-100 tmp/distrvars1-chr1-100-200 tmp/distrvars1-chr1_ tmp/distrvars1-chr2
#tmp/distrvars1-chr1-1-100
chromosome	begin	end	name
chr1	2	3	c1-2
chr1	80	81	c1-80
#tmp/distrvars1-chr1-100-200
chromosome	begin	end	name
chr1	180	181	c1-180
#tmp/distrvars1-chr1_
chromosome	begin	end	name
chr1_random1	90	91	c11-90
chr1_random2	99	100	c12-99
#tmp/distrvars1-chr2
chromosome	begin	end	name
chr2	1000	1001	c2-1000
}

test distrreg {chr1_ and unaligned} {
	test_cleantmp
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	name
		chr1	2	3	c1-2
		chr1	80	81	c1-80
		chr1	180	181	c1-180
		chr1_random1	90	91	c11-90
		chr1_random2	99	100	c12-99
		chr2	1000	1001	c2-1000
		*	0	0	*-0
	}]]\n
	exec distrreg tmp/distrvars1- {} 1 \
		{chr1-1-100 chr1-100-200 chr1_ chr2 unaligned} \
		0 1 2 1 \# < tmp/test.tsv
	file delete tmp/test.tsv
	set files [bsort [glob tmp/*]]
	set result $files\n
	foreach file $files {
		append result \#$file\n[file_read $file]
	}
	set result
} {tmp/distrvars1-chr1-1-100 tmp/distrvars1-chr1-100-200 tmp/distrvars1-chr1_ tmp/distrvars1-chr2 tmp/distrvars1-unaligned
#tmp/distrvars1-chr1-1-100
chromosome	begin	end	name
chr1	2	3	c1-2
chr1	80	81	c1-80
#tmp/distrvars1-chr1-100-200
chromosome	begin	end	name
chr1	180	181	c1-180
#tmp/distrvars1-chr1_
chromosome	begin	end	name
chr1_random1	90	91	c11-90
chr1_random2	99	100	c12-99
#tmp/distrvars1-chr2
chromosome	begin	end	name
chr2	1000	1001	c2-1000
#tmp/distrvars1-unaligned
chromosome	begin	end	name
*	0	0	*-0
}

test distrreg {no header} {
	test_cleantmp
	exec distrreg tmp/distrvars1- .bed 0 \
		{chr1-1000-5000 chr1-5000-10000 chr1-10000-20000 chr2} \
		0 1 2 1 \# < data/vars1.sft
	list [bsort [glob tmp/*]] \
		[file_read tmp/distrvars1-chr1-1000-5000.bed] \
		[file_read tmp/distrvars1-chr1-5000-10000.bed] \
		[lindex [exec wc -l tmp/distrvars1-chr1-10000-20000.bed] 0] \
		[lindex [exec wc -l  tmp/distrvars1-chr2.bed] 0]
} {{tmp/distrvars1-chr1-1000-5000.bed tmp/distrvars1-chr1-5000-10000.bed tmp/distrvars1-chr1-10000-20000.bed tmp/distrvars1-chr2.bed} {chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2
} {chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
} 0 9} 

test distrreg {cg distrreg} {
	test_cleantmp
	file copy data/vars1.sft tmp/vars1.tsv
	cg distrreg -stack 1 tmp/vars1.tsv tmp/distrvars1- .tsv {chr1-1000-5000 chr1-5000-10000 chr1-10000-20000 chr2}
	list [bsort [glob tmp/*]] \
		[file_read tmp/distrvars1-chr1-1000-5000.tsv] \
		[file_read tmp/distrvars1-chr1-5000-10000.tsv] \
		[lindex [cg select -g all tmp/distrvars1-chr1-10000-20000.tsv] end] \
		[lindex [cg select -g all tmp/distrvars1-chr2.tsv] end]
} {{tmp/distrvars1-chr1-1000-5000.tsv tmp/distrvars1-chr1-5000-10000.tsv tmp/distrvars1-chr1-10000-20000.tsv tmp/distrvars1-chr2.tsv tmp/vars1.tsv} {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2
} {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
} 0 9} 

test distrreg {distrreg_job} {
	test_cleantmp
	file copy data/vars1.sft tmp/vars1.tsv
	job_init
	set results [distrreg_job tmp/vars1.tsv tmp/distrvars1- .tsv {chr1-1000-5000 chr1-5000-10000 chr1-10000-20000 chr2}]
	job_wait
	list [bsort $results] \
		[file_read tmp/distrvars1-chr1-1000-5000.tsv] \
		[file_read tmp/distrvars1-chr1-5000-10000.tsv] \
		[lindex [cg select -g all tmp/distrvars1-chr1-10000-20000.tsv] end] \
		[lindex [cg select -g all tmp/distrvars1-chr2.tsv] end]
} {{tmp/distrvars1-chr1-1000-5000.tsv tmp/distrvars1-chr1-5000-10000.tsv tmp/distrvars1-chr1-10000-20000.tsv tmp/distrvars1-chr2.tsv} {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	4100	snp	C	T	T	T	47	v	T	T	35	v	1,2
} {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	sequenced-sample2	list
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v			41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
} 0 9} 

test distrreg {cg distrreg chr1 vs 1} {
	test_cleantmp
	file_write tmp/vars1.tsv [deindent {
		chromosome	begin	end	name
		1	10	20	a
		1_random	20	30	b
		1_random2	20	30	b
		2	30	40	c
	}]\n
	cg distrreg -stack 1 tmp/vars1.tsv tmp/distrvars1- .tsv {chr1 chr1_ chr2}
	list [bsort [glob tmp/*]] \
		[file_read tmp/distrvars1-chr1.tsv] \
		[file_read tmp/distrvars1-chr1_.tsv]
} {{tmp/distrvars1-chr1.tsv tmp/distrvars1-chr1_.tsv tmp/distrvars1-chr2.tsv tmp/vars1.tsv} {chromosome	begin	end	name
1	10	20	a
} {chromosome	begin	end	name
1_random	20	30	b
1_random2	20	30	b
}}

test distrreg {cg distrreg 1 vs chr1} {
	test_cleantmp
	file_write tmp/vars1.tsv [deindent {
		chromosome	begin	end	name
		chr1	10	20	a
		chr1_random	20	30	b
		chr1_random2	20	30	b
		chr2	30	40	c
	}]\n
	cg distrreg -stack 1 tmp/vars1.tsv tmp/distrvars1- .tsv {1 1_ 2}
	list [bsort [glob tmp/*]] \
		[file_read tmp/distrvars1-1.tsv] \
		[file_read tmp/distrvars1-1_.tsv]
} {{tmp/distrvars1-1.tsv tmp/distrvars1-1_.tsv tmp/distrvars1-2.tsv tmp/vars1.tsv} {chromosome	begin	end	name
chr1	10	20	a
} {chromosome	begin	end	name
chr1_random	20	30	b
chr1_random2	20	30	b
}}

test file_absolute {relative path} {
	cd /tmp
	file_absolute abc
} /tmp/abc

test file_absolute {relative path with ..} {
	cd /usr/bin
	file_absolute ../abc
} /usr/abc

test file_absolute {absolute with .. and .} {
	file_absolute /abc/def/../ghi/./j
} /abc/ghi/j

test file_absolute {empty} {
	file_absolute /abc//def
} /abc/def

test file_absolute {error} {
	file_absolute /abc/def/../../../ghi/./j
} {file_absolute error: cannot .. past root} error

test distrreg_regs {chr} {
	distrreg_regs chr $::refseqdir/hg19
} {chr1 chr1_ chr2 chr3 chr4 chr4_ chr5 chr6 chr7 chr7_ chr8 chr8_ chr9 chr9_ chr10 chr11 chr11_ chr12 chr13 chr14 chr15 chr16 chr17 chr17_ chr18 chr18_ chr19 chr19_ chr20 chr21 chr21_ chr22 chrM chrUn_ chrX chrY unaligned}

test distrreg_regs {s100000000} {
	distrreg_regs s100000000 $::refseqdir/hg19
} {chr1-0-29953082 chr1-29953082-142756022 chr1-142756022-235217211 chr1-235217211-249250621 chr1_ chr2-0-93826171 chr2-93826171-149740582 chr2-149740582-243199373 chr3-0-92004854 chr3-92004854-194044606 chr3-194044606-198022430 chr4-0-75439829 chr4-75439829-191154276 chr4_ chr5-0-91661128 chr5-91661128-180915260 chr6-0-95755543 chr6-95755543-171115067 chr7-0-74740724 chr7-74740724-159138663 chr7_ chr8-0-86651451 chr8-86651451-146364022 chr8_ chr9-0-92603796 chr9-92603796-141213431 chr9_ chr10-0-51423845 chr10-51423845-135534747 chr11-0-96362584 chr11-96362584-135006516 chr11_ chr12-0-36356694 chr12-36356694-133851895 chr13-0-86835324 chr13-86835324-115169878 chr14-0-107349540 chr15-0-85009473 chr15-85009473-102531392 chr16-0-90354753 chr17-0-81195210 chr17_ chr18-0-78077248 chr18_ chr19-0-59128983 chr19_ chr20-0-63025520 chr21-0-48129895 chr21_ chr22-0-51304566 chrM-0-16571 chrUn_ chrX-0-76678692 chrX-76678692-155270560 chrY-0-59363566 unaligned}

test distrreg_regs {100000000} {
	distrreg_regs 100000000 $::refseqdir/hg19
} {chr1-0-100000000 chr1-100000000-200000000 chr1-200000000-249250621 chr1_ chr2-0-100000000 chr2-100000000-200000000 chr2-200000000-243199373 chr3-0-100000000 chr3-100000000-198022430 chr4-0-100000000 chr4-100000000-191154276 chr4_ chr5-0-100000000 chr5-100000000-180915260 chr6-0-100000000 chr6-100000000-171115067 chr7-0-100000000 chr7-100000000-159138663 chr7_ chr8-0-100000000 chr8-100000000-146364022 chr8_ chr9-0-100000000 chr9-100000000-141213431 chr9_ chr10-0-100000000 chr10-100000000-135534747 chr11-0-100000000 chr11-100000000-135006516 chr11_ chr12-0-100000000 chr12-100000000-133851895 chr13-0-100000000 chr13-100000000-115169878 chr14-0-100000000 chr14-100000000-107349540 chr15-0-100000000 chr15-100000000-102531392 chr16-0-90354753 chr17-0-81195210 chr17_ chr18-0-78077248 chr18_ chr19-0-59128983 chr19_ chr20-0-63025520 chr21-0-48129895 chr21_ chr22-0-51304566 chrM-0-16571 chrUn_ chrX-0-100000000 chrX-100000000-155270560 chrY-0-59373566 unaligned}

test gnusort {basic -N} {
	file_write tmp/text "chr10\nchr9\nchr8\nchr2\nchr1\n"
	exec gnusort8 -N tmp/text
} {chr1
chr2
chr8
chr9
chr10}

test gnusort {-N --header-lines 1} {
	file_write tmp/text "chr10\nchr9\nchr8\nchr2\nchr1\n"
	exec gnusort8 -N --header-lines 1 < tmp/text
} {chr10
chr1
chr2
chr8
chr9}

test gnusort {-N --header-lines 2} {
	file_write tmp/text "chr10\nchr9\nchr8\nchr2\nchr1\n"
	exec gnusort8 -N --header-lines 2 < tmp/text
} {chr10
chr9
chr1
chr2
chr8}

testsummarize
