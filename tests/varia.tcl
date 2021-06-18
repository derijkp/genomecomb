#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# cg select -q 'oneof($locus,1224,1226,1598520,2816907,2816911,2818387,3038185,3042129,5019241,5054443,5211957,9605719,10679396,19364044,21233191,22632054,21542413)' cgNA19240/fannotvar-cgNA19240.tsv > testvars.tsv
# exit

test zyg {t} {
	zyg v A C A C
} t

test zyg {t} {
	zyg v C A A C
} t

test zyg {m} {
	zyg v C C A C
} m

test zyg {c} {
	zyg v C G A C
} c

test zyg {c} {
	zyg v G C A C
} c

test zyg {o} {
	zyg v G A A C
} o

test zyg {r} {
	zyg v A A A C
} r

test oargs {basic} {
	oargs test {a b} {1 2}
	list $a $b
} {1 2}

test oargs {missing arg} {
	oargs test {a b} {1}
	list $a $b
} {missing arg(s): b, should be: test a b} error

test oargs {extra arg} {
	oargs test {a b} {1 2 3}
	list $a $b
} {too many args: 3, should be: test a b} error

test oargs {default not used} {
	oargs test {a {b 5}} {1 2}
	list $a $b
} {1 2}

test oargs {default used} {
	oargs test {a {b 5}} {1}
	list $a $b
} {1 5}

test oargs {option} {
	oargs test {a {b 5}} {1 -b 2}
	list $a $b
} {1 2}

test oargs {option reorder} {
	oargs test {a {b 5}} {-b 2 1}
	list $a $b
} {1 2}

test oargs {option --} {
	oargs test {a {b 5} {c 10}} {-b 2 -- -1 -2}
	list $a $b $c
} {-1 2 -2}

test oargs {wrong option} {
	oargs test {a {b 5}} {1 2 -d 2}
	list $a $b $c
} {unknown option -d, cmd should be: test a ?b?} error

test oargs {args} {
	oargs test {a b args} {1 2 3 4}
	list $a $b $args
} {1 2 {3 4}}

test oargs {empty values} {
	oargs test {a b {c {}} {d {}}} {1 2 {} 4}
	list $a $b $c $d
} {1 2 {} 4}

test oargs {error end option} {
	oargs test {a {b {}}} {1 -b}
	list $a $b
} {option -b without value, should be: test a ?b?} error

test oargs {args option} {
	oargs test {a {b {}} args} {1 -b 2}
	list $a $b
} {1 2}

test oargs {args option} {
	oargs test {a {b {}} args} {1 -b 2 -c 3}
	list $a $b $args
} {1 2 {-c 3}}

test cg_options {basic} {
	set args {-opt o 1 2 3 4}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2}
	list $opt $p1 $p2 $args
} {o 1 2 {3 4}}

test cg_options {args error (not enough)} {
	set args {-opt o}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2}
} {
ERROR: Wrong number of arguments, correct format is:
cg test ?options? p1 p2 ...
  with options: -opt} error

test cg_options {args option} {
	set args {-opt o 1 2 3 4}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2} 0 4
	list $opt $p1 $p2 $args
} {o 1 2 {3 4}}

test cg_options {option error} {
	set args {-unknown o 1 2 3 4}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2} 0 4
} {unknown option "-unknown", must be one of: -opt} error

test cg_options {args error (not enough)} {
	set args {-opt o}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2} 1 4
} {
ERROR: Wrong number of arguments, correct format is:
cg test ?options? p1 ?p2? ...
  with options: -opt} error

test cg_options {args error (too much)} {
	set args {-opt o 1 2 3 4}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2} 1 2
} {
ERROR: Wrong number of arguments, correct format is:
cg test ?options? p1 ?p2?
  with options: -opt} error

test cg_options {args option} {
	set args {-opt o 1 2 3 4}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2 p3 p4} 0 4
	list $opt $p1 $p2 $p3 $p4 $args
} {o 1 2 3 4 {}}

test cg_options {args no parameters given} {
	set args {-opt o 1 2}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {} 1
	set args
} {1 2}

test cg_options {do not set parameters not given} {
	set args {-opt o 1}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {a b} 1 2
	list $a [info exists b]
} {1 0}

test tsvdiff {basic} {
	write_tab tmp/file1.tsv {
		chromosome	begin	end
		chr1	1	2
		chr2	1	2
	}
	cg zst tmp/file1.tsv
	write_tab tmp/file2.tsv {
		chromosome	begin	end	test
		chr1	1	2	t1
		chr2	2	3	t2
		chr3	3	4	t3
	}
	cg tsvdiff tmp/file1.tsv.zst tmp/file2.tsv
} {diff tmp/file1.tsv.zst tmp/file2.tsv
header diff
<extrafields: 
---
>extrafields: test
header
  chromosome	begin	end
3c3,4
< chr2	1	2
---
> chr2	2	3
> chr3	3	4
child process exited abnormally} error

test tsvdiff {dir} {
	file mkdir tmp/d1
	file mkdir tmp/d2
	write_tab tmp/d1/file1.tsv {
		chromosome	begin	end
		chr1	1	2
		chr2	1	2
	}
	file copy -force tmp/d1/file1.tsv tmp/d1/same.tsv
	file_write tmp/d1/only1 ""
	# cg zst tmp/d1/file1.tsv
	write_tab tmp/d2/file1.tsv {
		chromosome	begin	end	test
		chr1	1	2	t1
		chr2	2	3	t2
		chr3	3	4	t3
	}
	file copy -force tmp/d1/file1.tsv tmp/d2/same.tsv
	file_write tmp/d2/only2 ""
	cg tsvdiff tmp/d1 tmp/d2
} {diff tmp/d1/file1.tsv tmp/d2/file1.tsv
header diff
<extrafields: 
---
>extrafields: test
header
  chromosome	begin	end
3c3,4
< chr2	1	2
---
> chr2	2	3
> chr3	3	4
Only in 1: tmp/d1/only1
Only in 2: tmp/d2/only2
child process exited abnormally} error

test tsvdiff {dir with sam vs bam} {
	test_cleantmp
	file mkdir tmp/d1
	file mkdir tmp/d2
	#
	file copy data/bwa.sam tmp/d1/same.sam
	exec samtools view -b data/bwa.sam > tmp/d2/same.bam
	exec samtools view -b data/bwa.sam > tmp/d1/same2.bam
	exec samtools view -b data/bwa.sam > tmp/d2/same2.bam
	exec gzip -c data/bwa.sam > tmp/d1/same3.sam.gz
	exec samtools view -b data/bwa.sam > tmp/d2/same3.bam
	exec samtools view -b data/bwa.sam > tmp/d1/same4.bam
	exec gzip -c data/bwa.sam > tmp/d2/same4.sam.gz
	#
	file copy data/bwa.sam tmp/d1/diff.sam
	exec samtools view -b data/rdsbwa.sam > tmp/d2/diff.bam
	exec samtools view -b data/bwa.sam > tmp/d1/diff2.bam
	exec samtools view -b data/rdsbwa.sam > tmp/d2/diff2.bam
	exec gzip -c data/bwa.sam > tmp/d1/diff3.sam.gz
	exec samtools view -b data/rdsbwa.sam > tmp/d2/diff3.bam
	exec samtools view -b data/rdsbwa.sam > tmp/d1/diff4.bam
	exec gzip -c data/bwa.sam > tmp/d2/diff4.sam.gz
	#
	write_tab tmp/d1/file1.tsv {
		chromosome	begin	end
		chr1	1	2
		chr2	1	2
	}
	file_write tmp/d1/only1 ""
	# cg zst tmp/d1/file1.tsv
	write_tab tmp/d2/file1.tsv {
		chromosome	begin	end	test
		chr1	1	2	t1
		chr2	2	3	t2
		chr3	3	4	t3
	}
	file copy -force tmp/d1/file1.tsv tmp/d1/same.tsv
	file copy -force tmp/d1/file1.tsv tmp/d2/same.tsv
	file_write tmp/d2/only2 ""
	cg tsvdiff -brief 1 tmp/d1 tmp/d2
} {Files differ: tmp/d1/diff.sam tmp/d2/diff.bam
Files differ: tmp/d1/diff2.bam tmp/d2/diff2.bam
Files differ: tmp/d1/diff3.sam.gz tmp/d2/diff3.bam
Files differ: tmp/d1/diff4.bam tmp/d2/diff4.sam.gz
Files differ: tmp/d1/file1.tsv tmp/d2/file1.tsv
Only in 1: tmp/d1/only1
Only in 2: tmp/d2/only2
child process exited abnormally} error

test cg_extracthomopolymers {basic} {
	file_write tmp/genome_test.ifas [deindent {
		>chr1
		GCCGAAAAAAAAAGCC
		>chr1_test
		GCGAGGGGGGGCTGTGCAAAAAAAA
	}]
	cg extracthomopolymers tmp/genome_test.ifas
} {chromosome	begin	end	base	size
1	4	13	A	9
1_test	4	11	G	7
1_test	17	25	A	8}

test link {find_link} {
	file_write tmp/1.txt 1
	mklink tmp/1.txt tmp/2.txt
	mklink tmp/2.txt tmp/3.txt
	find_link tmp/3.txt
} {*tmp/1.txt} match

test link {find_link} {
	file_write tmp/1.txt 1
	mklink tmp/1.txt tmp/2.txt
	mklink tmp/2.txt tmp/3.txt
	find_link tmp/3.txt 1
} {*tmp/2.txt} match

test error {catch exec testerror 1 1 1} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	exec data/testerror.sh 1 1 1
} {written to stdout
written to stderr} error

test error {catch exec testerror 1 1 0 (catch sees error if output to stderr, even if no error exit)} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	exec data/testerror.sh 1 1 0
} {written to stdout
written to stderr} error

test error {catch exec testerror 0 0 1 (allways error on error exit} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	exec data/testerror.sh 0 0 1
} {child process exited abnormally} error

test error {catch_exec testerror 1 1 1} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	catch_exec data/testerror.sh 1 1 1
} {written to stdout
written to stderr} error

test error {catch_exec testerror 1 1 0 (catch_exec only sees error if exit error} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	catch_exec data/testerror.sh 1 1 0
} {written to stdout
written to stderr}

test error {catch_exec testerror 0 0 1 (allways error on error exit)} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	catch_exec data/testerror.sh 0 0 1
} {child process exited abnormally} error

test error {catch catch_exec testerror 1 1 1} {
	# args to data/testerror.sh: writestdout writestderr exit_with_error
	set error [catch {
		catch_exec data/testerror.sh 1 1 1
	} msg]
} 1

test analysisinfo_write {analysisinfo_write} {
	file_write tmp/src.analysisinfo [deindent {
		source
		test
	}]\n
	analysisinfo_write tmp/src tmp/dest dest test2
	file_read tmp/dest.analysisinfo
} {source	dest
test	test2
}

test analysisinfo_write {analysisinfo_write no duplicates} {
	file_write tmp/src.analysisinfo [deindent {
		source	dest
		test	test1
	}]\n
	analysisinfo_write tmp/src tmp/dest dest testnew dest2 test2
	file_read tmp/dest.analysisinfo
} {source	dest	dest2
test	test1	test2
}

test analysisinfo_write {analysisinfo_write no src.analysisinfo} {
	analysisinfo_write tmp/src tmp/dest dest test2
	file_read tmp/dest.analysisinfo
} {dest
test2
}

test analysisinfo_write {analysisinfo_write empty src.analysisinfo} {
	file_write tmp/src.analysisinfo {}
	analysisinfo_write tmp/src tmp/dest dest test2
	file_read tmp/dest.analysisinfo
} {dest
test2
}

test analysisinfo_write {analysisinfo_write add} {
	file_write tmp/src.analysisinfo [deindent {
		annotate_cg_version
		0.101.0
	}]\n
	analysisinfo_write tmp/src tmp/dest.zst sample cg-cg-tmp
	file_read tmp/dest.analysisinfo
} {sample	annotate_cg_version
cg-cg-tmp	0.101.0
}

test analysisinfo_write {analysisinfo_write from empty} {
	analysisinfo_write {} tmp/dest dest test2
	file_read tmp/dest.analysisinfo
} {dest
test2
}

testsummarize
