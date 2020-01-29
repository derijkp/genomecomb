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

proc checkbsort {list} {
	set pos 1
	foreach a [lrange $list 0 end-1] {
		foreach b [lrange $list $pos end] {
			# putsvars a b
			set sorted [list $a $b]
			set sort [bsort $sorted]
			# putsvars sort
			if {$sort ne $sorted} {
				error "error on [list bsort $sorted]\n -> $sort"
			}
			set sort [bsort [list $b $a]]
			# putsvars sort
			if {$sort ne $sorted} {
				error "error2 on [list bsort [list $b $a]]\n -> $sort"
			}
		}
		incr pos
	}
	# check gnusort
	set tempfile [tempfile]
	set list [list {*}$list]
	file_write $tempfile [join [list_reverse $list] \n]\n
	set temp [split [exec gnusort8 -N $tempfile] \n]
	if {$temp ne $list} {
		error "error sorting list using gnusort8: result is\n$temp\n and should be\n$list"
	}
}

test bsort {basic} {
	bsort {chr1 chr10 chr2 chrX *}
} {* chr1 chr2 chr10 chrX}

test bsort {basic 2} {
	bsort {a1 b10 b2 a2 b1 1 *}
} {* 1 a1 a2 b1 b2 b10}

test bsort {-index 1} {
	bsort -index 1 {{1 a1} {2 b10} {3 b2} {4 a2} {5 b1} {6 1} {7 *}}
} {{7 *} {6 1} {1 a1} {4 a2} {5 b1} {3 b2} {2 b10}}

test bsort {- in string} {
	set list {-10 -2 -1 0 1e-10 1e-1 1 1.0e2 1e2 120 1.0e3 1e3 {a -2} {a -1} a-1 a-2 a-10 de-2 de-10 e-2 e-10}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {- in string} {
	set list  {-10 -2 -1 0 1e-10 1e-1 1 1.0e2 1e2 120 1.0e3 1e3 {a -10} {a -2} {a -1} a-1 a-2 a-10 de-2 de-10 e-2 e-10}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}


test bsort {v-chr1-1000-5000} {
	bsort {v-chr1-10000-20000 v-chr1-1000-5000}
} {v-chr1-1000-5000 v-chr1-10000-20000}

test bsort {a-1a a-10} {
	bsort {a-1a a-10}
} {a-1a a-10}

test bsort {a-1- a-10} {
	bsort {a-1- a-10}
} {a-1- a-10}

test bsort {1 +1} {
	set list {-2 -1 1 +1 2 +2 - -a + +a a}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {various combinations} {
	set list {-10 -2 -1.2 -1.10 -1.1 -1 1 +1 1.1 1.10 1.2 2 +2 10 +10 a-1 a-2 a-10 a1 a2 a10 a12 a101 a+1 a+2 a+10}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {number separated by dashes} {
	set list {a-1-1 a-1-2 a-1-10 a-2-1 a-10-1 a-10-2 a-10-10}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {number separated by dashes} {
	set list {-10 -2 -1 1 +1 2 +2 10 +10 a-1 a-2 a-10 a1 a2 a10 a12 a101 a+1 a+2 a+10}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {all chars} {
	set list {-10 -2 -1 0 1e-10 1e-1 1 1.0e2 1e2 120 1.0e3 1e3 {a -2} {a -1} a-1 a-2 a-10 de-2 de-10 e-2 e-10}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {scientific notation and mixed} {
	set numbers {
		-1e4 -10e3 -1e3 -10e2 -1000 -2e2 -200 -101 -1e+2 -1e2 -10e1 -100
		-2e1 -20 -12 -1e+1 -1e1 -10 -2 -1.2 -1.10 -1.1 -1 -0.2 -1e-1
		-0.10 -0.1 -0.09 -0.011 -1e-2 -0.01 -0 0 
		10e-10 0.01 +0.01 1e-2 0.011 +0.011 0.09 0.1 +0.1 0.10 +0.10 1e-1 +1e-1 0.18 19e-2 +19e-2 0.2
		1 +1 1.1 +1.1 1.10 +1.10 1.2 1.20 2 +2
		10 +10 1e1 1e+1 12 19 20 +20 2e1 +2e1 100 +100 10e1 1e2 1e+2 +1e+2 101 200 +200
		2e2 1000 +1000 10e2 1e3 19e2 10e3 1e4 1e10 10e10 1e20 1e100
	}
	checkbsort $numbers
	set list $numbers
	foreach number $numbers {
		lappend list "a$number"
		lappend list "a-$number"
		lappend list "${number}a"
		lappend list "${number}e"
		lappend list "${number}e+"
		lappend list "a $number"
		lappend list "a\t$number"
		lappend list "a\t$number\tb"
		lappend list "a $number b"
		lappend list $number "$number a"
		lappend list "$number\ta"
		lappend list "a $number"
		lappend list "a\t$number"
	}
	checkbsort [bsort $list]
} {}

test bsort {special cases of simplenum} {
	set list {
		-10-2 -10 -2 -1.2.1 -1.10.1 -1.10.2 -1.1.1 -1-2 -1a -0.2 -0.2a -0.18a -0.1.1e10 
		-0.1 0.1 0.1.1e10 0.18a 0.2 0.2a 1a 1-2 1.1.1 1.10.1 1.10.2 1.2.1 10-2 
		e-0 e-0.1 e-1 e-1.1 e-1.2 e-1.10 e-2 e-2.1 e-10 e-10.1 e0 e0.1 e1 e1.1 e1.2 e1.10 
		e1e-0 e1e-0.1 e1e-1 e1e-1.1 e1e-1.2 e1e-1.10 e1e-2 e1e-2.1 e1e-10 e1e-10.1 
		e1e0 e1e0.1 e1e1 e1e1.1 e1e1.2 e1e1.10 e1e2 e1e2.1 e1e10 e1e10.1 e2 e2.1 e10 e10.1
	}
	set test [bsort [list_reverse $list]]
	set list [list {*}$list]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {how to deal with combinations of simplenum, nonum and complex num} {
	set list {
		-10 -10a -2 -2a -1 -1a -0.2 -0.10 -0.1 0 0a 0.1 0.1a 0.10 0.10a 0.2 0.2a 
		1 +1 1a +1a 2 +2 2a +2a 10 +10 10a +10a 
		-a +a +a1 +a2 +a10 +a+1 +a+2 +a+10 a a-1 a-2 a-10 a1 a2 a10 a+1 a+2 a+10
	}
	set test [bsort [list_reverse $list]]
	set list [list {*}$list]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {all chars} {
	set list [list \t \n { } * 0 1 2 3 4 5 6 7 8 9 - + ! \" \# $ % & \' ( ) , . / : \; < = > ? @ \[ \\ \] ^ _ ` \{ | \} ~ A a B b C c D d E e F f G g H h I i J j K k L l M m N n O o P p Q q R r S s T t U u V v W w X x Y y Z z]
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
}

test bsort {adding to numbers} {
	set list {-1e+1 -1e1 -10a -2a -1 -1a -1e -1e+ 1 +1 1a +1a 1a1 1a1a 1a2 1a2a 1a10 1a10a 1aa 1e 1e+ 2a +2a 10a +10a 1e1 1e+1 -a +a +a1 +a2 +a10 a a-1 a-2 a-10 a-a a1 a2 a10 aa}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {empties, blanks and *} {
	set list [list {} * 1 - + a a\t\ta a\t*\ta a\t1\ta a\t2\ta a\t10\ta a* a1 aa]
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {looking for trouble} {
	set list {-1 0.1 1 +1 -a .1 /1 /2 /10 /a +a}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

testsummarize
