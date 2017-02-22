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

test compression {lz4} {
	test_cleantmp
	file_write tmp/test1.txt a
	file_write tmp/test2.txt b
	cg lz4 {*}[glob tmp/test*.txt]
	catch {exec lz4c -d tmp/test1.txt.lz4 2> /dev/null} c1
	catch {exec lz4c -d tmp/test2.txt.lz4 2> /dev/null} c2
	list [lsort -dict [glob tmp/test*]] $c1 $c2
} {{tmp/test1.txt.lz4 tmp/test2.txt.lz4} a b}

test compression {lz4 -o} {
	test_cleantmp
	file_write tmp/test1.txt a
	cg lz4 -o tmp/out.lz4 tmp/test1.txt
	catch {exec lz4c -d tmp/out.lz4 2> /dev/null} c1
	list [lsort -dict [glob tmp/out*]] $c1
} {tmp/out.lz4 a}

test compression {lz4 -i} {
	test_cleantmp
	file_write tmp/test1.txt a
	file_write tmp/test2.txt b
	cg lz4 -i 1 {*}[glob tmp/test*.txt]
	catch {exec lz4c -d tmp/test1.txt.lz4 2> /dev/null} c1
	catch {exec lz4c -d tmp/test2.txt.lz4 2> /dev/null} c2
	list [lsort -dict [glob tmp/test*]] $c1 $c2
} {{tmp/test1.txt.lz4 tmp/test1.txt.lz4.lz4i tmp/test2.txt.lz4 tmp/test2.txt.lz4.lz4i} a b}

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
cg test ?options? p1 p2 ...} error

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
cg test ?options? p1 ?p2? ...} error

test cg_options {args error (too much)} {
	set args {-opt o 1 2 3 4}
	cg_options test args {
		-opt {
			set opt $value
		}
	} {p1 p2} 1 2
} {
ERROR: Wrong number of arguments, correct format is:
cg test ?options? p1 ?p2?} error

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
	cg lz4 tmp/file1.tsv
	write_tab tmp/file2.tsv {
		chromosome	begin	end	test
		chr1	1	2	t1
		chr2	2	3	t2
		chr3	3	4	t3
	}
	cg tsvdiff tmp/file1.tsv.lz4 tmp/file2.tsv
} {diff tmp/file1.tsv.lz4 tmp/file2.tsv
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
	# cg lz4 tmp/d1/file1.tsv
	write_tab tmp/d2/file1.tsv {
		chromosome	begin	end	test
		chr1	1	2	t1
		chr2	2	3	t2
		chr3	3	4	t3
	}
	file copy -force tmp/d1/file1.tsv tmp/d2/same.tsv
	file_write tmp/d2/only2 ""
	cg tsvdiff tmp/d1 tmp/d2
} {Only in tmp/d1: tmp/d1/only1
Only in tmp/d2: tmp/d2/only2
diff tmp/d1/file1.tsv tmp/d2/file1.tsv
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

testsummarize
