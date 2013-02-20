#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test_cleantmp

test bcol_index {basic} {
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv tmp/temp.sft
	exec cg size tmp/temp.sft
} {46} 

test bcol_index {basic} {
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv tmp/temp.sft
	exec cg index tmp/temp.sft
	exec cg bcol get tmp/temp.sft.index/lines.bcol 0 3
} {76 103 228} 

test bcol_make {basic} {
	test_cleantmp
	exec cg bcol make tmp/temp coverage < data/cov.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f coverage data/cov.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< value
---
> coverage}

test bcol_make {wide} {
	test_cleantmp
	exec cg bcol make -t w tmp/temp coverage < data/cov.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f coverage data/cov.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< value
---
> coverage}

test bcol_make {c met -p} {
	test_cleantmp
	cg select -q {$chromosome == "chr2" && $coverage < 129} data/cov.tsv tmp/temp.tsv
	exec cg bcol make -p pos -t c tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 10 51 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage
9,24d8
< 17	0
< 18	0
< 19	0
< 20	0
< 21	0
< 22	0
< 23	0
< 24	0
< 25	0
< 26	0
< 27	0
< 28	0
< 29	0
< 30	0
< 31	0
< 32	0
28,31d11
< 36	0
< 37	0
< 38	0
< 39	0
43d22
<}

test bcol_make {-p} {
	test_cleantmp
	cg select -q {$chromosome == "chr1"} data/cov.tsv tmp/temp.tsv
	cg bcol make -p pos tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 50 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage
52d51
<}

test bcol_make {-p with gaps} {
	test_cleantmp
	cg select -q {$chromosome == "chr1" && $pos != 10} data/cov.tsv tmp/temp.tsv
	cg bcol make -p pos tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 50 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage
12d11
< 10	0
52d50
<}

test bcol_make {-p with gaps and -d} {
	test_cleantmp
	cg select -q {$chromosome == "chr1" && $pos != 10} data/cov.tsv tmp/temp.tsv
	cg bcol make -p pos -d 10 tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 50 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage
12d11
< 10	10
52d50
<}

test bcol_make {-c} {
	test_cleantmp
	catch {exec cg bcol make -c chromosome tmp/temp- coverage < data/cov.tsv}
	exec cg bcol table tmp/temp-chr2.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f coverage -q {$chromosome == "chr2"} data/cov.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< value
---
> coverage}

test bcol_make {-p and -c} {
	test_cleantmp
	catch {exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv}
	exec cg bcol table tmp/temp-chr2.bcol > tmp/temp.test
	exec cg select -f {pos	coverage} -q {$chromosome == "chr2"} data/cov.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1,11c1
< pos	value
< 0	0
< 1	0
< 2	0
< 3	0
< 4	0
< 5	0
< 6	0
< 7	0
< 8	0
< 9	0
---
> pos	coverage
19,34d8
< 17	0
< 18	0
< 19	0
< 20	0
< 21	0
< 22	0
< 23	0
< 24	0
< 25	0
< 26	0
< 27	0
< 28	0
< 29	0
< 30	0
< 31	0
< 32	0
53d26
<}

test bcol_make {-p and -c allchr: check bug} {
	test_cleantmp
	catch {exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/allchr.tsv} e
	set files [glob tmp/temp-*.bcol]
	foreach file $files {
		if {[file size $file] == 0} {error "file $file has size 0"}
	}
	set result {}
} {}

test bcol_regextract {basic} {
	test_cleantmp
	catch {exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv} e
	# exec cg bcol table tmp/temp-chr2.bcol > tmp/temp.test
	cg regextract -above 1 10 {*}[lsort -dict [glob tmp/temp-*.bcol]] 2> /dev/null
} {chromosome	begin	end
chr1	22	49
chr2	10	17
chr2	33	51}

test bcol_regextract {su} {
	test_cleantmp
	catch {exec cg bcol make -p pos -c chromosome -t su tmp/temp- coverage < data/cov.tsv} e
	# exec cg bcol table tmp/temp-chr2.bcol > tmp/temp.test
	cg regextract -above 1 10 {*}[lsort -dict [glob tmp/temp-*.bcol]] 2> /dev/null
} {chromosome	begin	end
chr1	22	49
chr2	10	17
chr2	33	51}

test bcol_regextract {s error} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome -t s tmp/temp- coverage < data/cov.tsv
} {conversion error for type s: value 60000 too large.
} error

test bcol_histo {basic} {
	test_cleantmp
	catch {exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv}
	file_write tmp/reg.tsv "chromosome\tbegin\tend\tname\nchr1\t4\t20\tt1\nchr1\t40\t60\tt2\n"
	cg bcol_histo tmp/reg.tsv tmp/temp-chr1.bcol {1 8 10}
} {name	r<1	r1<8	r8<10	r10<
t1	0	11	5	0
t2	10	0	0	10
----------
Total	10	11	5	10
Totalpercent	27.78	30.56	13.89	27.78}

test bcol_make {c too large error} {
	test_cleantmp
	exec cg bcol make -t c tmp/temp coverage < data/cov.tsv
} {conversion error for type c: value 129 too large.
} error

test bcol_make {c too large error limits.tsv} {
	test_cleantmp
	exec cg bcol make -t c tmp/temp coverage < data/limits.tsv
} {conversion error for type c: value 128 too large.
} error

test bcol_make {c too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t c tmp/temp coverage < tmp/temp.tsv
} {conversion error for type c: value -128 too small.
} error

test bcol_make {cu too large error limits.tsv} {
	test_cleantmp
	exec cg bcol make -t cu tmp/temp coverage < data/limits.tsv
} {conversion error for type cu: value 256 too large.
} error

test bcol_make {cu too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t cu tmp/temp coverage < tmp/temp.tsv
} {conversion error for type cu: value -1 too small.
} error

test bcol_make {s too large error} {
	test_cleantmp
	exec cg bcol make -t s tmp/temp coverage < data/limits.tsv
} {conversion error for type s: value 32768 too large.
} error

test bcol_make {s too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t s tmp/temp coverage < tmp/temp.tsv
} {conversion error for type s: value -32769 too small.
} error

test bcol_make {su too large error} {
	test_cleantmp
	exec cg bcol make -t su tmp/temp coverage < data/limits.tsv
} {conversion error for type su: value 65536 too large.
} error

test bcol_make {su too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t su tmp/temp coverage < tmp/temp.tsv
} {conversion error for type su: value -1 too small.
} error

test bcol_make {i too large error} {
	test_cleantmp
	exec cg bcol make -t i tmp/temp coverage < data/limits.tsv
} {conversion error for type i, value 2147483648: Numerical result out of range.
} error

test bcol_make {i too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t i tmp/temp coverage < tmp/temp.tsv
} {conversion error for type i, value -2147483649: Numerical result out of range.
} error

test bcol_make {iu too large error} {
	test_cleantmp
	exec cg bcol make -t iu tmp/temp coverage < data/limits.tsv
} {conversion error for type iu, value 4294967296: Numerical result out of range.
} error

test bcol_make {iu too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t iu tmp/temp coverage < tmp/temp.tsv
} {conversion error for type iu, value -1: Numerical result out of range.
} error

test bcol_make {w too large error} {
	test_cleantmp
	exec cg bcol make -t w tmp/temp coverage < data/limits.tsv
} {conversion error for type w, value 9223372036854775808: Numerical result out of range.
} error

test bcol_make {w too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t w tmp/temp coverage < tmp/temp.tsv
} {conversion error for type w, value -9223372036854775809: Numerical result out of range.
} error

test bcol_make {wu too large error} {
	test_cleantmp
	exec cg bcol make -t wu tmp/temp coverage < data/limits.tsv
} {conversion error for type wu, value 18446744073709551616: Numerical result out of range.
} error

test bcol_make {wu too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t wu tmp/temp coverage < tmp/temp.tsv
} {conversion error for type wu, value -1: Numerical result out of range.
} error

test bcol_make {types c} {
	test_cleantmp
	cg select -q {$coverage <= 127 && $coverage >= -127} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t c tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types cu} {
	test_cleantmp
	cg select -q {$coverage <= 255 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t cu tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types s} {
	test_cleantmp
	cg select -q {$coverage <= 32767 && $coverage >= -32768} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t s tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types su} {
	test_cleantmp
	cg select -q {$coverage <= 65535 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t su tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types i} {
	test_cleantmp
	cg select -q {$coverage <= 2147483647 && $coverage >= -2147483648} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t i tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types iu} {
	test_cleantmp
	cg select -q {$coverage <= 4294967295 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t iu tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types w} {
	test_cleantmp
	cg select -q {$coverage <= 9223372036854775807 && $coverage >= -9223372036854775808} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t w tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types wu} {
	test_cleantmp
	cg select -q {$coverage <= 18446744073709551615 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t wu tmp/temp coverage < tmp/temp.tsv
	exec cg bcol table tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
