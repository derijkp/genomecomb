#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test bcol_index {basic} {
	test_cleantmp
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv tmp/temp.sft
	exec cg size tmp/temp.sft
} {47}

test bcol_index {basic} {
	test_cleantmp
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv tmp/temp.sft
	exec cg index tmp/temp.sft
	exec cg bcol get tmp/temp.sft.index/lines.bcol 0 3
} {76 103 228} 

test bcol_index {basic sib index not writable} {
	test_cleantmp
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv tmp/temp.sft
	file mkdir tmp/temp.sft.index
	file_write tmp/temp.sft.index/lines.bcol error
	file attributes tmp/temp.sft.index -permissions ugo-xw
	exec cg index tmp/temp.sft
	set bcolfile [indexdir_file tmp/temp.sft lines.bcol]
	exec cg bcol get $bcolfile 0 3
} {76 103 228} 

test bcol_make {basic} {
	test_cleantmp
	exec cg bcol make tmp/temp.bcol coverage < data/cov.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
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
	exec cg bcol make -t w tmp/temp.bcol coverage < data/cov.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
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
	exec cg bcol make -p pos -t c tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 10 51 > tmp/temp.test
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
< 39	0}

test bcol_make {-p} {
	test_cleantmp
	cg select -q {$chromosome == "chr1"} data/cov.tsv tmp/temp.tsv
	cg bcol make -p pos tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 50 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage}

test bcol_make {-p with gaps} {
	test_cleantmp
	cg select -q {$chromosome == "chr1" && $pos != 10} data/cov.tsv tmp/temp.tsv
	cg bcol make -p pos tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 50 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage
12d11
< 10	0}

test bcol_make {-p with gaps and -d} {
	test_cleantmp
	cg select -q {$chromosome == "chr1" && $pos != 10} data/cov.tsv tmp/temp.tsv
	cg bcol make -p pos -d 10 tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 50 > tmp/temp.test
	exec cg select -f {pos	coverage} tmp/temp.tsv tmp/temp.test2
	catch {exec diff tmp/temp.test tmp/temp.test2} e
	regsub {child process exited abnormally} $e {} e
	string trim $e
} {1c1
< pos	value
---
> pos	coverage
12d11
< 10	10}

test bcol_make {-c} {
	test_cleantmp
	exec cg bcol make -c chromosome tmp/temp.bcol coverage < data/cov.tsv
	exec cg bcol table -s 0 -c 2 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
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
	catch {exec cg bcol make -p pos -c chromosome tmp/temp.bcol coverage < data/cov.tsv}
	exec cg bcol table -s 0 -c 2 tmp/temp.bcol 0 > tmp/temp.test
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
< 32	0}

test bcol_make {-p and -c allchr} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp.bcol coverage < data/allchr.tsv
	exec cg bcol table -c all tmp/temp.bcol > tmp/temp.tsv
	cg select -f {chromosome=chr_clip($chromosome) pos value=$coverage} data/allchr.tsv tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test bcol_regextract {basic} {
	test_cleantmp
	cg select -q {$chromosome eq "chr2"} data/cov.tsv tmp/tempcov.tsv
	exec cg bcol make -p pos -c chromosome tmp/temp.bcol coverage < tmp/tempcov.tsv
	# exec cg bcol table -s 0 tmp/temp-chr2.bcol > tmp/temp.test
	cg regextract -above 1 10 tmp/temp.bcol
} {chromosome	begin	end
2	10	17
2	33	51}

test bcol_regextract {multiple chromosomes} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp.bcol coverage < data/cov.tsv
	cg regextract -above 1 10 tmp/temp.bcol
} {chromosome	begin	end
1	22	49
2	10	17
2	33	51}

test bcol_regextract {su} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome -t su tmp/temp.bcol coverage < data/cov.tsv
	# cg bcol table tmp/temp.bcol
	cg regextract -above 1 10 tmp/temp.bcol
} {chromosome	begin	end
1	22	49
2	10	17
2	33	51}

test bcol_regextract {s error} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome -t s tmp/temp- coverage < data/cov.tsv
} {conversion error for type s: value 60000 too large.
} error

test bcol_histo {basic} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp.bcol coverage < data/cov.tsv
	file_write tmp/reg.tsv "chromosome\tbegin\tend\tname\nchr1\t4\t20\tt1\nchr1\t40\t60\tt2\n"
	cg bcol_histo tmp/reg.tsv tmp/temp.bcol {1 8 10}
} {name	r<1	r1<8	r8<10	r10<	size	avg	min	max
t1	0	11	5	0	16	5.81	2	9
t2	10	0	0	10	20	6.25	0	14
----------
Total	10	11	5	10	36	6.06	0	14
Totalpercent	27.78	30.56	13.89	27.78}

test bcol_histo {2 chromosomes} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp.bcol coverage < data/cov.tsv
	file_write tmp/reg.tsv "chromosome\tbegin\tend\tname\nchr1\t4\t20\tt1\nchr1\t40\t60\tt2\nchr2\t30\t40\tt3\n"
	cg bcol_histo tmp/reg.tsv tmp/temp.bcol {1 8 10}
} {name	r<1	r1<8	r8<10	r10<	size	avg	min	max
t1	0	11	5	0	16	5.81	2	9
t2	10	0	0	10	20	6.25	0	14
t3	3	0	0	7	10	7749.80	0	60000
----------
Total	13	11	5	17	46	1689.48	0	60000
Totalpercent	28.26	23.91	10.87	36.96}

test bcol_histo {old multifile format} {
	test_cleantmp
	file copy {*}[glob data/old*bcol*] tmp
	file_write tmp/reg.tsv "chromosome\tbegin\tend\tname\nchr1\t4\t20\tt1\nchr1\t40\t60\tt2\n"
	cg bcol_histo tmp/reg.tsv tmp/old {1 8 10}
} {name	r<1	r1<8	r8<10	r10<	size	avg	min	max
t1	0	11	5	0	16	5.81	2	9
t2	10	0	0	10	20	6.25	0	14
----------
Total	10	11	5	10	36	6.06	0	14
Totalpercent	27.78	30.56	13.89	27.78}

test bcol_histo {old multifile format 2 chromosomes} {
	test_cleantmp
	file copy {*}[glob data/old*bcol*] tmp
	file_write tmp/reg.tsv "chromosome\tbegin\tend\tname\nchr1\t4\t20\tt1\nchr1\t40\t60\tt2\nchr2\t30\t40\tt3\n"
	cg bcol_histo tmp/reg.tsv tmp/old {1 8 10}
} {name	r<1	r1<8	r8<10	r10<	size	avg	min	max
t1	0	11	5	0	16	5.81	2	9
t2	10	0	0	10	20	6.25	0	14
t3	3	0	0	7	10	7749.80	0	60000
----------
Total	13	11	5	17	46	1689.48	0	60000
Totalpercent	28.26	23.91	10.87	36.96}

test bcol_make {c too large error} {
	test_cleantmp
	exec cg bcol make -t c tmp/temp.bcol coverage < data/cov.tsv
} {conversion error for type c: value 129 too large.
} error

test bcol_make {c too large error limits.tsv} {
	test_cleantmp
	exec cg bcol make -t c tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type c: value 128 too large.
} error

test bcol_make {c too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t c tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type c: value -128 too small.
} error

test bcol_make {cu too large error limits.tsv} {
	test_cleantmp
	exec cg bcol make -t cu tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type cu: value 256 too large.
} error

test bcol_make {cu too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t cu tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type cu: value -1 too small.
} error

test bcol_make {s too large error} {
	test_cleantmp
	exec cg bcol make -t s tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type s: value 32768 too large.
} error

test bcol_make {s too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t s tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type s: value -32769 too small.
} error

test bcol_make {su too large error} {
	test_cleantmp
	exec cg bcol make -t su tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type su: value 65536 too large.
} error

test bcol_make {su too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t su tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type su: value -1 too small.
} error

test bcol_make {i too large error} {
	test_cleantmp
	exec cg bcol make -t i tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type i, value 2147483648: Numerical result out of range.
} error

test bcol_make {i too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t i tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type i, value -2147483649: Numerical result out of range.
} error

test bcol_make {iu too large error} {
	test_cleantmp
	exec cg bcol make -t iu tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type iu, value 4294967296: Numerical result out of range.
} error

test bcol_make {iu too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t iu tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type iu, value -1: Numerical result out of range.
} error

test bcol_make {w too large error} {
	test_cleantmp
	exec cg bcol make -t w tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type w, value 9223372036854775808: Numerical result out of range.
} error

test bcol_make {w too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t w tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type w, value -9223372036854775809: Numerical result out of range.
} error

test bcol_make {wu too large error} {
	test_cleantmp
	exec cg bcol make -t wu tmp/temp.bcol coverage < data/limits.tsv
} {conversion error for type wu, value 18446744073709551616: Numerical result out of range.
} error

test bcol_make {wu too small error} {
	test_cleantmp
	cg select -q {$coverage < 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t wu tmp/temp.bcol coverage < tmp/temp.tsv
} {conversion error for type wu, value -1: Numerical result out of range.
} error

test bcol_make {types c} {
	test_cleantmp
	cg select -q {$coverage <= 127 && $coverage >= -127} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t c tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types cu} {
	test_cleantmp
	cg select -q {$coverage <= 255 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t cu tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types s} {
	test_cleantmp
	cg select -q {$coverage <= 32767 && $coverage >= -32768} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t s tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types su} {
	test_cleantmp
	cg select -q {$coverage <= 65535 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t su tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types i} {
	test_cleantmp
	cg select -q {$coverage <= 2147483647 && $coverage >= -2147483648} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t i tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types iu} {
	test_cleantmp
	cg select -q {$coverage <= 4294967295 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t iu tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types w} {
	test_cleantmp
	cg select -q {$coverage <= 9223372036854775807 && $coverage >= -9223372036854775808} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t w tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {types wu} {
	test_cleantmp
	cg select -q {$coverage <= 18446744073709551615 && $coverage >= 0} data/limits.tsv tmp/temp.tsv
	exec cg bcol make -t wu tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_make {-compress 1} {
	test_cleantmp
	cg select -q {$coverage <= 9223372036854775807 && $coverage >= -9223372036854775808} data/limits.tsv tmp/temp.tsv
	exec cg bcol make --compress 1 -t w tmp/temp.bcol coverage < tmp/temp.tsv
	exec cg bcol table -s 0 tmp/temp.bcol 0 | cg select -f value > tmp/temp.test
	exec cg select -f {value=$coverage} tmp/temp.tsv tmp/temp.test2
	exec diff tmp/temp.test tmp/temp.test2
} {}

test bcol_update {update and join old bcols} {
	test_cleantmp
	exec cg bcol update tmp/new.bcol data/old-chr1.bcol data/old-chr2.bcol
	exec cg bcol table -s 1 -c all tmp/new.bcol > tmp/cov.tsv.new
	exec cg select -f {{chromosome=chr_clip($chromosome)} pos value=$coverage} data/cov.tsv tmp/cov.tsv.ori
	exec diff tmp/cov.tsv.new tmp/cov.tsv.ori
} {59,74d58
< 2	17	0
< 2	18	0
< 2	19	0
< 2	20	0
< 2	21	0
< 2	22	0
< 2	23	0
< 2	24	0
< 2	25	0
< 2	26	0
< 2	27	0
< 2	28	0
< 2	29	0
< 2	30	0
< 2	31	0
< 2	32	0
child process exited abnormally} error

test bcol_make_multi {var-annot type d} {
	test_cleantmp
	cg bcol make -t d --multicol alt --multilist A,C,T,G -p begin -c chromosome tmp/temp.bcol score < data/var-annot.tsv
	cg bcol table -c all tmp/temp.bcol > tmp/result.tsv
	cg splitalleles tmp/result.tsv tmp/split.tsv
	cg select -f {chromosome begin=$pos end=$pos+1 alt {score=regsub($value,"\\.0","")}} -q {$value ne "0.0"} tmp/split.tsv tmp/temp.tsv
	cg collapsealleles tmp/temp.tsv > tmp/test.tsv
	exec diff tmp/test.tsv data/var-annot.tsv
} {}

test bcol_make_multi {var-annot type f} {
	test_cleantmp
	cg bcol make -t f --multicol alt --multilist A,C,T,G --compress 0 -p begin -c chromosome tmp/temp.bcol score < data/var-annot.tsv
	cg bcol table -p 1 -c all tmp/temp.bcol > tmp/result.tsv
	cg splitalleles tmp/result.tsv tmp/split.tsv
	cg select -f {chromosome begin=$pos end=$pos+1 alt {score=regsub($value,"\\.0","")}} -q {$value ne "0.0"} tmp/split.tsv tmp/temp.tsv
	cg collapsealleles tmp/temp.tsv > tmp/test.tsv
	exec diff tmp/test.tsv data/var-annot.tsv
} {2c2
< 1	0	1	C,G	0.1,1.1
---
> 1	0	1	A,C,G	1e-7,0.1,1.1
13c13
< 2	29	30	A,C,T	30,3000000,30.2
---
> 2	29	30	A,C,T	30,3000000.1,30.21
*} error match

test bcol_make {--precision} {
	test_cleantmp
	write_tab tmp/varannot.tsv {
		chromosome      begin   end     score
		1       0       1   0.401
	}
	cg bcol make --precision 2 -t f -p begin -c chromosome tmp/temp.bcol score < tmp/varannot.tsv
	cg bcol table -c all tmp/temp.bcol
} {chromosome	pos	value
1	0	0.40}

test bcol_make_multi {--precision} {
	test_cleantmp
	write_tab tmp/varannot.tsv {
		chromosome	begin	end	alt	score
		1	0	1	A,C,G	0.401,0.1,0.999999
	}
	cg bcol make --precision 2 -t f --multicol alt --multilist A,C,T,G --compress 0 -p begin -c chromosome tmp/temp.bcol score < tmp/varannot.tsv
	cg bcol table -c all tmp/temp.bcol
} {chromosome	pos	alt	value
1	0	A,C,T,G	0.40,0.10,0.00,1.00}

test bcol_make {-e (endpos) option} {
	test_cleantmp
	write_tab tmp/reg.tsv {
		chromosome	begin	end	score
		1	5	8	1
		1	8	10	2
		2	1	3	3
	}
	exec cg bcol make -t c -c chromosome -p begin -e end tmp/temp.bcol score < tmp/reg.tsv
	cg bcol table tmp/temp.bcol
} {chromosome	pos	value
1	5	1
1	6	1
1	7	1
1	8	2
1	9	2
2	1	3
2	2	3}

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
