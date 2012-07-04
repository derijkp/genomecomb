#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

catch {file delete -force {*}[glob tmp/*]}

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

test bcol_make {-p} {
	file delete tmp/temp.bcol
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

test bcol_make {-c} {
	file delete tmp/temp.bcol
	catch {exec cg bcol make -c chromosome tmp/temp coverage < data/cov.tsv}
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
	file delete tmp/temp.bcol
	catch {exec cg bcol make -p pos -c chromosome tmp/temp coverage < data/cov.tsv}
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

test bcol_regextract {basic} {
	catch {file delete {*}[glob tmp/temp*.bcol]}
	catch {exec cg bcol make -p pos -c chromosome tmp/temp coverage < data/cov.tsv} e
	# exec cg bcol table tmp/temp-chr2.bcol > tmp/temp.test
	cg regextract -above 1 10 {*}[lsort -dict [glob tmp/temp-*.bcol]] 2> /dev/null
} {chromosome	begin	end
chr1	22	49
chr2	10	17
chr2	33	51}

test bcol_regextract {su} {
	catch {file delete {*}[glob tmp/temp*.bcol]}
	catch {exec cg bcol make -p pos -c chromosome -t su tmp/temp coverage < data/cov.tsv} e
	# exec cg bcol table tmp/temp-chr2.bcol > tmp/temp.test
	cg regextract -above 1 10 {*}[lsort -dict [glob tmp/temp-*.bcol]] 2> /dev/null
} {chromosome	begin	end
chr1	22	49
chr2	10	17
chr2	33	51}

test bcol_regextract {s error} {
	catch {file delete {*}[glob tmp/temp*.bcol]}
	exec cg bcol make -p pos -c chromosome -t s tmp/temp coverage < data/cov.tsv
} {value 60000 too large for type s} error

test bcol_histo {basic} {
	file delete tmp/temp.bcol
	catch {exec cg bcol make -p pos -c chromosome tmp/temp coverage < data/cov.tsv}
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

file delete -force {*}[glob tmp/*]

set ::env(PATH) $keeppath

testsummarize
