#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test multireg {basic} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg2.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg2.sft
} {}

test multireg {same} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg1b.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg1b.sft
} {}

test multireg {add empty} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1b.tsv data/regempty.tsv
	file_read tmp/temp.tsv
} {chromosome	begin	end	reg1b	regempty
1	10	20	1	0
1	50	60	1	0
}

test multireg {add empty first} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/regempty.tsv data/reg1b.tsv
	file_read tmp/temp.tsv
} {chromosome	begin	end	regempty	reg1b
1	10	20	0	1
1	50	60	0	1
}

test multireg {add fully empty} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1b.tsv data/empty.tsv
} {header error: fields \(or alternatives\) not found: chromosome begin end} error regexp

test multireg {add header error} {
	file delete tmp/temp.tsv
	file_write tmp/error.tsv "chromosome\ntest"
	exec cg multireg tmp/temp.tsv data/reg1b.tsv tmp/error.tsv
} {header error: fields \(or alternatives\) not found: chromosome begin end} error regexp

test multireg {3 adds} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg1b.tsv data/reg2.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg1b-reg2.sft
} {}

test multireg {different chromosome naming} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg4.sft
} {}

test multireg {sort error 1 in compar_file file} {
	file delete tmp/temp.tsv
	file copy data/vars_sorterror1.sft tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
} {*File (*tmp/temp.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 1 in added file from new} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror1.sft
} {*File (*data/vars_sorterror1.sft) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 1 in added file} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror1.sft
} {*File (*data/vars_sorterror1.sft) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 2 in database file} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror2.sft
} {*File (*data/vars_sorterror2.sft) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 3 in database file} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror3.sft
} {*File (*data/vars_sorterror3.sft) is not correctly sorted (sort correctly using "cg select -s -")*} error match

proc multiregtest {num} {
	set start 10
	set files {}
	for {set i 1} {$i <= $num} {incr i} {
		lappend files tmp/sreg$i.tsv
		set f [open tmp/sreg$i.tsv w]
		puts $f chromosome\tbegin\tend
		set pos $start
		for {set j 1} {$j <= 5} {incr j} {
			puts $f 1\t$pos\t[expr {$pos+10}]
			incr pos 10
		}
		close $f
		incr start 5
	}
	set f [open tmp/expected.tsv w]
	set header {chromosome begin end}
	for {set i 1} {$i <= $num} {incr i} {
		lappend header sreg$i
	}
	puts $f [join $header \t]
	set max [expr {$num*5+50}]
	for {set start 10} {$start <= $max} {incr start 5} {
		set line [list 1 $start [expr {$start+5}]]
		for {set i 1} {$i <= $num} {incr i} {
			set istart [expr {10 + ($i-1)*5}]
			set iend [expr {$istart + 50}]
			if {$start >= $istart && $start < $iend} {
				lappend line 1
			} else {
				lappend line 0
			}
		}
		puts $f [join $line \t]
	}
	close $f
	return $files
}

test multireg {10 files} {
	set files [multiregtest 10]
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multireg {10 files distribute (maxopenfiles)} {
	set files [multiregtest 10]
	file delete tmp/temp.tsv
	exec cg multireg -m 5 tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multireg {5 files 8 distribute 1 reg one not reg (maxopenfiles)} {
	set files [multiregtest 5]
	file delete tmp/temp.tsv
	# 8 maximum open files means 4 files can be processed at the same time
	# (-4 because stdout, etc. also count as "open files")
	exec cg multireg -m 8 tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multireg {6 files 8 distribute all reg (maxopenfiles)} {
	set files [multiregtest 6]
	file delete tmp/temp.tsv
	# 8 maximum open files means 4 files can be processed at the same time
	# (-4 because stdout, etc. also count as "open files")
	exec cg multireg -m 8 tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multireg {8 files 8 distribute all reg limit (maxopenfiles)} {
	set files [multiregtest 6]
	file delete tmp/temp.tsv
	# 8 maximum open files means 4 files can be processed at the same time
	# (-4 because stdout, etc. also count as "open files")
	exec cg multireg -m 8 tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multireg {100 files multilevel distribute (maxopenfiles)} {
	set files [multiregtest 100]
	file delete tmp/temp.tsv
	exec cg multireg -m 5 tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test regsubtract {basic} {
	exec cg regsubtract data/reg1.tsv data/reg2.tsv
} {chromosome	begin	end
1	10	15
1	55	60
2	100	150
2	160	170
2	180	200
3	2000	2100
Y	1000	1010
Y	1900	2000}

test regsubtract {bugfix: hang on file2 longer} {
	exec cg regsubtract data/reg1b.tsv data/reg4.tsv
} {chromosome	begin	end
1	10	15
1	55	60}

test covered {basic} {
	exec cg covered data/reg1.tsv
} {chromosome	bases
1	20
2	130
3	200
M	10
X	100
Y	1000

total	1460}

test covered {basic 2} {
	exec cg covered data/reg2.tsv
} {chromosome	bases
1	20
2	220
3	100
M	15
X	110
Y	890

total	1355}

test covered {name} {
	cg select -f {name=$chromosome test begin end} data/reg1.tsv temp.tsv
	exec cg covered -n name temp.tsv
} {chromosome	bases
1	20
2	130
3	200
M	10
X	100
Y	1000

total	1460}

test covered {name error} {
	cg select -f {name=$chromosome test begin end} data/reg1.tsv temp.tsv
	exec cg covered temp.tsv
} {header error: some fields (or alternatives) not found} error

test getregions {above} {
	exec getregions < data/coverage-chr1.tsv chr1 -1 0 1 10 1 0 1
} {chr1	25	27
chr1	28	31
chr1	40	42}

test getregions {below} {
	exec getregions < data/coverage-chr1.tsv chr1 -1 0 1 10 0 0 1
} {chr1	20	24
chr1	27	28
chr1	42	43}

test getregions {shift} {
	exec getregions < data/coverage-chr1.tsv chr1 -1 0 1 10 1 -1 1
} {chr1	24	26
chr1	27	30
chr1	39	41}

test cg_regextract {basic} {
	exec cg regextract 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	20	24
chr1	27	28
chr1	42	43}

test cg_regextract {above} {
	exec cg regextract --verbose 0 -qfields coverage -above 1 -shift 0 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	25	27
chr1	28	31
chr1	40	42}

test cg_regextract {above with chromosome field} {
	exec cg regextract --verbose 0 -qfields coverage -above 1 -shift 0 10 data/coverage.tsv
} {chromosome	begin	end
chr1	25	27
chr1	28	31
chr1	40	42
chr2	24	27
chr2	28	30}

test cg_regextract {below} {
	exec cg regextract --verbose 0 -qfields coverage -above 0 -shift 0 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	20	24
chr1	27	28
chr1	42	43}

test cg_regextract {shift} {
	exec cg regextract --verbose 0 -qfields coverage -above 1 -shift -1 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	24	26
chr1	27	30
chr1	39	41}

test regjoin {basic} {
	exec cg regjoin data/reg1.tsv data/reg2.tsv
} {chromosome	begin	end
1	10	25
1	45	60
2	100	200
2	300	500
3	1000	1100
3	2000	2100
M	10	25
X	90	200
Y	1000	2000}

test regjoin {basic} {
	exec cg regjoin data/reg1.tsv data/reg3.tsv
} {chromosome	begin	end
1	5	8
1	10	25
1	45	60
1	70	80
2	90	250
2	300	500
3	100	250
3	1000	2100
3	2500	2600
M	10	20
X	100	200
Y	1000	2000}

test regjoin {self} {
	exec cg regjoin data/reg3.tsv
} {chromosome	begin	end
1	5	8
1	15	18
1	19	25
1	45	50
1	70	80
2	90	100
2	200	250
2	300	500
3	100	250
3	1100	2000
3	2500	2600}

test regjoin {sort error chromosome} {
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		10	10	25
		1	45	60
	}
	exec cg regjoin tmp/sorterror.tsv > tmp/result.tsv
} {file1 is not correctly sorted (sort correctly using "cg select -s -")
child process exited abnormally} error

test regjoin {sort error chromosome file2} {
	write_tab tmp/f1.tsv {
		chromosome	begin	end
		1	5	8
	}
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		10	10	25
		1	45	60
	}
	exec cg regjoin tmp/f1.tsv tmp/sorterror.tsv > tmp/result.tsv
} {file2 is not correctly sorted (sort correctly using "cg select -s -")
child process exited abnormally} error

test regjoin {sort error begin} {
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		1	45	60
		1	10	25
	}
	exec cg regjoin tmp/sorterror.tsv > tmp/result.tsv
} {file1 is not correctly sorted (sort correctly using "cg select -s -")
child process exited abnormally} error

test regjoin {sort error begin file2} {
	write_tab tmp/f1.tsv {
		chromosome	begin	end
		1	5	8
	}
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		1	45	60
		1	10	25
	}
	exec cg regjoin tmp/f1.tsv tmp/sorterror.tsv > tmp/result.tsv
} {file2 is not correctly sorted (sort correctly using "cg select -s -")
child process exited abnormally} error

test regcollapse {basic} {
	exec cg regcollapse data/reg1.tsv data/reg2.tsv > tmp/test.tsv
	write_tab tmp/expected.tsv {
		chromosome	test	begin	end	test2
		1	t	10	15	{}
		1	t	15	20	,t2
		1	t	20	25	t2
		1	t	45	50	t2
		1	t	50	55	,t2
		1	t	55	60	{}
		2	t	100	150	{}
		2	t	150	160	,t2
		2	t	160	170	{}
		2	t	170	180	,t2
		2	t	180	200	{}
		2	t	300	450	t2
		2	t	450	480	,t2
		2	t	480	500	t2
		3	t	1000	1100	,t2
		3	t	2000	2100	{}
		M	t	10	20	,t2
		M	t	20	25	t2
		X	t	90	100	t2
		X	t	100	200	,t2
		Y	t	1000	1010	{}
		Y	t	1010	1900	,t2
		Y	t	1900	2000	{}
	}
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test regselect {basic} {
	exec cg regselect data/vars1.sft data/reg_annot.sft > tmp/temp.tsv
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} ../tests/data/expected-vars1-reg_annot.sft tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {basic2} {
	exec cg regselect data/vars1.sft data/reg_annot.sft > tmp/temp.tsv
	exec cg select -f {chromosome begin end type} tmp/temp.tsv
} {chromosome	begin	end	type
chr1	4000	4001	snp
chr1	4001	4002	snp
chr1	4099	4100	snp
chr1	5000	5010	del
chr2	5000	5000	ins
chr2	5000	5001	snp}

test regselect {near} {
	exec cg regselect data/vars1.sft data/reg_annot.sft 10 > tmp/temp.tsv
	exec cg select -f {chromosome begin end type} tmp/temp.tsv
} {chromosome	begin	end	type
chr1	4000	4001	snp
chr1	4001	4002	snp
chr1	4099	4100	snp
chr1	5000	5010	del
chr1	5020	5021	snp
chr2	5000	5000	ins
chr2	5000	5001	snp
chr2	5005	5006	snp
chr2	5010	5011	snp
chr2	5011	5012	snp}

set ::env(PATH) $keeppath

file delete -force tmp/temp.tsv tmp/temp.tsv.old

testsummarize
