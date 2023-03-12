#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

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

test multireg {basic} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg2.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg2.tsv
} {}

test multireg {basic compressed} {
	file delete tmp/temp.tsv
	file copy data/reg1.tsv data/reg2.tsv tmp
	cg zst {*}[glob tmp/*.tsv]
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv tmp/reg1.tsv.zst tmp/reg2.tsv.zst
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg2.tsv
} {}

test multireg {same} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg1b.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg1b.tsv
} {}

test multireg {(try to add) add existing} {
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg2.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg2.tsv
} {}

test multireg {(try to add) add existing compressed} {
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg2.tsv
	compress data/reg1.tsv tmp/reg1.tsv.zst
	exec cg multireg tmp/temp.tsv tmp/reg1.tsv.zst
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg2.tsv
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
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg1b-reg2.tsv
} {}

test multireg {different chromosome naming} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec diff tmp/temp.tsv data/expected-multireg-reg1-reg4.tsv
} {}

test multireg {sort error 1 in compar_file file} {
	file delete tmp/temp.tsv
	file copy data/vars_sorterror1.tsv tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
} {*File (*tmp/temp.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 1 in added file from new} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror1.tsv
} {*File (*data/vars_sorterror1.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 1 in added file} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror1.tsv
} {*File (*data/vars_sorterror1.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 2 in database file} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror2.tsv
} {*File (*data/vars_sorterror2.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test multireg {sort error 3 in database file} {
	file delete tmp/temp.tsv
	exec cg multireg tmp/temp.tsv data/reg1.tsv data/reg4.tsv
	exec cg multireg tmp/temp.tsv data/vars_sorterror3.tsv
} {*File (*data/vars_sorterror3.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

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

test multireg {-limireg} {
	set files [multiregtest 5]
	file delete tmp/temp.tsv
	file_write tmp/regfile.tsv [deindent {
		chromosome	begin	end
		1	30	50
		1	70	100
	}]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	sreg1	sreg2	sreg3	sreg4	sreg5
		1	30	35	1	1	1	1	1
		1	35	40	1	1	1	1	1
		1	40	45	1	1	1	1	1
		1	45	50	1	1	1	1	1
		1	70	75	0	0	0	1	1
		1	75	80	0	0	0	0	1
		1	80	100	0	0	0	0	0
	}]\n
	exec cg multireg -limitreg tmp/regfile.tsv tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multireg {5 files 8 distribute 1 reg one not reg (maxopenfiles) -limireg} {
	set files [multiregtest 5]
	# 8 maximum open files means 4 files can be processed at the same time
	# (-4 because stdout, etc. also count as "open files")
	file_write tmp/regfile.tsv [deindent {
		chromosome	begin	end
		1	30	50
		1	70	100
	}]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	sreg1	sreg2	sreg3	sreg4	sreg5
		1	30	35	1	1	1	1	1
		1	35	40	1	1	1	1	1
		1	40	45	1	1	1	1	1
		1	45	50	1	1	1	1	1
		1	70	75	0	0	0	1	1
		1	75	80	0	0	0	0	1
		1	80	100	0	0	0	0	0
	}]\n
	file delete tmp/temp.tsv
	exec cg multireg -m 8 -limitreg tmp/regfile.tsv tmp/temp.tsv {*}$files
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

test multireg {. in name} {
	file delete tmp/temp.tsv
	file copy data/regempty.tsv tmp/reg.empty.tsv
	exec cg multireg tmp/temp.tsv data/reg1b.tsv tmp/reg.empty.tsv
	file_read tmp/temp.tsv
} {chromosome	begin	end	reg1b	reg.empty
1	10	20	1	0
1	50	60	1	0
}

test multireg {. in name compressed} {
	file delete tmp/temp.tsv
	file copy data/regempty.tsv tmp/reg.empty.tsv
	cg zst tmp/reg.empty.tsv
	exec cg multireg tmp/temp.tsv data/reg1b.tsv tmp/reg.empty.tsv.zst
	file_read tmp/temp.tsv
} {chromosome	begin	end	reg1b	reg.empty
1	10	20	1	0
1	50	60	1	0
}

test multireg {give error if there are overlapping ranges in file} {
	file_write tmp/sreg1.tsv [deindent {
		chromosome	begin	end
		chr9	90017364	90017384
		chr9	90017467	90017468
		chr9	90017470	90017472
		chr9	90017474	90017476
		chr9	90017478	90017480
	}]
	file_write tmp/sreg2.tsv [deindent {
		chromosome	begin	end
		chr9	90017435	90017486
		chr9	90017467	90017468
		chr9	90017471	90017472
		chr9	90017475	90017476
		chr9	90017478	90017481
	}]
	cg multireg tmp/result.tsv tmp/sreg1.tsv tmp/sreg2.tsv
} {*File (*) contains overlapping region(s) (correct this first using "cg regjoin")
region chr9:90017435-90017486 overlaps chr9:90017467-90017468*} error match

test multireg {bugfix: check if overlapping ranges in file} {
	file_write tmp/sreg1.tsv [deindent {
		chromosome	begin	end
		chr1	10039	10042
		chr1	10043	10050
	}]
	file_write tmp/sreg2.tsv [deindent {
		chromosome	begin	end
		chr1	10023	10024
		chr1	10025	10026
	}]
	cg multireg tmp/result.tsv tmp/sreg1.tsv tmp/sreg2.tsv
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

test regsubtract {compressed} {
	cg zst -o tmp/reg1.tsv.zst data/reg1.tsv
	cg zst -o tmp/reg2.tsv.zst data/reg2.tsv
	exec cg regsubtract tmp/reg1.tsv.zst tmp/reg2.tsv.zst
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

test regsubtract {sort error chromosome} {
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		10	10	25
		1	45	60
	}
	write_tab tmp/sub.tsv {
		chromosome	begin	end
		1	50	80
	}
	exec cg regsubtract tmp/sorterror.tsv tmp/sub.tsv > tmp/result.tsv
} {File (tmp/sorterror.tsv) is not correctly sorted (sort correctly using "cg select -s -")
10:10-10 came before 1:45-45
child process exited abnormally} error

test regsubtract {sort error chromosome, file2 overlap} {
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		10	10	25
		1	45	60
	}
	write_tab tmp/sub.tsv {
		chromosome	begin	end
		10	15	20
	}
	exec cg regsubtract tmp/sorterror.tsv tmp/sub.tsv > tmp/result.tsv
} {File (tmp/sorterror.tsv) is not correctly sorted (sort correctly using "cg select -s -")
10:10-10 came before 1:45-45
child process exited abnormally} error

test regsubtract {sort error chromosome, file2 later} {
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		10	10	25
		1	45	60
	}
	write_tab tmp/sub.tsv {
		chromosome	begin	end
		10	50	80
	}
	exec cg regsubtract tmp/sorterror.tsv tmp/sub.tsv > tmp/result.tsv
} {File (tmp/sorterror.tsv) is not correctly sorted (sort correctly using "cg select -s -")
10:10-10 came before 1:45-45
child process exited abnormally} error

test regsubtract {sort error chromosome file2} {
	write_tab tmp/f1.tsv {
		chromosome	begin	end
		1	5	8
		20	10	25
	}
	write_tab tmp/sorterror.tsv {
		chromosome	begin	end
		1	5	8
		10	10	25
		1	45	60
	}
	exec cg regsubtract tmp/f1.tsv tmp/sorterror.tsv > tmp/result.tsv
} {File (tmp/sorterror.tsv) is not correctly sorted (sort correctly using "cg select -s -")
10:10-10 came before 1:45-45
child process exited abnormally} error

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
	cg select -f {name=$chromosome test begin end} data/reg1.tsv tmp/temp.tsv
	exec cg covered -n name tmp/temp.tsv
} {chromosome	bases
1	20
2	130
3	200
M	10
X	100
Y	1000

total	1460}

test covered {name error} {
	cg select -f {name=$chromosome test begin end} data/reg1.tsv tmp/temp.tsv
	exec cg covered tmp/temp.tsv
} {header error: some fields (or alternatives) not found} error

test getregions {above} {
	exec getregions < data/coverage-chr1.tsv chr1 -1 0 1 11 {} 0 1
} {chr1	25	27
chr1	28	31
chr1	40	42}

test getregions {below} {
	exec getregions < data/coverage-chr1.tsv chr1 -1 0 1 0 9 0 1
} {chr1	20	24
chr1	27	28
chr1	42	43}

test getregions {shift} {
	exec getregions < data/coverage-chr1.tsv chr1 -1 0 1 11 {} -1 1
} {chr1	24	26
chr1	27	30
chr1	39	41}

test cg_regextract {basic} {
	exec cg regextract -max 9 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	20	24
chr1	27	28
chr1	42	43}

test cg_regextract {-min} {
	exec cg regextract --verbose 0 -qfields coverage -min 10 -shift 0 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	24	27
chr1	28	31
chr1	40	42}

test cg_regextract {above with chromosome field} {
	exec cg regextract --verbose 0 -qfields coverage -shift 0 -min 10 data/coverage.tsv
} {chromosome	begin	end
chr1	24	27
chr1	28	31
chr1	40	42
chr2	24	27
chr2	28	30}

test cg_regextract {-max} {
	exec cg regextract --verbose 0 -qfields coverage -shift 0 -max 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	20	25
chr1	27	28
chr1	42	43}

test cg_regextract {-min and -max} {
	exec cg regextract --verbose 0 -qfields coverage -shift 0 -min 10 -max 12 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	24	26}

test cg_regextract {shift} {
	exec cg regextract --verbose 0 -qfields coverage -shift -1 -min 11 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	24	26
chr1	27	30
chr1	39	41}

test cg_regextract {bam} {
	file copy data/bwa.bam tmp/bwa.bam
	exec samtools index tmp/bwa.bam
	exec cg regextract --verbose 0 -qfields coverage -shift -1 -min 8 tmp/bwa.bam
} {chromosome	begin	end
chr21	42735714	42735716
chr21	42735718	42735719
chr21	42735723	42735741
chr21	42775213	42775220
chr21	42775232	42775235
chr21	42775291	42775313
chr21	42779907	42779919
chr21	42779986	42780044
chr22	41923406	41923413}

test cg_regextract {bam -region} {
	file copy data/bwa.bam tmp/bwa.bam
	exec samtools index tmp/bwa.bam
	exec cg regextract --verbose 0 -region chr22 -qfields coverage -shift -1 -min 8 tmp/bwa.bam
} {chromosome	begin	end
chr22	41923406	41923413}

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

test regjoin {compressed} {
	cg zst -o tmp/reg1.tsv.zst data/reg1.tsv
	cg zst -o tmp/reg2.tsv.zst data/reg2.tsv
	exec cg regjoin tmp/reg1.tsv.zst tmp/reg2.tsv.zst
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

test regjoin {stdin} {
	cg cat -m 1 data/reg1.tsv data/reg2.tsv | cg select -s - > tmp/reg.tsv
	exec cg regjoin < tmp/reg.tsv
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

test regjoin {-fields *} {
	write_tab tmp/test.tsv { 
		chromosome	begin	end f1 f2
		1	10	15	1	2
		1	15	20	1	2
		1	20	25	3	4
		1	25	30	3	4
		1	30	35	3	5
		1	40	50	3	4
	}
	cg regjoin -fields * tmp/test.tsv
} {chromosome	begin	end	f1	f2
1	10	20	1	2
1	20	30	3	4
1	30	35	3	5
1	40	50	3	4}

test regjoin {-fields f1} {
	write_tab tmp/test.tsv { 
		chromosome	begin	end f1 f2
		1	10	15	1	2
		1	15	20	1	2
		1	20	25	3	4
		1	25	30	3	4
		1	30	35	3	5
		1	40	50	3	4
	}
	cg regjoin -fields {f1} tmp/test.tsv
} {chromosome	begin	end	f1
1	10	20	1
1	20	35	3
1	40	50	3}

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

test regjoin {empty files} {
	write_tab tmp/f1.tsv {}
	exec cg regjoin tmp/f1.tsv > tmp/result.tsv
} {header error: fields (or alternatives) not found: chromosome begin end} error

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

test reg_subtract {basic} {
	write_tab tmp/reg1.tsv {
		chromosome	begin	end
		chr1	5	8
		chr1	20	30
	}
	write_tab tmp/reg2.tsv {
		chromosome	begin	end
		chr1	6	8
		chr1	10	20
		chr1	25	28
	}
	exec cg regsubtract tmp/reg1.tsv tmp/reg2.tsv
} {chromosome	begin	end
chr1	5	6
chr1	20	25
chr1	28	30}

test reg_subtract {diff chr notation} {
	write_tab tmp/reg1.tsv {
		chromosome	begin	end
		1	5	8
		1	20	30
		2	10	20
	}
	write_tab tmp/reg2.tsv {
		chromosome	begin	end
		chr1	6	8
		chr1	10	20
		chr1	25	28
		chr1	100	110
		chr2	16	30
	}
	exec cg regsubtract tmp/reg1.tsv tmp/reg2.tsv
} {chromosome	begin	end
1	5	6
1	20	25
1	28	30
2	10	16}

test regselect {basic} {
	exec cg regselect data/vars1.tsv data/reg_annot.tsv > tmp/temp.tsv
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {basic piped} {
	exec cat data/vars1.tsv | cg regselect -stack 1 - data/reg_annot.tsv > tmp/temp.tsv
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {basic piped without -} {
	exec cat data/vars1.tsv | cg regselect -stack 1 data/reg_annot.tsv > tmp/temp.tsv
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {regselect -o} {
	exec cg regselect -o tmp/temp.tsv data/vars1.tsv data/reg_annot.tsv
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {regselect -o from stdin} {
	exec cat data/vars1.tsv | cg regselect -stack 1 -o tmp/temp.tsv - data/reg_annot.tsv
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {regselect -o compressed} {
	exec cg regselect -o tmp/temp.tsv.zst data/vars1.tsv data/reg_annot.tsv
	exec cg select -rf {list} tmp/temp.tsv.zst tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {regselect -o compressed from stdin} {
	exec cat data/vars1.tsv | cg regselect -o tmp/temp.tsv.zst - data/reg_annot.tsv
	exec cg select -rf {list} tmp/temp.tsv.zst tmp/temp2.tsv
	exec cg select -q {$regtest != ""} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} data/expected-vars1-reg_annot.tsv tmp/tempexpected.tsv
	exec diff tmp/temp2.tsv tmp/tempexpected.tsv
} {}

test regselect {basic2} {
	exec cg regselect data/vars1.tsv data/reg_annot.tsv > tmp/temp.tsv
	exec cg select -f {chromosome begin end type} tmp/temp.tsv
} {chromosome	begin	end	type
chr1	4000	4001	snp
chr1	4001	4002	snp
chr1	4099	4100	snp
chr1	5000	5010	del
chr2	5000	5000	ins
chr2	5000	5001	snp}

test regselect {near} {
	exec cg regselect data/vars1.tsv data/reg_annot.tsv 10 > tmp/temp.tsv
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

test regions_skip {basic} {
	list \
	[regions_skip chrM {chrM M}] \
	[regions_skip M {chrM M}] \
	[regions_skip chrM:0-100 {chrM M}] \
	[regions_skip M:0-100 {chrM M}] \
	[regions_skip chr2 {chrM M}] \
	[regions_skip 2 {chrM M}] \
	[regions_skip chr2:0-100 {chrM M}] \
	[regions_skip 2:0-100 {chrM M}] \
} {1 1 1 1 0 0 0 0}

set ::env(PATH) $keeppath

file delete -force tmp/temp.tsv tmp/temp.tsv.old

testsummarize
