#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

test alignedseq {basic} {
	foreach {chromosome	begin end cigar	seq expected} {
		chr1	99	119	20M	ABCDEFGHIJKLMNOPQRST IJKLMNOPQRST
		chr1	99	119	2M1X17M	ABCDEFGHIJKLMNOPQRST IJKLMNOPQRST
		chr1	99	117	2S18M	ABCDEFGHIJKLMNOPQRST KLMNOPQRST
		chr1	99	119	2S2M2D16M	ABCDEFGHIJKLMNOPQRST IJKLMNOPQRST
		chr1	99	115	2S2M2I14M	ABCDEFGHIJKLMNOPQRST MNOPQRST
		chr1	101	115	4S2M2I12M	ABCDEFGHIJKLMNOPQRST MNOPQRST
		chr1	101	115	4S2M2I12M	* ?
	} {
		set reqchr chr1
		set reqbegin 107
		set reqend 131
		file_write tmp/ali.tsv "chromosome\tbegin\tend\tcigar\tseq\n$chromosome\t$begin\t$end\t$cigar\t$seq\n"
		set result [alignedseq $chromosome $begin $end $cigar $seq $reqchr $reqbegin $reqend]
		if {$result ne $expected} {
			error "result \"$result\" should have been \"$expected\" for: $chromosome	$begin $end $cigar $seq $reqbegin $reqend"
		}
		set result [lindex [split [cg select -stack 1 -f "ali=alignedseq(\"$reqchr:$reqbegin-$reqend\")" tmp/ali.tsv] \n] end]
		if {$result ne $expected} {
			putsvars reqchr reqbegin reqend expected
			error "$reqchr:$reqbegin-$reqend select error, result: \n$result\n instead of expected \n$expected"
		}
		set result [lindex [split [cg select -stack 1 -f "ali=alignedseq($reqchr,$reqbegin,$reqend)" tmp/ali.tsv] \n] end]
		if {$result ne $expected} {
			putsvars reqchr reqbegin reqend expected
			error "$reqchr:$reqbegin-$reqend select error, result: \n$result\n instead of expected \n$expected"
		}
	}
} {}

test alignedseq {20M} {
	# 99 is 100 in half-open -> range ali is 100-120
	# 100  105  110  115  120
	# ABCDEFGHIJKLMNOPQRST
	# ABCDEFGHIJKLMNOPQRST
	#   2  5    10   15   20
	foreach {chromosome	begin end cigar	seq} {
		chr1	100	120	20M	ABCDEFGHIJKLMNOPQRST
	} break
	file_write tmp/ali.tsv "chromosome\tbegin\tend\tcigar\tseq\n$chromosome\t$begin\t$end\t$cigar\t$seq\n"
	foreach {chr reqbegin reqend expected} {
		chr1 99 100 {}
		chr1 99 101 {A}
		chr1 99 110 {ABCDEFGHIJ}
		chr1 99 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 99 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 99 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 101 {A}
		chr1 100 110 {ABCDEFGHIJ}
		chr1 100 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 100 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 101 102 {B}
		chr1 101 110 {BCDEFGHIJ}
		chr1 101 119 {BCDEFGHIJKLMNOPQRS}
		chr1 101 120 {BCDEFGHIJKLMNOPQRST}
		chr1 101 121 {BCDEFGHIJKLMNOPQRST}
		chr2 99 100 {}
	} {
		set result [alignedseq $chromosome $begin $end $cigar $seq chr1 $reqbegin $reqend]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend had result: \n$result\n instead of expected \n$expected"
		}
		set result [lindex [split [cg select -stack 1 -f "ali=alignedseq(\"$chr:$reqbegin-$reqend\")" tmp/ali.tsv] \n] end]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend select error, result: \n$result\n instead of expected \n$expected"
		}
		set result [lindex [split [cg select -stack 1 -f "ali=alignedseq($chr,$reqbegin,$reqend)" tmp/ali.tsv] \n] end]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend select error, result: \n$result\n instead of expected \n$expected"
		}
	}
} {}

test alignedseq {2S20M} {
	# 99 is 100 in half-open -> range ali is 100-120
	#   100  105  110  115  120
	#   ABCDEFGHIJKLMNOPQRST
	# YZABCDEFGHIJKLMNOPQRST
	#   2  5    10   15   20
	foreach {chromosome	begin end cigar	seq} {
		chr1	100	120	2S20M	YZABCDEFGHIJKLMNOPQRST
	} break
	foreach {chr reqbegin reqend expected} {
		chr1 99 100 {}
		chr1 99 101 {A}
		chr1 99 110 {ABCDEFGHIJ}
		chr1 99 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 99 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 99 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 101 {A}
		chr1 100 110 {ABCDEFGHIJ}
		chr1 100 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 100 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 101 102 {B}
		chr1 101 110 {BCDEFGHIJ}
		chr1 101 119 {BCDEFGHIJKLMNOPQRS}
		chr1 101 120 {BCDEFGHIJKLMNOPQRST}
		chr1 101 121 {BCDEFGHIJKLMNOPQRST}
	} {
		set result [alignedseq $chromosome $begin $end $cigar $seq chr1 $reqbegin $reqend]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend had result: \n$result\n instead of expected \n$expected"
		}
	}
} {}

test alignedseq {2S20M2S} {
	# 99 is 100 in half-open -> range ali is 100-120
	#   100  105  110  115  120
	#   ABCDEFGHIJKLMNOPQRST
	# YZABCDEFGHIJKLMNOPQRSTXY
	#   2  5    10   15   20
	foreach {chromosome	begin end cigar	seq} {
		chr1	100	120	2S20M2S	YZABCDEFGHIJKLMNOPQRSTXY
	} break
	foreach {chr reqbegin reqend expected} {
		chr1 99 100 {}
		chr1 99 101 {A}
		chr1 99 110 {ABCDEFGHIJ}
		chr1 99 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 99 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 99 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 101 {A}
		chr1 100 110 {ABCDEFGHIJ}
		chr1 100 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 100 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 101 102 {B}
		chr1 101 110 {BCDEFGHIJ}
		chr1 101 119 {BCDEFGHIJKLMNOPQRS}
		chr1 101 120 {BCDEFGHIJKLMNOPQRST}
		chr1 101 121 {BCDEFGHIJKLMNOPQRST}
	} {
		set result [alignedseq $chromosome $begin $end $cigar $seq chr1 $reqbegin $reqend]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend had result: \n$result\n instead of expected \n$expected"
		}
	}
} {}

test alignedseq {2H20M2H} {
	# 99 is 100 in half-open -> range ali is 100-120
	#   100  105  110  115  120
	#   ABCDEFGHIJKLMNOPQRST
	# YZABCDEFGHIJKLMNOPQRSTXY
	#   2  5    10   15   20
	foreach {chromosome	begin end cigar	seq} {
		chr1	100	120	2H20M2H	ABCDEFGHIJKLMNOPQRST
	} break
	foreach {chr reqbegin reqend expected} {
		chr1 99 100 {}
		chr1 99 101 {A}
		chr1 99 110 {ABCDEFGHIJ}
		chr1 99 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 99 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 99 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 101 {A}
		chr1 100 110 {ABCDEFGHIJ}
		chr1 100 119 {ABCDEFGHIJKLMNOPQRS}
		chr1 100 120 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 121 {ABCDEFGHIJKLMNOPQRST}
		chr1 101 102 {B}
		chr1 101 110 {BCDEFGHIJ}
		chr1 101 119 {BCDEFGHIJKLMNOPQRS}
		chr1 101 120 {BCDEFGHIJKLMNOPQRST}
		chr1 101 121 {BCDEFGHIJKLMNOPQRST}
	} {
		set result [alignedseq $chromosome $begin $end $cigar $seq chr1 $reqbegin $reqend]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend had result: \n$result\n instead of expected \n$expected"
		}
	}
} {}

test alignedseq {2S2M2I16M2S} {
	# 99 is 100 in half-open -> range ali is 100-118 with insertion at 102
	#   100 102     110     118
	#   AB--EFGHIJKLMNOPQRST
	# YZABCDEFGHIJKLMNOPQRSTXY
	#   2 4 6   10   15   20
	foreach {chromosome	begin end cigar	seq} {
		chr1	100	118	2S2M2I16M2S	YZABCDEFGHIJKLMNOPQRSTXY
	} break
	foreach {chr reqbegin reqend expected} {
		chr1 99 100 {}
		chr1 99 101 {A}
		chr1 99 110 {ABCDEFGHIJKL}
		chr1 99 117 {ABCDEFGHIJKLMNOPQRS}
		chr1 99 118 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 101 {A}
		chr1 100 110 {ABCDEFGHIJKL}
		chr1 100 117 {ABCDEFGHIJKLMNOPQRS}
		chr1 100 118 {ABCDEFGHIJKLMNOPQRST}
		chr1 101 102 {B}
		chr1 101 110 {BCDEFGHIJKL}
		chr1 101 117 {BCDEFGHIJKLMNOPQRS}
		chr1 101 118 {BCDEFGHIJKLMNOPQRST}
		chr1 102 103 {E}
		chr1 102 110 {EFGHIJKL}
		chr1 102 117 {EFGHIJKLMNOPQRS}
		chr1 102 118 {EFGHIJKLMNOPQRST}
		chr1 117 118 {T}
		chr1 118 119 {}
		chr1 119 120 {}
	} {
		set result [alignedseq $chromosome $begin $end $cigar $seq chr1 $reqbegin $reqend]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend had result: \n$result\n instead of expected \n$expected"
		}
	}
	set result [alignedseq chr1 100 118 2S2M2I16M2S YZABCDEFGHIJKLMNOPQRSTXY chr1 101 102]
	if {$result ne "B"} {error "wrong result \"$result\", should have been \"B\""}
	set result [alignedseq chr1 100 118 2S2M2I16M2S YZABCDEFGHIJKLMNOPQRSTXY chr1 102 103]
	if {$result ne "E"} {error "wrong result \"$result\", should have been \"E\""}
	set result [alignedseq chr1 100 118 2S2M2I16M2S YZABCDEFGHIJKLMNOPQRSTXY chr1 101 102 1]
	if {$result ne "B"} {error "wrong result \"$result\", should have been \"B\""}
	set result [alignedseq chr1 100 118 2S2M2I16M2S YZABCDEFGHIJKLMNOPQRSTXY chr1 102 103 1]
	if {$result ne "E"} {error "wrong result \"$result\", should have been \"E\""}
} {}

test alignedseq {2S2M2D18M2S} {
	# 99 is 100 in half-open -> range ali is 100-120 with del at 102
	#   100 104   110  115  120
	#   ABxxCDEFGHIJKLMNOPQRST
	# YZAB--CDEFGHIJKLMNOPQRSTXY
	#   2   4     10   15   20
	foreach {chromosome	begin end cigar	seq} {
		chr1	100	122	2S2M2D18M2S	YZABCDEFGHIJKLMNOPQRSTXY
	} break
	foreach {chr reqbegin reqend expected} {
		chr1 99 100 {}
		chr1 99 101 {A}
		chr1 99 110 {ABCDEFGH}
		chr1 99 119 {ABCDEFGHIJKLMNOPQ}
		chr1 99 120 {ABCDEFGHIJKLMNOPQR}
		chr1 99 121 {ABCDEFGHIJKLMNOPQRS}
		chr1 99 122 {ABCDEFGHIJKLMNOPQRST}
		chr1 100 101 {A}
		chr1 100 110 {ABCDEFGH}
		chr1 100 119 {ABCDEFGHIJKLMNOPQ}
		chr1 100 120 {ABCDEFGHIJKLMNOPQR}
		chr1 100 121 {ABCDEFGHIJKLMNOPQRS}
		chr1 100 122 {ABCDEFGHIJKLMNOPQRST}
		chr1 101 102 {B}
		chr1 101 110 {BCDEFGH}
		chr1 101 119 {BCDEFGHIJKLMNOPQ}
		chr1 101 120 {BCDEFGHIJKLMNOPQR}
		chr1 101 121 {BCDEFGHIJKLMNOPQRS}
		chr1 101 122 {BCDEFGHIJKLMNOPQRST}
		chr1 104 105 {C}
		chr1 104 110 {CDEFGH}
		chr1 104 119 {CDEFGHIJKLMNOPQ}
		chr1 104 120 {CDEFGHIJKLMNOPQR}
		chr1 104 121 {CDEFGHIJKLMNOPQRS}
		chr1 104 122 {CDEFGHIJKLMNOPQRST}
		chr1 104 123 {CDEFGHIJKLMNOPQRST}
		chr1 120 121 {S}
		chr1 121 122 {T}
		chr1 102 103 {-}
		chr1 122 123 {}
		chr1 123 124 {}
	} {
		set result [alignedseq $chromosome $begin $end $cigar $seq chr1 $reqbegin $reqend]
		if {$result ne $expected} {
			putsvars chr reqbegin reqend expected
			error "$chr:$reqbegin-$reqend had result: \n$result\n instead of expected \n$expected"
		}
	}
} {}

test alignedseq {bugfix} {
	cg exec {package require genomecomb
	        alignedseq chr1 99 119 20M ABCDEFGHIJKLMNOPQRST chr1 107 131
	}
	set commands {package require genomecomb
	        alignedseq chr1 99 119 20M ABCDEFGHIJKLMNOPQRST chr1 107 131
	}
	eval $commands
} IJKLMNOPQRST

test alignedseq {next to insert} {
}

testsummarize
