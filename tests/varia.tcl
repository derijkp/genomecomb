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

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
