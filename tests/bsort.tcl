#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

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
} {chr1 chr2 chr10 chrX *}

test bsort {basic 2} {
	bsort {a1 b10 b2 a2 b1 1 *}
} {1 a1 a2 b1 b2 b10 *}

test bsort {-index 1} {
	bsort -index 1 {{1 a1} {2 b10} {3 b2} {4 a2} {5 b1} {6 1} {7 *}}
} {{6 1} {1 a1} {4 a2} {5 b1} {3 b2} {2 b10} {7 *}}

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

test bsort {nums and combinations} {
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
		# lappend list "a-$number"
		lappend list "${number}a"
		lappend list "${number}e"
		lappend list "${number}e+"
		lappend list "a $number"
		lappend list "a\t$number"
		lappend list "a\t$number\tb"
		lappend list "a $number b"
		lappend list "$number a"
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

	set list [list { } 0 1 2 3 4 5 6 7 8 9 ! \" \# $ % & \' ( ) , - . / + : \; < = > ? @ A a B b C c D d E e F f G g H h I i J j K k L l M m N n O o P p Q q R r S s T t U u V v W w X x Y y Z z \[ \\ \] ^ _ ` \{ | \} ~ *]
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	set olist [list { } ! \" \# $ % & \' ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : \; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z \[ \\ \] ^ _ ` a b c d e f g h i j k l m n o p q r s t u v w x y z \{ | \} ~]
	set test2 [bsort $olist]
	if {$test2 ne $list} {error "wrong sort"}
	checkbsort $list

} {}

test bsort {adding to numbers} {
	# set list {-1e+1 -1e1 -10a -2a -1 -1a -1e -1e+ 1 +1 1a +1a 1a1 1a1a 1a2 1a2a 1a10 1a10a 1aa 1e 1e+ 2a +2a 10a +10a 1e1 1e+1 -a +a +a1 +a2 +a10 a a-1 a-2 a-10 a1 a2 a10 a-a aa}
	set list {-1e+1 -1e1 -10a -2a -1 -1a -1e -1e+ 1 +1 1a +1a 1a1 1a1a 1a2 1a2a 1a10 1a10a 1aa 1e 1e+ 2a +2a 10a +10a 1e1 1e+1 -a +a +a1 +a2 +a10 a a-1 a-2 a-10 a1 a2 a10 a-a aa}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {empties, blanks and *} {
	set list [list {} 1 - + a a\t\ta a\t1\ta a\t2\ta a\t10\ta a\t*\ta a1 aa a* *]
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

test bsort {bugfix: undetermined when starting with 0} {
	# 86e1 = 860, so more than (0)445
	set list {0445a746 86e1c438 903c145d 1440c004}
	set test [bsort [list_reverse $list]]
	if {$test ne $list} {error "wrong sort"}
	checkbsort $list
} {}

test bsort {bugfix: undetermined when starting with 0 + various tests} {
	foreach list {
		{-010 -001 -01 -1 -00.1 -0.1 00 0 00.1 0.1 001 01 1 +001 +01 2 0010 010 10 011 0100 0101 01000 a-0 a-01 a-1 a0 a001 a01 a1 a2 a0010 a10}
		{-001 -01 -1 -00.1 -0.1 00 0 00.1 0.1 001 01 1 +001 +01 2 0010 010 10 011 a001 a01 a1 a2 a0010 a10}
		{{ } -001 -01 -1 -0.1 00 0 0.1 001 01 1 +001 +01 +1 2 0010 010 10 011 ! !a - a001 a01 a1 a2 a0010 a10 at1ta atta}
		{-10 -2 -1 0 1e-10 1e-1 1 1.0e2 1e2 120 1.0e3 1e3 {a -2} {a -1} a-1 a-2 a-10 de-2 de-10 e-2 e-10}
		{-2 -1 1 +1 2 +2 - -a + +a a}
		{-1e4 -10e3 -1e3 -10e2 -1000 -2e2 -200 -101 -1e+2 -1e2 -10e1 -100 -2e1 -20 -12 -1e+1 -1e1 -10 -2 -1.2 -1.10 -1.1 -1 -0.2 -1e-1 -0.10 -0.1 -0.09 -0.011 -1e-2 -0.01 -0 0 10e-10 0.01 +0.01 1e-2 0.011 +0.011 0.09 0.1 +0.1 0.10 +0.10 1e-1 +1e-1 0.18 19e-2 +19e-2 0.2 1 +1 1.1 +1.1 1.10 +1.10 1.2 1.20 2 +2 10 +10 1e1 1e+1 12 19 20 +20 2e1 +2e1 100 +100 10e1 1e2 1e+2 +1e+2 101 200 +200 2e2 1000 +1000 10e2 1e3 19e2 10e3 1e4 1e10 10e10 1e20 1e100}
		{01 1 010 10}
		{0 0a 1 1a - -a}
		{-1 -1a 0 0a 1 1a -a a}
		{0 -.1 .1}
		{-10 -1e-1}
		{-2e2 -2e1}
		{a+0.01 a+0.10}
		{a01 a10}
		{01 1 02 010}
	} {
		checkbsort $list
		set test [bsort [list_reverse $list]]
		if {$test ne $list} {error "wrong sort for: $list"}
	}
} {}

test bsort {bugfix: a00f a0e} {
	set list {a00 a00f a0e}
	set list {00f 0e 1e-1 01e 01f 1e a00 a00f a0e}
	checkbsort $list
} {}

test bsort {-chromosome} {
	bsort -sortchromosome {chr1 chr10 chr2 chrX 2 1 11 3 *}
} {chr1 1 chr2 2 3 chr10 11 chrX *}

if 0 {

	COPT="-g make ../bin/test_naturalcompare

	COPT="-g -DDEBUG=1" make ../bin/test_naturalcompare
	../bin/test_naturalcompare 00 00f 0 0e 0.1 1e-1 01e 01f 1e a00 a00f a0e 2> /dev/null | grep ' > '

	../bin/test_naturalcompare 1e-1 01e


	../bin/test_naturalcompare 0e 0f
	../bin/test_naturalcompare 01e 01f
	../bin/test_naturalcompare 01f 1e
	../bin/test_naturalcompare a01e a01f
	../bin/test_naturalcompare a1e a01f 
	../bin/test_naturalcompare a00f a0e
	../bin/test_naturalcompare 00f 0e
	../bin/test_naturalcompare 01f 1e-1
	../bin/test_naturalcompare 01f 01e
	../bin/test_naturalcompare a00f a0a
	../bin/test_naturalcompare 00     01 # on 0
	../bin/test_naturalcompare -010 -001 -01 -1 -00.1 -0.1 00 0 00.1 0.1 +001 001 +01 01 1 2 0010 010 10 011 0100 0101 01000 a-0 a-01 a-1 a0 a001 a01 a1 a2 a0010 a10 2> /dev/null | grep ' > '
	../bin/test_naturalcompare ' ' -001 -01 -1 -0.1 00 0 0.1 001 +001 01 +01 1 2 0010 010 10 011 '!' '!a' - a001 a01 a1 a2 a0010 a10 at1ta atta 2> /dev/null | grep ' > '
	../bin/test_naturalcompare -10 -2 -1 0 1e-10 1e-1 1 1.0e2 1e2 120 1.0e3 1e3 'a -2' 'a -1' a-1 a-2 a-10 de-2 de-10 e-2 e-10 2> /dev/null | grep ' > '
	../bin/test_naturalcompare -2 -1 1 +1 2 +2 - -a + +a a  2> /dev/null | grep ' > '
	../bin/test_naturalcompare -1e4 -10e3 -1e3 -10e2 -1000 -2e2 -200 -101 -1e+2 -1e2 -10e1 -100 -2e1 -20 -12 -1e+1 -1e1 -10 -2 -1.2 -1.10 -1.1 -1 -0.2 -1e-1 -0.10 -0.1 -0.09 -0.011 -1e-2 -0.01 -0 0 10e-10 0.01 +0.01 1e-2 0.011 +0.011 0.09 0.1 +0.1 0.10 +0.10 1e-1 +1e-1 0.18 19e-2 +19e-2 0.2 1 +1 1.1 +1.1 1.10 +1.10 1.2 1.20 2 +2 10 +10 1e1 1e+1 12 19 20 +20 2e1 +2e1 100 +100 10e1 1e2 1e+2 +1e+2 101 200 +200 2e2 1000 +1000 10e2 1e3 19e2 10e3 1e4 1e10 10e10 1e20 1e100 2> /dev/null | grep ' > '
# 	../bin/test_naturalcompare -1e+1 -1e1 -10a -2a -1 -1a -1e -1e+ 1 +1 1a +1a 1a1 1a1a 1a2 1a2a 1a10 1a10a 1aa 1e 1e+ 2a +2a 10a +10a 1e1 1e+1 -a +a +a1 +a2 +a10 a a-1 a-2 a-10 a1 a2 a10 a-a aa 2> /dev/null | grep ' > '

	../bin/test_naturalcompare 01 1 010 10
	../bin/test_naturalcompare - -a 0 1 0a 1a
	../bin/test_naturalcompare -1 -1a -a 0 0a 1 1a a | less
	../bin/test_naturalcompare  -.1 0 .1

	../bin/test_naturalcompare -10 -1e-1
	../bin/test_naturalcompare -2e2 -2e1

}

testsummarize
