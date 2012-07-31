#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test libext {lmin} {
	tcl::mathfunc::lmin 1,10 3,10
} 1

test libext {lmin} {
	tcl::mathfunc::lmin 1,- 3,10
} 1

test libext {lminpos} {
	tcl::mathfunc::lminpos 1,10 3,10
} 0

test libext {lminpos} {
	tcl::mathfunc::lminpos 1,10 3,0
} 1

test libext {lminpos} {
	tcl::mathfunc::lminpos -,- -,-
} NaN

test libext {lmax} {
	tcl::mathfunc::lmax 1,10 3,10
} 10

test libext {lmax} {
	tcl::mathfunc::lmax 1,- 3,10
} 10

test libext {lmaxpos} {
	tcl::mathfunc::lmaxpos 1,10 30,10
} 0

test libext {lmaxpos} {
	tcl::mathfunc::lmaxpos 1,10 3,0
} 1

test libext {lmaxpos} {
	tcl::mathfunc::lmaxpos -,- -,-
} NaN

test libext {lindex} {
	tcl::mathfunc::lindex 1,10,2 1
} 10

test libext {lindex} {
	expr {lindex("1,10,2",1)}
} 10

test libext {lone} {
	expr {lone("0,0,0")}
} 0

test libext {lone} {
	expr {lone("0,0,1")}
} 1

test libext {lall} {
	expr {lall("0,0,0")}
} 0

test libext {lall} {
	expr {lall("0,0,1")}
} 0

test libext {lall} {
	expr {lall("1,1,1")}
} 1

test libext {lall} {
	expr {lall("1,1,1","1,1,0")}
} 0

test libext {lall} {
	expr {lall("1,1,1","1,1,-")}
} 1

test libext {lall} {
	expr {lall("1,1,1","1,1,1")}
} 1

test libext {lcount} {
	expr {lcount("1,0,0","1,1,1")}
} 4

test libext {llen} {
	expr {llen("1,0,0","1,1,1")}
} 6

test libext {contains} {
	expr {contains("1,2,3","2")}
} 1

test libext {contains} {
	expr {contains("1,2,3","4")}
} 0

test libext {shares 1} {
	expr {shares("1,2,3","2")}
} 1

test libext {shares 2} {
	expr {shares("1,2,3","4,2")}
} 1

test libext {shares 3} {
	expr {shares("1,2,3","4")}
} 0

test libext {shares 4} {
	expr {shares("1,2,3","4,5")}
} 0

test libext {vabs} {
	tcl::mathfunc::vabs 1,-1,2
} {1,1,2}

test libext {vavg} {
	tcl::mathfunc::vavg 1,10 3,10
} {2.0,10.0}

test libext {vavg -} {
	tcl::mathfunc::vavg 1,- 2,-
} {1.5,NaN}

test libext {vmax} {
	expr {vmax("10,2","1,20")}
} {10,20}

test libext {vmin} {
	expr {vmin("10,2","1,20")}
} {1,2}

testsummarize
