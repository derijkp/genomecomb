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

test libext {vector} {
	expr {vector("1","2,3","4,5")}
} {1,2,3,4,5}

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

test libext {vand} {
	expr {vand("0,1,1","0,0,1")}
} {0,0,1}

test libext {vor} {
	expr {vor("0,1,1","0,0,1")}
} {0,1,1}

test libext {vin} {
	expr {vin("1,2,3","4,3,2")}
} {0,1,1}

test libext {vni} {
	expr {vni("1,2,3","4,3,2")}
} {1,0,0}

test libext {vif} {
	expr {vif("1,0,0","a,a,a","1,1,0","b,b,b","c,c,c")}
} {a,b,c}

test libext {regexp} {
	list [expr {regexp("abbb","ab+")}] [expr {regexp("abbb","aB+")}] [expr {regexp("-nocase","abbb","aB*")}]
} {1 0 1}

test libext {matches} {
	list [expr {matches("abbb","a*")}] [expr {matches("abbb","A*")}] [expr {matches("-nocase","abbb","A*")}]
} {1 0 1}

test libext {q1 1,2,3,4,5,6} {
	expr {q1(1,2,3,4,5,6)}
} 2

test libext {median 1,2,3,4,5,6} {
	expr {median(1,2,3,4,5,6)}
} 3.5

test libext {q3 1,2,3,4,5,6} {
	expr {q3(1,2,3,4,5,6)}
} 5

test libext {q1 1,2,3,4,5,6} {
	expr {q1(1,2,3,4,5,6)}
} 2

test libext {median 1,2,3,4,5,6} {
	expr {median(1,2,3,4,5,6)}
} 3.5

test libext {q3 1,2,3,4,5,6} {
	expr {q3(1,2,3,4,5,6)}
} 5

test libext {q1 1,2,3,4,5,6,7} {
	expr {q1(1,2,3,4,5,6,7)}
} 2

test libext {median 1,2,3,4,5,6,7} {
	expr {median(1,2,3,4,5,6,7)}
} 4

test libext {q3 1,2,3,4,5,6,7} {
	expr {q3(1,2,3,4,5,6,7)}
} 6

test libext {q1 1,2,3,4,5,6,7,8} {
	expr {q1(1,2,3,4,5,6,7,8)}
} 2.5

test libext {median 1,2,3,4,5,6,7,8} {
	expr {median(1,2,3,4,5,6,7,8)}
} 4.5

test libext {q3 1,2,3,4,5,6,7,8} {
	expr {q3(1,2,3,4,5,6,7,8)}
} 6.5

test libext {q1 1,2,3,4,5,6,7,8,9} {
	expr {q1(1,2,3,4,5,6,7,8,9)}
} 2.5

test libext {median 1,2,3,4,5,6,7,8,9} {
	expr {median(1,2,3,4,5,6,7,8,9)}
} 5

test libext {q3 1,2,3,4,5,6,7,8,9} {
	expr {q3(1,2,3,4,5,6,7,8,9)}
} 7.5

test libext {q1 1,2,3,4,5,6,7,8,9,10} {
	expr {q1(1,2,3,4,5,6,7,8,9,10)}
} 3

test libext {median 1,2,3,4,5,6,7,8,9,10} {
	expr {median(1,2,3,4,5,6,7,8,9,10)}
} 5.5

test libext {q3 1,2,3,4,5,6,7,8,9,10} {
	expr {q3(1,2,3,4,5,6,7,8,9,10)}
} 8

test libext {q1 1} {
	expr {q1(1)}
} 1.0

test libext {median 1} {
	expr {median(1)}
} 1

test libext {q3 1} {
	expr {q3(1)}
} 1.0

test libext {q1 1,2} {
	expr {q1(1,2)}
} 1

test libext {median 1,2} {
	expr {median(1,2)}
} 1.5

test libext {q3 1,2} {
	expr {q3(1,2)}
} 2

test libext {q1 1,2,3} {
	expr {q1(1,2,3)}
} 1

test libext {median 1,2,3} {
	expr {median(1,2,3)}
} 2

test libext {q3 1,2,3} {
	expr {q3(1,2,3)}
} 3

testsummarize
