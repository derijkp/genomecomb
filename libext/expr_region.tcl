proc tcl::mathfunc::samechr {chr1 chr2} {
	regsub {^[Cc]hr} $chr1 {} chr1
	regsub {^[Cc]hr} $chr2 {} chr2
	expr {$chr1 eq $chr2}
}
