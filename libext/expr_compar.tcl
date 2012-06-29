proc tcl::mathfunc::samegeno {a11 a12 a21 a22} {
	if {($a11 == $a21) && ($a12 == $a22)} {return 1}
	if {($a11 == $a22) && ($a12 == $a21)} {return 1}
	return 0
}
