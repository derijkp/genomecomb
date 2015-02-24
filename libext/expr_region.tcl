proc tcl::mathfunc::samechr {chr1 chr2} {
	set chr1 [chr_clip $chr1]
	set chr2 [chr_clip $chr2]
	expr {$chr1 eq $chr2}
}

proc tcl::mathfunc::chr_clip {chr1} {
	::chr_clip $chr1
}
