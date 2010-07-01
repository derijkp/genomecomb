package require Extral

proc covered regfile {
	set f [open $regfile]
	set poss [open_region $f]
	close $f
	exec covered $regfile {*}$poss >@ stdout 2>@ stderr
}
