# distribute to experiments
# -------------------------

proc cg_distribute {projectfile srcdir destdir {subdir {}} {postfix {}}} {
	set srcdir [file_absolute $srcdir]
	set destdir [file_absolute $destdir]
	set c [split [string trim [file_read $projectfile]] \n]
	set header [split [list_shift c] \t]
	if {[lrange $header 0 3] ne "project experiment sample barcode"} {
		error "project-$project.txt must have columns: project experiment sample barcode"
	}
	list_foreach {p experiment sample barcode} $c {
		if {$subdir ne ""} {
			set dest $destdir/$experiment/$sample$postfix/$subdir
		} else {
			set dest $destdir/$experiment/$sample$postfix
		}
		file mkdir $dest
		set files [glob -nocomplain $srcdir/*$barcode]
		set poss [list_find -regexp $files "^.*\[^0-9\]$barcode\$"]
		set files [list_sub $files $poss]
		foreach file $files {
			putslog "hardlinking $file/* to $dest"
			exec cp -al {*}[glob $file/*] $dest
		}
	}

}
