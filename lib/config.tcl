proc configdir {{app cg}} {
	global configdata
	if {[info exists configdata(configdir)]} {return $configdata(configdir)}
	if [info exists env(APPDATA)] {
		set basedir $env(APPDATA)
	} elseif [info exists env(HOME)] {
		set basedir $env(HOME)
	} else {
		set basedir [file_absolute ~]
	}
	set configdir $basedir/.$app
	set makedir_err [catch {file mkdir $configdir}]
	if {!$makedir_err} {return $configdir}
	if {($::tcl_platform(platform) eq "windows") && $makedir_err} {
		package require registry
		set mydocuments [registry get "HKEY_CURRENT_USER\\Software\\Microsoft\\Windows\\CurrentVersion\\Explorer\\User Shell Folders" Personal]
		string_change $mydocuments [list $ \\$ \\ \\\\]
		regsub -all {%([^%]+)%} $mydocuments {$env(\1)} mydocuments
		set mydocuments [subst -nobackslashes -nocommands $mydocuments]
		set basedir [file join $mydocuments {Application Data}]
		set configdir $basedir/.$app
		set makedir_err [catch {file mkdir $configdir}]
	}
	set configdata(configdir) $configdir
}
