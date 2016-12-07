proc cg_checksort {file} {
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [list_sub [tsv_basicfields $header 6 0] -exclude 4]
	gzclose $f
	if {[catch {
		exec {*}[gzcat $file] $file | check_sort {*}$poss
	} result]} {
		error "error in file $file: $result"
	}
}
