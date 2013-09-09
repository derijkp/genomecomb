proc cg_checksort {file} {
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [list_sub [tsv_basicfields $header 6 0] -exclude 4]
	catch {close $f}
	if {[catch {
		exec {*}[gzcat $file] $file | check_sort {*}$poss
	} result]} {
		puts stderr "error in file $file: $result"
		exit 1
	}
}
