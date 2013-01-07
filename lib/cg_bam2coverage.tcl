proc cg_bam2coverage {bamfile destprefix} {
	set dir [file dir $destprefix]
	set prefix [file tail $destprefix]
	file mkdir $dir
	puts "coverage: making $destprefix"
	exec samtools depth $bamfile | cg bcol make -p pos -c contig -t su $dir/tmp-$prefix coverage
	file mkdir $dir/tmp
	set files [glob -nocomplain $dir/tmp/$prefix-*.bcol]
	foreach file $files {
		# samtools depth is 1 based, this hack will correct this (without having to remap the entire coverage file)
		set c [split [string trim [file_read $file]] \n]
		set line [lindex $c end-1]
		set temp [lindex $line 0]
		lset line 0 [incr temp -1]
		lset c end-1 [join $line \t]
		set line [lindex $c end]
		set temp [lindex $line 0]
		lset line 0 [incr temp -1]
		lset c end [join $line \t]
		write_file $file [join $c \n]
	}
	file rename {*}[glob $dir/tmp/$prefix-*.bcol*] $dir
	write_file $destprefix.FINISHED ""
}
