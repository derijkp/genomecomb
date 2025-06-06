proc cg_bam2coverage {bamfile destprefix} {
	set dir [file dir $destprefix]
	set prefix [file tail $destprefix]
	file mkdir $dir
	file mkdir $destprefix.temp
	puts "coverage: making $destprefix"
	# catch_exec samtools depth -d1000000 $bamfile | cg bcol make -p pos -c contig -t su $dir/tmp-$prefix coverage
	#basequal threshold to 1 to avoid counting of primer bases
	catch_exec samtools depth -d1000000 -q 1 $bamfile | cg bcol make --header 0 --chromosomecol 0 --poscol 1 --type iu $destprefix.temp/${prefix}- 2
	file mkdir $dir/tmp
	set files [glob -nocomplain $destprefix.temp/$prefix-*.bcol]
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
		file_write $file [join $c \n]
	}
	file rename -force -- {*}[glob $destprefix.temp/$prefix*.bcol*] $dir
	file delete $destprefix.temp
	file_write $destprefix.FINISHED ""
}
