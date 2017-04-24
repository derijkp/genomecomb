proc fastqc {target args} {
	set dir [filetemp $target]
	file delete $dir
	file mkdir $dir
	exec fastqc -o $dir --extract {*}$args
	file rename [glob $dir/*_fastqc] $target
	file delete -force $dir
}
