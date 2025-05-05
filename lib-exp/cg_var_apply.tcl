proc var_apply_job {args} {
	upvar job_logdir job_logdir
	#
	set checkchromosome 1
	set skipunknowntype 0
	cg_options var_apply args {
		-checkchromosome {set checkchromosome $value}
		-skipunknowntype {set skipunknowntype $value}
	} {sequencefile varfile resultfile} 3 3
	set sequencefile [file_absolute $sequencefile]
	set varfile [file_absolute $varfile]
	set resultfile [file_absolute $resultfile]

	catch {close $f} ; catch {close $fv}
	set f [gzopen $sequencefile]
	set name [gets $f]
	regsub ^> $name {} name
	set seq [gets $f]
	gzclose $f
	# todo check for correct fields in this file
	set fv [open "[list | cg select -s {chromosome -begin -end} $varfile]"]
	set header [tsv_open $fv]
	set poss [tsv_basicfields $header 6]
	while 1 {
		if {[gets $fv line] == -1} break
		foreach {chromosome begin end type ref alt} [list_sub [split $line \t] $poss] break
		if {$checkchromosome && $chromosome ne $name} {
			error "chromosome \"$chromosome\" in varfile $varfile does not match name of sequence in FASTA file $sequencefile"
		}
		switch $type {
			snp {
				set seq [string replace $seq $begin [expr {$end-1}] $alt]
			}
			del {
				set seq [string replace $seq $begin [expr {$end-1}] $alt]
			}
			ins {
				set seq [string range $seq 0 [expr {$begin-1}]]$alt[string range $seq $begin end]
			}
			sub {
				set seq [string replace $seq $begin [expr {$end-1}] $alt]
			}
			default {
				if {$skipunknowntype} {
					puts stderr "skipping unsuported $type variant: $line"
				} else {
					error "unsuported $type variant: $line"
				}
			}
		}
	}
	close $fv
	set o [wgzopen $resultfile]
	puts $o >$name
	puts $o $seq
	gzclose $o
}

proc cg_var_apply {args} {
	set args [job_init {*}$args]
	var_apply_job {*}$args
	job_wait
}

