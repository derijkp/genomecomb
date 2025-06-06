proc cg_error_report args {
	set logfile {}
	set format 2
	cg_options error_report args {
		-format {set format $value}
	} logfile 0 1
	if {[file isdir $logfile]} {
		set dir [file_absolute $logfile]
		set logfile ""
	} else {
		set dir [pwd]
	}
	if {$logfile eq ""} {
		set logfiles [bsort [glob -nocomplain $dir/process_*.*.finished $dir/process_*.*.running $dir/process_*.*.error $dir/process_*.*.submitting]]
		if {![llength $logfiles]} {
			set logfiles [bsort [glob -nocomplain $dir/*.*.finished $dir/*.*.running $dir/*.*.error $dir/*.*.submitting]]
		}
		set logfile [lindex $logfiles end]
	}
	set ext [file extension $logfile]
	if {$ext eq ".submitting"} {
		set num [lindex [cg select -g all $logfile] end]
		puts "experiment still begin submitted ($num jobs submitted)"
		puts "logfile: $logfile"
	} elseif {$ext eq ".running"} {
		puts "experiment still running"
		puts "logfile: $logfile"
	} elseif {$ext eq ".finished"} {
		set num [lindex [cg select -g all $logfile] end]
		puts "experiment finished successfully ($num jobs)"
		puts "logfile: $logfile"
	} else {
		if {$format eq "1"} {
			set bold ""
			set underline ""
			set green ""
			set yellow ""
			set cyan ""
			set normal ""
			set foutput [open "| less -R" w]
		} elseif {$format eq "2"} {
			set bold "\033\[1;1m"
			set underline "\033\[1;4m"
			set green "\033\[1;32m"
			set yellow "\033\[1;33m"
			set cyan "\033\[1;36m"
			set normal "\033\[0m"
			set foutput [open "| less -R" w]
		} else {
			set bold ""
			set underline ""
			set green ""
			set yellow ""
			set cyan ""
			set normal ""
			set foutput stdout
		}
		set num [lindex [cg select -g all $logfile] end]
		set nrerrors [lindex [cg select -g all -q {$status eq "error"} $logfile] end]
		puts $foutput "logfile: $logfile"
		puts $foutput "Analysis ended with $nrerrors errors (out of $num jobs)"
		flush $foutput
		set f [gzopen $logfile]
		set header [tsv_open $f]
		set poss [list_cor $header {job jobid status submittime starttime endtime duration time_seconds targets msg run cores pct_cpu maxmem bytes_read bytes_written}]
		
		while 1 {
			if {[gets $f line] == -1} break
			foreach {job jobid status submittime starttime endtime duration time_seconds targets msg run cores pct_cpu maxmem bytes_read bytes_written} [split $line \t] break
			if {$status ne "error"} continue
			regsub -all {\\n} $msg \n omsg
			regsub -all {\\t} $omsg \t omsg
			set pre [file dir $job]
			set post [file tail $job]

			puts $foutput ""
			set title "job $jobid ($job)"
			puts $foutput $green$title$normal
			puts $foutput $green[string_fill - [string length $title]]$normal
			puts $foutput ""
			puts $foutput "${yellow}error message ($pre/err/$post.err):$normal"
			puts $foutput $omsg
			puts $foutput ""
			puts $foutput "${yellow}time: $starttime - $endtime ($duration)$normal"
			puts $foutput "${yellow}maxmem: $maxmem$normal"
			puts $foutput "${yellow}run_file: $pre/run/$post.run$normal"
			flush $foutput
		}
		close $f
		if {$format in "1 2"} {
			catch {close $foutput}
		}
	}
}
