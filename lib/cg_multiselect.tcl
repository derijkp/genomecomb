proc cg_multiselect {args} {
	set query {}; set qfields {}; set sortfields {}; set oldheader {}
	set inverse 0; set group {}; set groupcols {} ; set samplingskip 0; set db {} ; set removecomment 0
	set samples {} ; set rf {}; set hc 0
	set outfile {} ; set split 1 ; set combine multicompar
	cg_options multiselect args {
		-q {set query $value}
		-qf {
			set f [gzopen $value]
			set header [tsv_open $f]
			set data [csv_file $f \t]
			gzclose $f
			set query {}
			foreach line $data {
				set el ""
				foreach field $header v $line {
					lappend el "\$$field == \"$v\""
				}
				lappend query "\( [join $el " && "] \)"
			}
			set query [join $query " || "]
		}
		-f {set qfields $value}
		-samples {set samples $value}
		-rf {set rf $value}
		-g {set group $value}
		-db {set db $value}
		-gc {lappend groupcols $value}
		-hc {set hc $value}
		-s {set sortfields $value}
		-sr {set sortfields -$value}
		-samplingskip {set samplingskip $value}
		-rc {set removecomment $value}
		-o - -outfile {
			set outfile $value
		}
		-split {
			set split $value
		}
		-combine {
			if {$value ni {multicompar cat files}} {
				error "unsupported value $value for option -combine, should be one of: multicompar cat files"
			}
			set combine $value
		}
	} {} 1
	if {$combine eq "files"} {
		if {$outfile eq ""} {error "outfile (-o) must be given when -combine files is used"}
	}
	set tempdir [tempdir]
	set resultfiles {}
	foreach file $args {
		if {$combine eq "files"} {
			set resultfile $outfile[file tail $file]
		} else {
			set resultfile $tempdir/[file tail $file]
		}
		if {$resultfile in $resultfiles} {
			error "files cannot have the same name: [file tail $file]"
		}
		cg select -q $query -f $qfields -s $sortfields -samples $samples \
			-rf $rf -g $group -db $db -gc $groupcols \
			-hc $hc -rc $removecomment \
			-s $sortfields -sr $sortfields -samplingskip $samplingskip \
			$file $resultfile
		lappend resultfiles $resultfile
	}
	if {$combine eq "multicompar"} {
		set tempfile [tempfile]
		cg multicompar -split $split $tempfile {*}$resultfiles
		if {$outfile ne ""} {
			file rename -- $tempfile $outfile
		} else {
			set f [open $tempfile]
			fcopy $f stdout
			close $f
		}
	} elseif {$combine eq "cat"} {
		if {$outfile ne ""} {
			exec cg cat -m 1 -c 0 {*}$resultfiles > $outfile
		} else {
			exec cg cat -m 1 -c 0 {*}$resultfiles >@ stdout
		}
	}
}
