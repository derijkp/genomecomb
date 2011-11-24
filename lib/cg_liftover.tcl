proc cg_liftover {args} {
#	set pos 0
#	set sumfields {}
#	foreach {key value} $args {
#		switch -- $key {
#			-sumfields {
#				set sumfields $value
#			}
#			-- break
#			default {
#				break
#			}
#		}
#		incr pos 2
#	}
#	set args [lrange $args $pos end]
	if {([llength $args] < 3)} {
		errorformat liftover
		exit 1
	}
	foreach {varfile liftoverfile resultfile} $args break
	if {[file exists $resultfile]} {
		error "file $resultfile already exists"
	}
	set unmappedfile $resultfile.unmapped
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set line [split [gets $f] \t]
	if {![regexp ^chr [lindex $line 0]]} {set addchr 1} else {set addchr 0}
	close $f
	set poss [tsv_basicfields $header 6]
	if {$poss ne "0 1 2 3 4 5"} {error "Rearranged header not supported yet, start header should be: chromosome begin end type ref alt"}
	#
	# make input file ($resultfile.temp) for liftover
	#
	set fields [lrange $header 0 5]
	set id {}
	foreach field $fields {
		lappend id \$$field
	}
	if {$addchr} {
		set fields [lreplace $fields 0 0 "chromosome=\"chr\" \$chromosome"]
	}
	lappend fields id=[join $id { "-" }]
	cg select -f $fields -sh $resultfile.temph $varfile $resultfile.temp
	#
	# do liftover -> $resultfile.temp2
	#
	# set dir [file dir [exec which liftOver]]
	if {[catch {exec liftOver -bedPlus=3 -tab $resultfile.temp $liftoverfile $resultfile.temp2 $unmappedfile} errmsg]} {
		puts "$errmsg"
	}
	#
#	# sort liftovered file ion id -> $resultfile.temp3
#	#
#	file_write $resultfile.temph [join {chrom begin end type ref alt id} \t]\n
#	cg select -s {chrom start end type alt} -hf $resultfile.temph $resultfile.temp2 $resultfile.temp3
	#
	# add original data back to liftovered file -> $resultfile.temp3
	#
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set fl [open $resultfile.temp2]
	set temp [tsv_open $fl]
	set o [open $resultfile.temp3 w]
	lappend header beforeliftover
	puts $o [join $header \t]
	while {![eof $fl]} {
		set lline [split [gets $fl] \t]
		if {![llength $lline]} continue
		set lname [lindex $lline 6]
		while 1 {
			set line [split [gets $f] \t]
			set name [join [list_sub $line $poss] -]
			if {$name eq $lname} break
			if {[eof $f]} {error "$lname not found"}
		}
		set rline [lrange $lline 0 5]
		lappend rline {*}[lrange $line 6 end] [join [lrange [split $lname -] 0 2] -]
		puts $o [join $rline \t]
	}
	close $o
	close $fl
	close $f
	#
	# sort result -> $resultfile.temp4
	#
	cg select -s [lrange $header 0 5] $resultfile.temp3 $resultfile.temp4
	#
	# rename result, cleanup
	#
	file rename $resultfile.temp4 $resultfile
	file delete $resultfile.temph $resultfile.temp $resultfile.temp2 $resultfile.temp3
}

if 0 {
	cd ~/dev/completegenomics/tests
	# liftOver -bedPlus=3 -tab data/vars1.sft ~/bin/hg18ToHg19.over.chain temp.sft temp.unmapped
	set varfile data/vars1.sft
	set liftoverfile /home/peter/bin/hg18ToHg19.over.chain
	set resultfile temp.sft
	set unmappedfile temp.unmapped
	cg liftover data/vars1.sft /home/peter/bin/hg18ToHg19.over.chain temp.sft temp.unmapped
}
