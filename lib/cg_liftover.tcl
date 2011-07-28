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
	set unmappedfile $resultfile.unmapped
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set line [split [gets $f] \t]
	if {![regexp ^chr [lindex $line 0]]} {set addchr 1} else {set addchr 0}
	close $f
	set poss [tsv_basicfields $header 3]
	if {$poss ne "0 1 2"} {error "Rearranged header not supported yet"}
	if {$addchr} {
		set fields [cg select -h $varfile]
		set fields [lreplace $fields 0 0 "chromosome=\"chr\" \$chromosome"]
		cg select -f $fields -sh $resultfile.temph $varfile $resultfile.temp
	} else {
		cg select -sh $resultfile.temph $varfile $resultfile.temp
	}
	# set dir [file dir [exec which liftOver]]
	if {[catch {exec liftOver -bedPlus=3 -tab $resultfile.temp $liftoverfile $resultfile.temp2 $unmappedfile} errmsg]} {
		puts "$errmsg"
	}
	exec cat $resultfile.temph $resultfile.temp2 > $resultfile
#	file delete $resultfile.temph $resultfile.temp $resultfile.temp2
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
