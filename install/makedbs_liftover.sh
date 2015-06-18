#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
if {[llength $argv]} {
	set dest [lindex $argv 0]
} else {
	set dest /complgen/refseq
}
set dest [file join $dest liftover]

# liftchanges
# ===========
#
file mkdir ${dest}
cd ${dest}

job_logdir log_jobs

# hg18ToHg19.over.tsv
foreach base {hg18ToHg19 hg19ToHg18 hg38ToHg19} {
	job liftchain2tsv-$base -deps {$base.over.chain} -targets {$base.over.tsv} -code {
		cg liftchain2tsv $dep > $target.temp
		file rename -force $target.temp $target
	}
}

# lift refchanges hg18Tohg19
job refchanges_hg18Tohg19 \
-deps {../hg18/genome_hg18.ifas ../hg19/genome_hg19.ifas hg18ToHg19.over.tsv} \
-targets {hg18ToHg19.over.refchanges.tsv} -code {
	cg liftfindchanges {*}$deps > $target.temp
	file rename -force $target.temp $target
}
