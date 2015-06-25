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

# get chain files
foreach {base src} {hg18ToHg19 hg18 hg19ToHg18 hg19 hg38ToHg19 hg38 hg19ToHg38 hg19} {
	job getchain-$base -vars {src base} -targets {$base.over.chain} -code {
		wgetfile http://hgdownload.soe.ucsc.edu/goldenPath/$src/liftOver/$base.over.chain.gz
	}
}


# make liftover files
foreach base {hg18ToHg19 hg19ToHg18 hg38ToHg19 hg19ToHg38} {
	job liftchain2tsv-$base -deps {$base.over.chain} -targets {$base.over.tsv} -code {
		cg liftchain2tsv $dep > $target.temp
		file rename -force $target.temp $target
	}
}

# lift refchanges hg18Tohg19
foreach {base src dest} {
	hg18ToHg19 hg18 hg19
	hg19ToHg18 hg19 hg18
	hg38ToHg19 hg38 hg19
	hg19ToHg38 hg19 hg38
} {
	job refchanges_$base \
	-deps {../$src/genome_$src.ifas ../$dest/genome_$dest.ifas $base.over.tsv} \
	-targets {$base.over.refchanges.tsv} -code {
		cg liftfindchanges {*}$deps > $target.temp
		file rename -force $target.temp $target
	}
}
