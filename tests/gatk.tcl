#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test gatkexec {version} {
	regexp {^4\.[0-9]\.[0-9]\.[0-9]$} [gatkexec version]
} 1

test gatkexec {error} {
	gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller -R bla -I bla -O bla
} {*The specified fasta file * does not exist.*} error match

test gatkexec {error -finishedpattern} {
	catch {
		gatkexec -finishedpattern {does not exist} {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller -R bla -I bla -O bla
	} msg
} 0

testsummarize

