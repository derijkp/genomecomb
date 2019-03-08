#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test compress {basic combinations} {
	set result {}
	foreach file {sample.bed annot_compar-exomes_yri_parts.tsv} {
		foreach startmethod {{} lz4 rz bgz bz2} {
			foreach method {lz4 rz bgz bz2} {
				test_cleantmp
				file copy -force data/$file tmp/test
				if {$startmethod ne ""} {
					cg $startmethod tmp/test
					set src tmp/test.$startmethod
				} else {
					set src tmp/test
				}
				append result [list $startmethod [file size $src] $method]
				cg $method -stack 1 $src
				append result " [file size tmp/test.$method] $file - [lsort -dict [glob tmp/*]]\n"
				cg zcat tmp/test.$method > tmp/test.res
				exec diff tmp/test.res data/$file
			}
		}
	}
	set result
} {{} 521 lz4 357 sample.bed - tmp/test.lz4
{} 521 rz 287 sample.bed - tmp/test.rz
{} 521 bgz 294 sample.bed - tmp/test.bgz
{} 521 bz2 281 sample.bed - tmp/test.bz2
lz4 357 lz4 357 sample.bed - tmp/test.lz4
lz4 357 rz 287 sample.bed - tmp/test.rz
lz4 357 bgz 294 sample.bed - tmp/test.bgz
lz4 357 bz2 281 sample.bed - tmp/test.bz2
rz 287 lz4 357 sample.bed - tmp/test.lz4
rz 287 rz 287 sample.bed - tmp/test.rz
rz 287 bgz 294 sample.bed - tmp/test.bgz
rz 287 bz2 281 sample.bed - tmp/test.bz2
bgz 294 lz4 357 sample.bed - tmp/test.lz4
bgz 294 rz 287 sample.bed - tmp/test.rz
bgz 294 bgz 294 sample.bed - tmp/test.bgz
bgz 294 bz2 281 sample.bed - tmp/test.bz2
bz2 281 lz4 357 sample.bed - tmp/test.lz4
bz2 281 rz 287 sample.bed - tmp/test.rz
bz2 281 bgz 294 sample.bed - tmp/test.bgz
bz2 281 bz2 281 sample.bed - tmp/test.bz2
{} 618974 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4
{} 618974 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.rz
{} 618974 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz
{} 618974 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2
lz4 209708 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4
lz4 209708 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.rz
lz4 209708 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz
lz4 209708 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2
rz 178075 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4
rz 178075 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.rz
rz 178075 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz
rz 178075 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2
bgz 173897 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4
bgz 173897 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.rz
bgz 173897 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz
bgz 173897 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2
bz2 133936 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4
bz2 133936 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.rz
bz2 133936 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz
bz2 133936 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2
}

test compress {basic combinations -keep 1 -index 1} {
	set result {}
	foreach file {sample.bed annot_compar-exomes_yri_parts.tsv} {
		foreach startmethod {{} lz4 rz bgz bz2} {
			foreach method {lz4 rz bgz bz2} {
				test_cleantmp
				file copy -force data/$file tmp/test
				if {$startmethod ne ""} {
					cg $startmethod tmp/test
					set src tmp/test.$startmethod
				} else {
					set src tmp/test
				}
				append result [list $startmethod [file size $src] $method]
				cg $method -keep 1 -index 1 $src
				append result " [file size tmp/test.$method] $file - [lsort -dict [glob tmp/*]]\n"
				cg zcat tmp/test.$method > tmp/test.res
				exec diff tmp/test.res data/$file
			}
		}
	}
	set result
} {{} 521 lz4 357 sample.bed - tmp/test tmp/test.lz4 tmp/test.lz4.lz4i
{} 521 rz 287 sample.bed - tmp/test tmp/test.rz
{} 521 bgz 294 sample.bed - tmp/test tmp/test.bgz
{} 521 bz2 281 sample.bed - tmp/test tmp/test.bz2
lz4 357 lz4 357 sample.bed - tmp/test.lz4 tmp/test.lz4.lz4i
lz4 357 rz 287 sample.bed - tmp/test.lz4 tmp/test.rz
lz4 357 bgz 294 sample.bed - tmp/test.bgz tmp/test.lz4
lz4 357 bz2 281 sample.bed - tmp/test.bz2 tmp/test.lz4
rz 287 lz4 357 sample.bed - tmp/test.lz4 tmp/test.lz4.lz4i tmp/test.rz
rz 287 rz 287 sample.bed - tmp/test.rz
rz 287 bgz 294 sample.bed - tmp/test.bgz tmp/test.rz
rz 287 bz2 281 sample.bed - tmp/test.bz2 tmp/test.rz
bgz 294 lz4 357 sample.bed - tmp/test.bgz tmp/test.lz4 tmp/test.lz4.lz4i
bgz 294 rz 287 sample.bed - tmp/test.bgz tmp/test.rz
bgz 294 bgz 294 sample.bed - tmp/test.bgz
bgz 294 bz2 281 sample.bed - tmp/test.bgz tmp/test.bz2
bz2 281 lz4 357 sample.bed - tmp/test.bz2 tmp/test.lz4 tmp/test.lz4.lz4i
bz2 281 rz 287 sample.bed - tmp/test.bz2 tmp/test.rz
bz2 281 bgz 294 sample.bed - tmp/test.bgz tmp/test.bz2
bz2 281 bz2 281 sample.bed - tmp/test.bz2
{} 618974 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test tmp/test.lz4 tmp/test.lz4.lz4i
{} 618974 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test tmp/test.rz
{} 618974 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test tmp/test.bgz
{} 618974 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test tmp/test.bz2
lz4 209708 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4 tmp/test.lz4.lz4i
lz4 209708 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4 tmp/test.rz
lz4 209708 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz tmp/test.lz4
lz4 209708 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2 tmp/test.lz4
rz 178075 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.lz4 tmp/test.lz4.lz4i tmp/test.rz
rz 178075 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.rz
rz 178075 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz tmp/test.rz
rz 178075 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2 tmp/test.rz
bgz 173897 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz tmp/test.lz4 tmp/test.lz4.lz4i
bgz 173897 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz tmp/test.rz
bgz 173897 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz
bgz 173897 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz tmp/test.bz2
bz2 133936 lz4 209708 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2 tmp/test.lz4 tmp/test.lz4.lz4i
bz2 133936 rz 178075 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2 tmp/test.rz
bz2 133936 bgz 173897 annot_compar-exomes_yri_parts.tsv - tmp/test.bgz tmp/test.bz2
bz2 133936 bz2 133936 annot_compar-exomes_yri_parts.tsv - tmp/test.bz2
}

test compress {multiple} {
	set result {}
	foreach method {lz4 rz bgz bz2} {
		test_cleantmp
		file copy -force data/sample.bed tmp/test1
		file copy -force data/annot_compar-exomes_yri_parts.tsv tmp/test2
		cg $method tmp/test1
		cg $method tmp/test2
		cg zcat tmp/test1.$method > tmp/test1-2
		cg zcat tmp/test2.$method > tmp/test2-2
		exec diff tmp/test1-2 data/sample.bed
		exec diff tmp/test2-2 data/annot_compar-exomes_yri_parts.tsv
		append result "$method [file size tmp/test1.$method] [file size tmp/test2.$method]\n"
	}
	set result
} {lz4 357 209708
rz 287 178075
bgz 294 173897
bz2 281 133936
}

test compress {-o} {
	set result {}
	foreach file {sample.bed annot_compar-exomes_yri_parts.tsv} {
		foreach method {lz4 rz bgz bz2} {
			test_cleantmp
			file copy -force data/$file tmp/test
			cg $method -o tmp/new.$method tmp/test
			cg zcat tmp/new.$method > tmp/test2
			exec diff --brief tmp/test2 data/$file
			append result "$method [file size data/$file] [file size tmp/new.$method] $file\n"
		}
	}
	set result
} {lz4 521 357 sample.bed
rz 521 287 sample.bed
bgz 521 294 sample.bed
bz2 521 281 sample.bed
lz4 618974 209708 annot_compar-exomes_yri_parts.tsv
rz 618974 178075 annot_compar-exomes_yri_parts.tsv
bgz 618974 173897 annot_compar-exomes_yri_parts.tsv
bz2 618974 133936 annot_compar-exomes_yri_parts.tsv
}

test compress {-o from stdin} {
	set result {}
	foreach file {sample.bed annot_compar-exomes_yri_parts.tsv} {
		foreach method {lz4 rz bgz bz2} {
			test_cleantmp
			file copy -force data/$file tmp/test
			exec cg $method -o tmp/new.$method < tmp/test
			cg zcat tmp/new.$method > tmp/test2
			exec diff tmp/test2 data/$file
			append result "$method [file size data/$file] [file size tmp/new.$method] $file\n"
		}
	}
	set result
} {lz4 521 357 sample.bed
rz 521 287 sample.bed
bgz 521 294 sample.bed
bz2 521 281 sample.bed
lz4 618974 209708 annot_compar-exomes_yri_parts.tsv
rz 618974 178075 annot_compar-exomes_yri_parts.tsv
bgz 618974 173897 annot_compar-exomes_yri_parts.tsv
bz2 618974 133936 annot_compar-exomes_yri_parts.tsv
}

test compress {in pipe} {
	set result {}
	foreach file {sample.bed annot_compar-exomes_yri_parts.tsv} {
		foreach method {lz4 rz bgz bz2} {
			test_cleantmp
			file copy -force data/$file tmp/test
			exec cat tmp/test | cg $method > tmp/new.$method
			cg zcat tmp/new.$method > tmp/test2
			exec diff tmp/test2 data/$file
			append result "$method [file size data/$file] [file size tmp/new.$method] $file\n"
		}
	}
	set result
} {lz4 521 357 sample.bed
rz 521 287 sample.bed
bgz 521 294 sample.bed
bz2 521 281 sample.bed
lz4 618974 209708 annot_compar-exomes_yri_parts.tsv
rz 618974 178075 annot_compar-exomes_yri_parts.tsv
bgz 618974 173897 annot_compar-exomes_yri_parts.tsv
bz2 618974 133936 annot_compar-exomes_yri_parts.tsv
}

test compress {compress command} {
	set result {}
	foreach file {sample.bed} {
		foreach startmethod {{} lz4 rz bgz bz2} {
			foreach method {lz4 rz bgz bz2} {
				test_cleantmp
				file copy -force data/$file tmp/test
				if {$startmethod ne ""} {
					cg $startmethod tmp/test
					set src tmp/test.$startmethod
				} else {
					set src tmp/test
				}
				append result [list $startmethod [file size $src] $method]
				set destfile tmp/test.$method
				compress $src $destfile 0 0 1 1
				append result " [file size $destfile] $file - [lsort -dict [glob tmp/*]]\n"
				cg zcat $destfile > tmp/test.res
				exec diff tmp/test.res data/$file
			}
		}
	}
	set result
} {{} 521 lz4 364 sample.bed - tmp/test.lz4
{} 521 rz 287 sample.bed - tmp/test.rz
{} 521 bgz 294 sample.bed - tmp/test.bgz
{} 521 bz2 281 sample.bed - tmp/test.bz2
lz4 357 lz4 364 sample.bed - tmp/test.lz4
lz4 357 rz 287 sample.bed - tmp/test.rz
lz4 357 bgz 294 sample.bed - tmp/test.bgz
lz4 357 bz2 281 sample.bed - tmp/test.bz2
rz 287 lz4 364 sample.bed - tmp/test.lz4
rz 287 rz 287 sample.bed - tmp/test.rz
rz 287 bgz 294 sample.bed - tmp/test.bgz
rz 287 bz2 281 sample.bed - tmp/test.bz2
bgz 294 lz4 364 sample.bed - tmp/test.lz4
bgz 294 rz 287 sample.bed - tmp/test.rz
bgz 294 bgz 294 sample.bed - tmp/test.bgz
bgz 294 bz2 281 sample.bed - tmp/test.bz2
bz2 281 lz4 364 sample.bed - tmp/test.lz4
bz2 281 rz 287 sample.bed - tmp/test.rz
bz2 281 bgz 294 sample.bed - tmp/test.bgz
bz2 281 bz2 281 sample.bed - tmp/test.bz2
}

testsummarize
