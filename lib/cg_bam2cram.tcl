proc cg_bam2cram {args} {
	set bamfile -
	set cramfile -
	set index 1
	set refseq {}
	set threads 1
	set handlebam {}
	set links ignore
	cg_options cg_bam2cram args {
		-refseq {
			set refseq [refseq $value]
		}
		-handlebam {
			if {$value ni "rm old keep"} {
				error "wrong value for -handlebam, must be one of: rm, old, keep"
			}
			set handlebam $value
		}
		-links {
			if {$value ni "ignore rename convert error"} {
				error "wrong value for -links, must be one of: ignore rename convert error"
			}
			set links $value
		}
		-threads {set threads $value}
		-index {set index $value}
	} {bamfile cramfile} 1
	set bamfile [file_absolute $bamfile]
	set outformat cram
	if {$cramfile eq "-"} {
		set cramfile [file root $bamfile].cram
		if {$handlebam eq ""} {set handlebam old}
	} else {
		if {$handlebam eq ""} {set handlebam keep}
	}
	set cramfile [file_absolute $cramfile]
	if {![catch {file link $bamfile} link]} {
		if {$links eq "rename"} {
			set link [file root $link].cram
			if {[file pathtype $link] ne "absolute"} {
				set link [file dir $bamfile]/$link
			}
			mklink $link $cramfile
			mklink $link.crai $cramfile.crai
			if {[file exists $link.analysisinfo]} {
				mklink $link.analysisinfo $cramfile.analysisinfo
			}
			if {![catch {file lstat $bamfile a}]} {
				set mtime $a(mtime)
				exec touch -h -d [clock format $mtime] $cramfile
			}
			if {![catch {file lstat $bamfile.bai a}]} {
				set mtime $a(mtime)
				exec touch -h -d [clock format $mtime] $cramfile.crai
			}
			if {$handlebam eq "old"} {
				catch {file rename -force $bamfile $bamfile.old}
				catch {file rename -force $bamfile.bai $bamfile.bai.old}
			} elseif {$handlebam eq "rm"} {
				file delete $bamfile
			}
			return
		} elseif {$links eq "ignore"} {
			return
		} elseif {$links eq "error"} {
			error "$bam is a symlink"
		}
	}
	exec samtools view --threads $threads --no-PG -h -C -T $refseq --no-PG $bamfile > $cramfile.temp.cram
	exec samtools index $cramfile.temp.cram
	file rename -force $cramfile.temp.cram.crai $cramfile.crai
	file rename -force $cramfile.temp.cram $cramfile
	exec touch -r $bamfile $cramfile
	if {$handlebam eq "old"} {
		catch {file rename -force $bamfile $bamfile.old}
		catch {file rename -force $bamfile.bai $bamfile.bai.old}
	} elseif {$handlebam eq "rm"} {
		file delete $bamfile
	}
}
