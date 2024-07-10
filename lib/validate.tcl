proc validate_distrreg {distrreg refseq} {
# putsvars validate_distrreg distrreg refseq
	putslog "validating distrreg $distrreg $refseq"
	if {$distrreg eq ""} return
	distrreg_regs $distrreg $refseq
}

proc validate_map {cmd refseq preset} {
# putsvars validate_map cmd refseq preset
	putslog "validating map $cmd $refseq $preset"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_map_${cmd}]} {
		validate_map_${cmd} $refseq $preset
		return
	}
	if {[catch {exec which $cmd}]} {
		error "command \"$cmd\" not available, make sure it is installed, e.g. using \"cg install\""
	}
	if {[auto_load refseq_${cmd}]} {
		refseq_${cmd} $refseq $preset
	}
}

proc validate_var {cmd refseq distrreg datatype} {
# putsvars validate_var cmd refseq distrreg datatype
	putslog "validating var $cmd $refseq $distrreg $datatype"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_var_${cmd}]} {
		validate_var_${cmd} $refseq $distrreg $datatype
		return
	}
	if {[catch {exec which $cmd}]} {
		error "command \"$cmd\" not available, make sure it is installed, e.g. using \"cg install\""
	}
	
}

proc validate_sv {cmd refseq distrreg} {
# putsvars validate_sv cmd refseq distrreg
	putslog "validating sv $cmd $refseq $distrreg"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_sv_${cmd}]} {
		validate_sv_${cmd} $refseq $distrreg
		return
	}
	if {[catch {exec which $cmd}]} {
		error "command \"$cmd\" not available, make sure it is installed, e.g. using \"cg install\""
	}
	
}

proc validate_meth {cmd refseq preset distrreg} {
# putsvars validate_meth cmd refseq preset distrreg
	putslog "validationg meth $cmd $refseq $distrreg"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_meth_${cmd}]} {
		validate_meth_${cmd} $refseq $preset $distrreg
		return
	}
	if {[catch {exec which $cmd}]} {
		error "command \"$cmd\" not available, make sure it is installed, e.g. using \"cg install\""
	}

}

proc validate_count {cmd refseq} {
# putsvars validate_count cmd refseq
	putslog "validating count $cmd $refseq"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_count_${cmd}]} {
		validate_count_${cmd} $refseq
		return
	}
	if {[catch {exec which $cmd}]} {
		error "command \"$cmd\" not available, make sure it is installed, e.g. using \"cg install\""
	}
	
}

proc validate_iso {cmd refseq preset reftranscripts organelles distrreg} {
# putsvars validate_iso cmd refseq preset reftranscripts organelles distrreg
	putslog "validating iso $cmd $refseq $reftranscripts $organelles $distrreg"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_iso_${cmd}]} {
		validate_iso_${cmd} $refseq $preset $reftranscripts $organelles $distrreg
		return
	}
	if {[catch {exec which $cmd}]} {
		error "command \"$cmd\" not available, make sure it is installed, e.g. using \"cg install\""
	}
	if {$reftranscripts eq ""}  {
		ref_gtftranscripts $refseq
	} elseif {$reftranscripts in {none -}} {
		# ignore
	} elseif {![file exists $reftranscripts]} {
		error "reference transcripts file does not exist ($reftranscripts)"
	}
	ref_transcripts_convert $reftranscripts tsvreftranscripts gtfreftranscripts
}

proc validate_sc_filter {cmd refseq reftranscripts organelles} {
# putsvars validate_sc_filter cmd refseq reftranscripts organelles
	putslog "validating sc_filter $cmd $refseq $reftranscripts $organelles"
	if {![file exists $refseq]} {
		error "reference sequence does not exist ($refseq)"
	}
	if {[auto_load validate_sc_filter_${cmd}]} {
		validate_sc_filter_${cmd} $refseq $reftranscripts $organelles
		return
	}
	if {[catch {exec which dirR}]} {
		error "command \"dirR\" not available (which i needed for \"$cmd\" sc_filter), make sure it is installed, e.g. using \"cg install dirR\""
	}
	if {$reftranscripts eq ""}  {
		ref_gtftranscripts $refseq
	} elseif {$reftranscripts in {none -}} {
		# ignore
	} elseif {![file exists $reftranscripts]} {
		error "reftranscripts file does not exist ($reftranscripts)"
	}
	ref_transcripts_convert $reftranscripts tsvreftranscripts gtfreftranscripts
}

proc validate_sc_celltyper {cmd cellmarkerfile tissue} {
# putsvars validate_sc_celltyper cmd cellmarkerfile tissue
	putslog "validating sc_celltyper $cmd $cellmarkerfile $tissue"
	if {[auto_load validate_sc_celltyper_${cmd}]} {
		validate_sc_celltyper_${cmd} $cellmarkerfile $tissue
		return
	}
	if {[catch {exec which dirR}]} {
		error "command \"dirR\" not available (which should have package \"$cmd\"), make sure it is installed, e.g. using \"cg install dirR\""
	}
	if {$cellmarkerfile eq ""}  {
		ref_gtftranscripts $refseq
	} elseif {$cellmarkerfile in {none -}} {
		# ignore
	} elseif {![file exists $cellmarkerfile]} {
		error "cellmarkerfile does not exist ($cellmarkerfile)"
	}
	set header [cg select -header $cellmarkerfile]
	if {"celltype" ni $header} {
		error "cellmarkerfile has no column named \"celltype\" ($cellmarkerfile)"
	}
	if {$tissue ne "" && "tissue" ni $header} {
		error "cellmarkerfile has no column named \"tissue\" while tissue is given as a parameter ($cellmarkerfile)"
	}
}
