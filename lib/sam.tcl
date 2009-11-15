#1	QNAME	Query (pair) NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (\u2018=\u2019 if same as RNAME)
#8	MPOS	1-based Mate POSistion
#9	ISIZE	Inferred insert SIZE
#10	SEQ	query SEQuence on the same strand as the reference
#11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
#12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

#0x0001	the read is paired in sequencing
#0x0002	the read is mapped in a proper pair
#0x0004	the query sequence itself is unmapped
#0x0008	the mate is unmapped
#0x0010	strand of the query (1 for reverse)
#0x0020	strand of the mate
#0x0040	the read is the first read in a pair
#0x0080	the read is the second read in a pair
#0x0100	the alignment is not primary
#0x0200	the read fails platform/vendor quality checks
#0x0400	the read is either a PCR or an optical duplicate

foreach {mask var} {
	0x001	paired
	0x002	mapped_pair
	0x004	unmapped
	0x008	unmapped_mate
	0x010	reverse
	0x020	reverse_mate
	0x040	pair_first
	0x080	pair_second
	0x100	notprimary
	0x200	failed
	0x400	duplicate
} {
	if {[expr {$flags & $mask}]} {set a($var) 1} else {set a($var) 0}
}

foreach {QNAME FLAG RNAME POS MAPQ CIGAR MRNM MPOS ISIZE SEQ QUAL OPT} $line break


