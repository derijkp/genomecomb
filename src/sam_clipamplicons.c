/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "debug.h"

inline int min ( int a, int b ) { return a < b ? a : b; }

inline int max ( int a, int b ) { return a > b ? a : b; }

#define AMPLICONSBUFFERSIZE 4

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract suplementary alignment (part of chimeric alignment) */
#define BAM_FSUPPL      2048

typedef struct Amplicon {
	DString *chr2;
	int start2;
	int end2;
	int outerstart2;
	int outerend2;
} Amplicon;

typedef struct Cigar {
	int size;
	int memsize;
	int *num;
	char *action;
} Cigar;

int parse_cigar(Cigar *cigar,char *string) {
	if (cigar->memsize == 0) {
		cigar->size = 0;
		cigar->memsize = 5;
		cigar->num = malloc(cigar->memsize*sizeof(int));
		cigar->action = malloc(cigar->memsize*sizeof(char));
	}
	cigar->size = 0;
	while (*string && *string != '\t') {
		if (cigar->size >= cigar->memsize) {
			cigar->memsize += 5;
			cigar->num = realloc(cigar->num,cigar->memsize*sizeof(int));
			cigar->action = realloc(cigar->action,cigar->memsize*sizeof(char));
		}
		cigar->num[cigar->size] = (int)strtol(string,&string,10);
		cigar->action[cigar->size++] = *string++;
	}
	return(cigar->size);
}

int cigar_refend(Cigar *cigar,int begin) {
	int count = cigar->size;
	int *num = cigar->num;
	char *action = cigar->action;
	while(count--) {
		if (*action == 'M' || *action == 'D' || *action == 'N' || *action == '=' || *action == 'X') {
			begin += *num;
		}
		action++;
		num++;
	}
	return begin;
}

int cigar_ref2seq(Cigar *cigar,int begin, int pos) {
	int count = cigar->size;
	int *num = cigar->num, cur = 0, prev = 0, prevbegin = begin;
	char *action = cigar->action;
	if (begin > pos) {
		return 0;
	}
	while(count--) {
		if (*action == 'M' || *action == '=' || *action == 'X') {
			begin += *num;
			cur += *num;
		} else if (*action == 'D' || *action == 'N') {
			begin += *num;
		} else if (*action == 'I' || *action == 'S') {
			cur += *num;
		}
		if (begin > pos) break;
		prevbegin = begin;
		prev = cur;
		action++;
		num++;
	}
	if (action > cigar->action+cigar->size-1 || *action == 'M' || *action == '=' || *action == 'X') {
		return prev+(pos-prevbegin);
	} else {
		return prev;
	}
}

int amp_clip_read(DString *seq,DString *qual,int from, int to) {
	char *stringseq = seq->string+from;
	char *stringqual = qual->string+from;
	int count;
	if (from >= seq->size) {return 0;}
	if (to > seq->size) {to = seq->size;}
	count = (seq->size<to) ? seq->size : to - from;
	while (count--) {
		*stringseq++ = 'N';
		*stringqual++ = '!';
	}
	return 1;
}

int main(int argc, char *argv[]) {
	FILE *f1=stdin,*f2;
	Amplicon amplicon[AMPLICONSBUFFERSIZE];
	int close,maxclose,closepos;
	Cigar cigar;
	DStringArray *result1=NULL,*result2=NULL,*resultkeep=NULL,*resulttemp=NULL;
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*empty=NULL,*qual=NULL,*seq=NULL;
	DString *chromosome1 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	ssize_t read;
	unsigned int curpos=0, flag;
	int comp,chr2pos,start2pos,end2pos,outerstart2pos,outerend2pos,max2,i;
	unsigned int numfields2,numfields,pos1,pos2;
	int count;
	int start1,end1,amppos=0,ampnum=0,ampcur=0,found;
	int prevstart1 = -1, prevend1 = -1, prevstart2 = -1,prevend2 = -1;
	int error2,reverse;
	if ((argc != 7)) {
		fprintf(stderr,"Format is: sam_clipamplicons ampliconsfile chrpos startpos endpos outerstartpos outerendpos");
		exit(EXIT_FAILURE);
	}
	for (i = 0 ; i < AMPLICONSBUFFERSIZE ; i++) {
		amplicon[i].chr2 = DStringNew();
	}
	cigar.size = 0;
	cigar.memsize = 0;
	f2 = fopen64_or_die(argv[1],"r");
	chr2pos = atoi(argv[2]);
	start2pos = atoi(argv[3]);
	end2pos = atoi(argv[4]);
	outerstart2pos = atoi(argv[5]);
	outerend2pos = atoi(argv[6]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (outerstart2pos > max2) {max2 = outerstart2pos;} ; if (outerend2pos > max2) {max2 = outerend2pos;} ;
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	curchromosome = DStringEmtpy();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(10+2);
	result2 = DStringArrayNew(max2+2);
	resultkeep = DStringArrayNew(max2+2);
	read = DStringGetLine(line1, f1);
	while (read != -1) {
		if (line1->string[0] != '@') break;
		fprintf(stdout,"%s\n",line1->string);
		read = DStringGetLine(line1, f1); curpos++;
		pos1++;
	}
	skip_header(f2,line2,&numfields2,&pos2);
	error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields);	pos2++;
	if (!error2) {
		check_numfieldserror(numfields,numfields2,line2,argv[1],&pos2);
		DStringCopy(amplicon[0].chr2,result2->data+chr2pos);
		sscanf(result2->data[start2pos].string,"%d",&(amplicon[0].start2));
		sscanf(result2->data[end2pos].string,"%d",&(amplicon[0].end2));
		sscanf(result2->data[outerstart2pos].string,"%d",&(amplicon[0].outerstart2));
		sscanf(result2->data[outerend2pos].string,"%d",&(amplicon[0].outerend2));
		amplicon[0].start2++; amplicon[0].end2++; amplicon[0].outerstart2++; amplicon[0].outerend2++;
		chromosomekeep = amplicon[ampcur].chr2;
		ampnum = 1;
	}
	DStringSplitTab(line1,10,result1,0,&numfields);
	while (1) {
		pos1++;
		/* check_numfieldserror(numfields,15,line1,"stdin",&pos1); */
		sscanf(result1->data[1].string,"%d",&(flag));
		chromosome1 = result1->data+2;
		seq = result1->data+9;
		qual = result1->data+10;
		if (chromosome1->string[0] == '*' || seq->string[0] == '*') {
			/* skip unmapped reads */
			NODPRINT("unmapped")
			fprintf(stdout,"%s\n",line1->string); 
			if (DStringGetTab(line1,f1,10,result1,0,&numfields)) break;
			continue;
		}
		sscanf(result1->data[3].string,"%d",&start1);
		checksortreg(curchromosome,&prevstart1,&prevend1,chromosome1,start1,start1,"input sam file");
		parse_cigar(&cigar,result1->data[5].string);
		end1 = cigar_refend(&cigar,start1);
		reverse = flag & BAM_FREVERSE;
		ampcur = amppos+ampnum-1; if (ampcur >= AMPLICONSBUFFERSIZE) {ampcur -= AMPLICONSBUFFERSIZE;}
		while (!error2) {
			comp = DStringLocCompare(amplicon[ampcur].chr2, chromosome1);
			if (comp > 0) break;
			if (comp == 0) {
				if (reverse) {
					if (end1 < amplicon[ampcur].outerend2) break;
				} else {
					if (start1 < amplicon[ampcur].start2) break;
				}
			}
			/* add new ones to end */
			ampcur++; if (ampcur >= AMPLICONSBUFFERSIZE) {ampcur -= AMPLICONSBUFFERSIZE;}
			/* add next amplicon to stack */
			/* keep data of previous */
			/* to avoid allocating new memory everytime, reuse linekeep and associated data */
			linetemp = linekeep;
			linekeep = line2;
			line2 = linetemp;
			resulttemp = resultkeep;
			resultkeep = result2;
			result2 = resulttemp;
			/* get new line */
			error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields); pos2++;
			if (error2)  {
				comp = -1;
				break;
			} else {
				check_numfieldserror(numfields,numfields2,line2,argv[1],&pos2);
			}
			DStringCopy(amplicon[ampcur].chr2,result2->data+chr2pos);
			sscanf(result2->data[start2pos].string,"%d",&(amplicon[ampcur].start2));
			sscanf(result2->data[end2pos].string,"%d",&(amplicon[ampcur].end2));
			sscanf(result2->data[outerstart2pos].string,"%d",&(amplicon[ampcur].outerstart2));
			sscanf(result2->data[outerend2pos].string,"%d",&(amplicon[ampcur].outerend2));
			amplicon[ampcur].start2++; amplicon[ampcur].end2++; amplicon[ampcur].outerstart2++; amplicon[ampcur].outerend2++;
			if (ampnum < AMPLICONSBUFFERSIZE) {ampnum++;} else {
				amppos++; if (amppos >= AMPLICONSBUFFERSIZE) {amppos = 0;}
			}
			
			comp = DStringLocCompare(amplicon[ampcur].chr2, chromosomekeep);
			if (comp < 0 || (comp == 0 && (amplicon[ampcur].start2 < prevstart2 || (amplicon[ampcur].start2 == prevstart2 && amplicon[ampcur].end2 < prevend2)))) {
				fprintf(stderr,"Cannot annotate because the amplicon file is not correctly sorted (sort correctly using \"cg select -s -\"):\n");
				fprintf(stderr,"%*.*s:%d-%d came before %*.*s:%d-%d\n",
					chromosomekeep->size,chromosomekeep->size,chromosomekeep->string,
					prevstart2,prevend2,
					amplicon[ampcur].chr2->size,amplicon[ampcur].chr2->size,amplicon[ampcur].chr2->string,
					amplicon[ampcur].start2,amplicon[ampcur].end2
				);
				exit(1);
			}
			prevstart2 = amplicon[ampcur].start2; prevend2 = amplicon[ampcur].end2;
			chromosomekeep = amplicon[ampcur].chr2;
		}
		while (ampnum) {
			comp = DStringLocCompare(amplicon[amppos].chr2, chromosome1);
			if (comp > 0) break;
			if (comp == 0 && amplicon[amppos].outerend2 >= start1) break;
			amppos++; if (amppos >= AMPLICONSBUFFERSIZE) {amppos = 0;}
			ampnum--;
		}
		ampcur = amppos;
		count = ampnum;
		if (reverse) {
			found = 0;
			close = -1; maxclose = INT_MAX; closepos = -1;
			while (count--) {
				comp = DStringLocCompare(amplicon[ampcur].chr2, chromosome1);
				if (comp == 0) {
					if (end1 > amplicon[ampcur].end2 && end1 <= amplicon[ampcur].outerend2) {
						found=1;
						break;
					}
					if (start1 >= amplicon[ampcur].outerstart2 && end1 <= amplicon[ampcur].start2) {
						amp_clip_read(seq,qual,0,cigar_ref2seq(&cigar,start1,amplicon[ampcur].start2));
						found = 0; closepos = -1;
						break;
					}
					close = abs(amplicon[ampcur].outerend2 - end1);
					if (close < maxclose) {maxclose = close; closepos = ampcur;}
				}
				ampcur++; if (ampcur >= AMPLICONSBUFFERSIZE) {ampcur -= AMPLICONSBUFFERSIZE;}
			}
			if (!found && closepos != -1) {
				found=1; ampcur = closepos;
			}
			if (found) {
				if (amplicon[ampcur].start2 > start1 && amplicon[ampcur].start2 < end1) {
					amp_clip_read(seq,qual,0,cigar_ref2seq(&cigar,start1,amplicon[ampcur].start2));
				}
				if (amplicon[ampcur].outerend2 > start1 && amplicon[ampcur].end2 < end1) {
					amp_clip_read(seq,qual,cigar_ref2seq(&cigar,start1,amplicon[ampcur].end2),seq->size);
				}
			}
		} else {
			found = 0;
			close = -1; maxclose = INT_MAX; closepos = -1;
			while (count--) {
				comp = DStringLocCompare(amplicon[ampcur].chr2, chromosome1);
				if (comp == 0) {
					if (start1 < amplicon[ampcur].start2 && start1 >= amplicon[ampcur].outerstart2) {
						found = 1;
						break;
					}
					if (start1 >= amplicon[ampcur].end2 && end1 <= amplicon[ampcur].outerend2) {
						amp_clip_read(seq,qual,0,cigar_ref2seq(&cigar,start1,amplicon[ampcur].outerend2));
						found = 0; closepos = -1;
						break;
					}
					close = abs(amplicon[ampcur].outerstart2-start1);
					if (close < maxclose) {maxclose = close; closepos = ampcur;}
				}
				ampcur++; if (ampcur >= AMPLICONSBUFFERSIZE) {ampcur -= AMPLICONSBUFFERSIZE;}
			}
			if (!found && closepos != -1) {
				found=1; ampcur = closepos;
			}
			if (found) {
				if (amplicon[ampcur].start2 > start1 && amplicon[ampcur].outerstart2 < end1) {
					amp_clip_read(seq,qual,0,cigar_ref2seq(&cigar,start1,amplicon[ampcur].start2));
				}
				if (amplicon[ampcur].end2 > start1 && end1 >= amplicon[ampcur].end2) {
					amp_clip_read(seq,qual,cigar_ref2seq(&cigar,start1,amplicon[ampcur].end2),seq->size);
				}
			}
		}
		fprintf(stdout,"%s\n",line1->string);
		if (DStringGetTab(line1,f1,10,result1,0,&numfields)) break;
	}
	fclose(f1);
	fclose(f2);
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (empty) {DStringDestroy(empty);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	if (resultkeep) {DStringArrayDestroy(resultkeep);}
	exit(EXIT_SUCCESS);
}
