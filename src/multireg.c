/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "tools.h"
#include "tools_var.h"
#include "gztools.h"
#include "debug.h"

typedef struct RegFile {
	char *file;
	GZFILE *f;
	DString *line;
	DStringArray *result;
	VariantPos varpos;
	int max;
	unsigned int numfields;
	unsigned int pos;
	unsigned int error;
	int isreg;
	DString *chr;
	int start;
	int end;
} RegFile;

RegFile *OpenRegfile(char *filename,char *isreg) {
	unsigned int numfields;
	RegFile *regfile = malloc(sizeof(RegFile));
	DStringArray *result;
	if (isreg[0] == '1') {regfile->isreg = 1;} else {regfile->isreg = 0;}
	regfile->file = filename;
	regfile->f = gz_open(filename);
	regfile->line = DStringNew();
	gz_skip_header(regfile->f,regfile->line,&regfile->numfields,&regfile->pos);
	regfile->result = DStringArrayNew(regfile->numfields+2);
	DStringSplitTab(regfile->line, regfile->numfields, regfile->result, 0,NULL);
	varpos_fromheader(&regfile->varpos,regfile->result);
	if (regfile->varpos.chr == -1 || regfile->varpos.start == -1 || regfile->varpos.end == -1) {
		fprintf(stderr,"header error: fields (or alternatives) not found: chromosome begin end\n");
		exit(1);
	}
	if (regfile->isreg) {
		regfile->max = regfile->numfields;
	} else {
		regfile->max = regfile->varpos.chr;
		if (regfile->max < regfile->varpos.start) {regfile->max = regfile->varpos.start;}
		if (regfile->max < regfile->varpos.end) {regfile->max = regfile->varpos.end;}
	}
	result = regfile->result;
	if (regfile->isreg) {
		register int i;
		for (i = 0; i < regfile->numfields ; i++) {
			if (i == regfile->varpos.chr || i == regfile->varpos.start || i == regfile->varpos.end) continue;
			putc_unlocked('\t',stdout);
			DStringputs(result->data+i,stdout);
		}
	} else {
		char *c = filename, *start = filename,*end = filename,*pend = filename;
		putc_unlocked('\t',stdout);
		while (*c != '\0') {
			if (*c == '/') {start = c+1; end = start; pend = start;}
			if (*c == '.') {pend = end; end = c;}
			c++;
		}
		if (regfile->f->type != UNCOMPRESSED && pend != start) {
			end = pend;
		}
		if (end == start) {end = c;}
		while (*start != '\0') {
			if (start == end) break;
			putc_unlocked(*start++,stdout);
		}
	}
	regfile->error = gz_DStringGetTab(regfile->line,regfile->f,regfile->max,regfile->result,1,&numfields);
	if (!regfile->error) {
		regfile->chr = regfile->result->data+regfile->varpos.chr;
		sscanf(regfile->result->data[regfile->varpos.start].string,"%d",&(regfile->start));
		sscanf(regfile->result->data[regfile->varpos.end].string,"%d",&(regfile->end));
	}
	return regfile;
}

#define OPEN 1
#define OPEN 1

int main(int argc, char *argv[]) {
	RegFile **regfiles = NULL, *regfile = NULL;
	int numtodo, numfiles;
	DString *curchr=DStringNew(),*nextchr=DStringNew();
	DString *maxchr = DStringNew();
	int curpos=-1,nextpos=-1,anymatch=0, comp;
	unsigned int numfields;
	register unsigned int i,pos;

	DStringSetS(maxchr,"Z",1);
/*	maxchr->string[0] = 254; */
	if ((argc < 2)) {
		fprintf(stderr,"Format is: multireg file1 isreg file2 isreg ...\n");
		exit(EXIT_FAILURE);
	}
	numtodo = argc -1;
	numfiles = numtodo/2;
	regfiles = (RegFile **)malloc(numtodo*sizeof(RegFile*));
	pos = 1;
	fprintf(stdout,"chromosome\tbegin\tend");
	for (i = 0 ; i < numfiles ; i++) {
		regfiles[i] = OpenRegfile(argv[pos],argv[pos+1]);
		pos+=2;
	}
	putc_unlocked('\n',stdout);
	while(1) {
		nextchr = maxchr; nextpos = 2147483647;
		anymatch = 0;
		for (i = 0 ; i < numfiles ; i++) {
			regfile = regfiles[i];
			if (regfile->error) continue;
			/* move to next region if needed */
			comp = DStringLocCompare(regfile->chr,curchr);
			if (comp < 0 || (comp == 0 && regfile->end <= curpos)) {
				regfile->error = gz_DStringGetTab(regfile->line,regfile->f,regfile->max,regfile->result,1,&numfields);
				if (!regfile->error) {
					regfile->chr = regfile->result->data+regfile->varpos.chr;
					sscanf(regfile->result->data[regfile->varpos.start].string,"%d",&(regfile->start));
					sscanf(regfile->result->data[regfile->varpos.end].string,"%d",&(regfile->end));
				}
				comp = DStringLocCompare(regfile->chr,curchr);
				if (comp < 0) {
					fprintf(stderr,"File (%s) is not correctly sorted (sort correctly using \"cg select -s -\")\n",regfile->file);
					exit(1);
				}
			}
			if (comp == 0) {
				/* on same chromosome, put same object so easier comparison later */
				regfile->chr = curchr;
				if (nextchr != curchr) {
					nextchr = curchr;
					nextpos = 2147483647;
				}
				if (regfile->start > curpos) {
					if (regfile->start < nextpos) {nextpos = regfile->start;}
				} else if (regfile->end < curpos) {
					fprintf(stderr,"File (%s) is not correctly sorted (sort correctly using \"cg select -s -\")\n",regfile->file);
					exit(1);
				} else if (regfile->end > curpos) {
					/* no need to check end: curpos should not be anle to past end*/
					anymatch = 1;
					if (regfile->end < nextpos) {nextpos = regfile->end;}
				}
			} else if (nextchr != curchr) {
				/* if nextchr is already curchr, this loc cannot be closer anyway */
				comp = DStringLocCompare(regfile->chr,nextchr);
				if (comp == 0) {
					if (regfile->start < nextpos) {nextpos = regfile->start;}
				} else if (comp < 0) {
					nextchr = regfile->chr;
					nextpos = regfile->start;
				}
			}
		}
		if (anymatch) {
			/* write only region starting from curpos if any match was found */
			fprintf(stdout,"%s\t%d\t%d",curchr->string,curpos,nextpos);
			for (i = 0 ; i < numfiles ; i++) {
				regfile = regfiles[i];
				if (regfile->error || regfile->chr != curchr || regfile->start >= nextpos) {
					if (regfile->isreg) {
						register int i;
						i = regfile->numfields;
						for (i = 0 ; i < regfile->numfields ; i++) {
							if (i == regfile->varpos.chr || i == regfile->varpos.start || i == regfile->varpos.end) continue;
							putc_unlocked('\t',stdout);
							putc_unlocked('0',stdout);
						}
					} else {
						putc_unlocked('\t',stdout);
						putc_unlocked('0',stdout);
					}
				} else {
					if (regfile->isreg) {
						register int i;
						for (i = 0 ; i < regfile->numfields ; i++) {
							if (i == regfile->varpos.chr || i == regfile->varpos.start || i == regfile->varpos.end) continue;
							putc_unlocked('\t',stdout);
							DStringputs(regfile->result->data+i,stdout);
						}
					} else {
						putc_unlocked('\t',stdout);
						putc_unlocked('1',stdout);
					}
				}
			}
			putc_unlocked('\n',stdout);
		}
		if (nextchr == maxchr) break;
		if (nextchr != curchr) {
			DStringCopy(curchr,nextchr);
		}
		curpos = nextpos;
	}
	for (i = 0 ; i < numfiles ; i++) {
		regfile = regfiles[i];
		gz_close(regfile->f);
	}
	exit(EXIT_SUCCESS);
}
