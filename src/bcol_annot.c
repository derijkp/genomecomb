/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

/* #define DEBUG 1 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "tools_bcol.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	BCol *bcol;
	BCol_table *table;
	int tablepos = 0,tablesize;
	FILE *f1;
	DStringArray *mvalues = NULL;
	DStringArray *result1=NULL;
	DString *defaultvalue;
	DString *line1 = NULL;
	DString *prevchromosome1 = DStringNew(), *prevchromosome2 = DStringNew();
	DString *prevtype1 = DStringNew();
	DString *prevalt1 = DStringNew();
	DString *chromosome1=NULL,*chromosome2=NULL,*type1 = NULL,*alt1 = NULL;
	unsigned int numfields1,numfields,pos1 = 0;
	off_t binpos,curstartpos = 0;
	int prevstart1 = -1,prevend1 = -1;
	int chr1pos,start1pos,end1pos,type1pos,alt1pos,max1,precision;
	int comp;
	int start1,end1;
	int start2,end2;
	int nextpos=0;
	if ((argc != 9)) {
		fprintf(stderr,"Format is: bcol_annot file1 chrpos1 startpos1 endpos1 type1pos alt1pos bcolfile precision");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	type1pos = atoi(argv[5]);
	alt1pos = atoi(argv[6]);
	precision = atoi(argv[8]);
	if (type1pos == -1) {
		type1 = DStringNew();
	}
	if (alt1pos == -1) {
		alt1 = DStringNew();
	}
	max1 = chr1pos;
	if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (type1pos > max1) {max1 = type1pos;} ; if (alt1pos > max1) {max1 = alt1pos;} ;
	line1 = DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(max1+2);
	bcol = bcol_open(argv[7]);
	if (precision == -1) {precision = bcol->precision;}
	if (alt1pos != -1) {
		if (bcol->multi->size == 0) {
			fprintf(stderr,"bcol_annot var analysis needs a multivalue bcol file");
			exit(EXIT_FAILURE);
		}
		mvalues = DStringArrayFromChar(bcol->multi->string,',');
	}
	table = bcol->table;
	tablesize = bcol->tablesize;
	defaultvalue = bcol->def;
	NODPRINT("bcol_annot %s %d %d %d %d %d %s %d %d %d %s ...",
		argv[1],chr1pos,start1pos,end1pos,type1pos,alt1pos,
		argv[7],chr2pos,start2pos,end2pos,binfile
	);
	skip_header(f1,line1,&numfields1,&pos1);
	tablepos = 0;
	chromosome2 = table[tablepos].chr;
	start2 = table[tablepos].begin;
	end2 = table[tablepos].end;
	curstartpos = table[tablepos].pos;
	NODPRINT("line2 %s,%d,%d",Loc_ChrString(chromosome2),start2,end2)
	while (!DStringGetTab(line1,f1,max1,result1,1,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,argv[1],&pos1);
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
		if (type1pos != -1) {
			type1 = result1->data+type1pos;
		}
		if (alt1pos != -1) {
			alt1 = result1->data+alt1pos;
		}
		NODPRINT("line1 %s,%d,%d (type=%3.3s) %s",chromosome1->string,start1,end1,type1->string,alt1->string)
		NODPRINT("line2 %s,%d,%d",chromosome2->string,start2,end2)
		checksort(prevchromosome1,&prevstart1,&prevend1,prevtype1,prevalt1,chromosome1,start1,end1,type1,alt1,argv[1],&nextpos,1);
		comp = DStringLocCompare(chromosome2, chromosome1);
		while (tablepos < tablesize && ((comp < 0) || ((comp == 0) && ((end2 < start1) || (end2 == start1 && start1 != end1))))) {
			tablepos += 1;
			chromosome2 = table[tablepos].chr;
			start2 = table[tablepos].begin;
			end2 = table[tablepos].end;
			curstartpos = table[tablepos].pos;
			comp = DStringLocCompare(chromosome2,chromosome1);
		}
		if (tablepos >= tablesize || (comp > 0) || (alt1pos != -1 && end1 - start1 != 1) || ((comp == 0) && ((end1 < start2) || (end1 == start2 && start1 != end1)))) {
			NODPRINT("no overlap or ref != 1")
			if (alt1pos == -1) {
				fprintf(stdout,"%s\n",defaultvalue->string);
			} else {
				char *str = alt1->string;
				int cursize = 0,c;
				while (1) {
					c = str[cursize];
					if (c == ',' || c =='\0' || c =='\t') {
						if (str != alt1->string) {
							fprintf(stdout,",");
						}
						fprintf(stdout,"%s",defaultvalue->string);
						if (c =='\0' || c =='\t') break;
						str += cursize + 1;
						cursize = 0;
					} else {
						cursize++;
					}
				}
				fprintf(stdout,"\n");
			}
		} else { /* start1 = start position of query, start2 = start position of db block */
			if (alt1pos == -1) {
				binpos = curstartpos + start1 - start2;
				bcol_getbin(bcol,binpos,binpos);
				bcol_printtext(stdout,bcol->reverse,bcol->isunsigned,bcol->type,bcol->buffer,precision);
			} else {
				char *str = alt1->string;
				int cursize=0,found,c,i;
				binpos = mvalues->size * (curstartpos + start1 - start2);
				while (1) {
					c = str[cursize];
					if (c == ',' || c =='\0' || c =='\t') {
						found = -1;
						for (i=0 ; i < mvalues->size ; i++) {
							DString *m = mvalues->data+i;
							if (cursize == m->size && strncmp(m->string,str,cursize) == 0) {
								found = i; break;
							}
						}
						if (str != alt1->string) {
							fprintf(stdout,",");
						}
						if (found == -1) {
							fprintf(stdout,"%s",defaultvalue->string);
						} else {
							bcol_getbin(bcol,binpos + found,binpos + found);
							bcol_printtext(stdout,bcol->reverse,bcol->isunsigned,bcol->type,bcol->buffer,precision);
						}
						if (c =='\0' || c =='\t') break;
						str += cursize + 1;
						cursize = 0;
					} else {
						cursize++;
					}
				}
			}
			fprintf(stdout,"\n");
		}
	}
	fclose(f1);
	bcol_close(bcol);
	DStringDestroy(prevchromosome1); DStringDestroy(prevtype1); DStringDestroy(prevalt1);
	DStringDestroy(prevchromosome2);
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
