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
#include "tools.h"
#include "debug.h"
#include "khash-dstring.h"
#include "tools_bcol.h"

typedef struct Dest {
	FILE *f;
	FILE *rf;
	uint64_t start;
	uint64_t lastpos;
} Dest;

KHASH_MAP_INIT_DSTR(DSTRING, Dest*);

/*#include "dstring-khash.h"*/

Dest *bcol_make_getout(khash_t(DSTRING) *hashtable,char *pre,DString *chromosome) {
	DString *buffer = NULL;
	Dest *o;
	khiter_t k;
	int ret;
	k = kh_put(DSTRING,hashtable, chromosome, &ret);
	if (ret == 0) {
		/* key was already present in the hashtable */
		o = kh_value(hashtable, k);
	} else {
		o = (Dest *)malloc(sizeof(Dest));
		o->start = 0;
		o->lastpos = -1;
		buffer = DStringNew();
		DStringAppend(buffer,pre);
		if (chromosome->size > 0) {
			DStringAppendS(buffer,chromosome->string,chromosome->size);
		}
		DStringAppend(buffer,".bcol");
		o->rf = fopen64(buffer->string,"w");
		if (o->rf == NULL) {
			fprintf(stderr,"Error opening file %s: %s.\n", buffer->string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		DStringAppend(buffer,".bin");
		o->f = fopen64(buffer->string,"w");
		if (o->f == NULL) {
			fprintf(stderr,"Error opening file %s: %s.\n", buffer->string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		kh_value(hashtable, k) = o;
		DStringSetS(buffer,chromosome->string,chromosome->size);
		kh_key(hashtable, k) = buffer;
		buffer = NULL;
	}
	return o;
}

int main(int argc, char *argv[]) {
	khash_t(DSTRING) *hashtable;
	Dest *o,*po = NULL;
	DString *result = NULL;
	DString *line = NULL,*chromosome = NULL;
	khiter_t k;
	char *pre,*type,*defaultvalue;
	uint64_t offset, poffset = -1, size;
	int reverse = 0, isunsigned = 0;
	int col = 0,max = 0,offsetcol = -1,chrcol = -1;
	if ((argc < 2)||(argc > 8)) {
		fprintf(stderr,"Format is: bcol_make output_pre type ?col? ?chromosomecol? ?offsetcol? ?default?\n");
		exit(EXIT_FAILURE);
	}
	pre = argv[1];
	if (argc >= 3) {
		type = argv[2];
	}
	if (argc >= 4) {
		col = atoi(argv[3]);
		max = col;
	}
	if (argc >= 5) {
		chrcol = atoi(argv[4]);
		if (chrcol > max) {max = chrcol;}
	}
	if (argc >= 6) {
		offsetcol = atoi(argv[5]);
		if (offsetcol > max) {max = offsetcol;}
	}
	if (argc >= 7) {
		defaultvalue = argv[6];
	}
NODPRINT("bcol_make %s %s %d %d %d\n",pre,type,col,chrcol,offsetcol)
	line = DStringNew();
	reverse = bcol_NeedReversing((int)type[0]);
	if (type[1] == 'u') {isunsigned = 1;}
	result = DStringArrayNew(max+1);
	hashtable = kh_init(DSTRING);
	if (chrcol == -1) {
		o = bcol_make_getout(hashtable,pre,DStringEmtpy());
		poffset = -1;
	}
	while (!DStringGetTab(line,stdin,max,result,0)) {
		if (chrcol != -1) {
			chromosome = result+chrcol;
			o = bcol_make_getout(hashtable,pre,chromosome);
			if (o != po) {
				po = o;
				poffset = -1;
			}
		}
		if (offsetcol != -1) {
			offset = atoll(result[offsetcol].string);
			if (poffset == -1) {
				o->start = offset;
			} else if (poffset != offset) {
				size = offset - poffset;
				if (size < 0) {
					fprintf(stderr, "error: cannot make position based bcol on unsorted data ($offset < $poffset sort on position first)\n");
					exit(EXIT_FAILURE);
				}
				while (poffset < offset) {
					bcol_printbin(o->f,reverse,isunsigned,type,defaultvalue);
					poffset++;
				}
				o->lastpos = o->lastpos + size;
			}
			poffset = offset+1;
		}
		NODPRINT("s=%s\n",result[col].string)
		bcol_printbin(o->f,reverse,isunsigned,type,result[col].string);
		o->lastpos ++;
	}
	for (k = kh_begin(hashtable); k != kh_end(hashtable); ++k) {
		if (kh_exist(hashtable, k)) {
			DStringDestroy(kh_key(hashtable, k));
			o = kh_value(hashtable, k);
			fclose(o->f);
			fprintf(o->rf,"# binary column\n");
			fprintf(o->rf,"# type %s\n",type);
			fprintf(o->rf,"# default %s\n","0");
			fprintf(o->rf,"begin\ttype\toffset\n");
			fprintf(o->rf,"%llu\t%s\t%d\n",o->start,type,0);
			fprintf(o->rf,"%llu\tend\t%d\n",o->start + o->lastpos,0);
			fclose(o->rf);
			free(o);
		}
	}
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	kh_destroy(DSTRING,hashtable);
	exit(EXIT_SUCCESS);
}

