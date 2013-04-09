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
#include "hash.h"
#include "tools_bcol.h"

typedef struct Dest {
	FILE *f;
	FILE *rf;
	uint64_t start;
	uint64_t lastpos;
} Dest;

Dest *bcol_make_getout(Hash_table *hashtable,char *pre,DString *chromosome) {
	DString *buffer = NULL;
	Dest *o;
	Hash_bucket *bucket;
	int new;
	bucket = hash_get(hashtable, (void *)chromosome, hash_Dstring_hash, hash_Dstring_compare, &new);
	if (new == 0) {
		/* key was already present in the hashtable */
		o = hash_getvalue(bucket);
	} else {
		o = (Dest *)malloc(sizeof(Dest));
		NODPRINT("new %s: %p",chromosome->string,o)
		o->start = 0;
		o->lastpos = -1;
		buffer = DStringNew();
		DStringAppend(buffer,pre);
		if (chromosome->size > 0) {
			DStringAppendS(buffer,chromosome->string,chromosome->size);
		}
		DStringAppend(buffer,".bcol");
		o->rf = fopen64_or_die(buffer->string,"w");
		DStringAppend(buffer,".bin");
		o->f = fopen64_or_die(buffer->string,"w");
		hash_setvalue(bucket,o)
		DStringSetS(buffer,chromosome->string,chromosome->size);
		bucket->key = (void *)buffer;
		hash_setkey(bucket,buffer)
		buffer = NULL;
	}
	return o;
}

int main(int argc, char *argv[]) {
	Hash_table *hashtable;
	Dest *o,*po = NULL;
	DStringArray *result = NULL;
	DString *line = NULL,*chromosome = NULL;
	Hash_iter iter;
	Hash_bucket *bucket;
	char *pre,*type,*defaultvalue = "";
	uint64_t offset, poffset = -1, size;
	int reverse = 0, isunsigned = 0, header = 0;
	int col = 0,max = 0,offsetcol = -1,chrcol = -1;
	if ((argc < 2)||(argc > 9)) {
		fprintf(stderr,"Format is: bcol_make output_pre type ?col? ?chromosomecol? ?offsetcol? ?default? ?header?\n");
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
	line = DStringNew();
	if (argc >= 8) {
		header = atoi(argv[7]);
		if (header) {
			skip_header(stdin,line,NULL);
		}
	}
	NODPRINT("bcol_make %s %s %d %d %d\n",pre,type,col,chrcol,offsetcol)
	reverse = bcol_NeedReversing((int)type[0]);
	if (type[1] == 'u') {isunsigned = 1;}
	result = DStringArrayNew(max+2);
	hashtable = hash_init();
	if (chrcol == -1) {
		o = bcol_make_getout(hashtable,pre,DStringEmtpy());
		poffset = -1;
	}
	while (!DStringGetTab(line,stdin,max,result,0,NULL)) {
		if (chrcol != -1) {
			chromosome = result->data+chrcol;
			o = bcol_make_getout(hashtable,pre,chromosome);
			if (o != po) {
				po = o;
				poffset = -1;
			}
		}
		if (offsetcol != -1) {
			offset = atoll(result->data[offsetcol].string);
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
		NODPRINT("s=%s\n",result->data[col].string)
		bcol_printbin(o->f,reverse,isunsigned,type,result->data[col].string);
		o->lastpos ++;
	}
	bucket = hash_first(hashtable,&iter);
	while(bucket != NULL) {
		DString *ds = hash_getkey(bucket);
		o = hash_getvalue(bucket);
		NODPRINT("close %s: %p",ds->string,o)
		fclose(o->f);
		fprintf(o->rf,"# binary column\n");
		fprintf(o->rf,"# type %s\n",type);
		fprintf(o->rf,"# default %s\n","0");
		fprintf(o->rf,"begin\ttype\toffset\n");
		fprintf(o->rf,"%llu\t%s\t%d\n",(long long int)o->start,type,0);
		fprintf(o->rf,"%llu\tend\t%d\n",(long long int)(o->start + o->lastpos),0);
		fclose(o->rf);
		DStringDestroy(ds);
		free(o);
		bucket = hash_next(&iter);
	}
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	hash_destroy(hashtable);
	exit(EXIT_SUCCESS);
}

