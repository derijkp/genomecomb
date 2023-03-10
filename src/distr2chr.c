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
#include "tools.h"
#include "debug.h"
#include "hash.h"

int main(int argc, char *argv[]) {
	Hash_table *hashtable;
	FILE *o;
	DStringArray *result = NULL;
	DString *buffer = NULL;
	DString *line = NULL,*chromosome = NULL;
	Hash_bucket *bucket;
	Hash_iter iter;
	int col = 0,new,header=0;
	if ((argc < 2)&&(argc > 4)) {
		fprintf(stderr,"Format is: distr2chr output_pre ?col? ?header?");
		exit(EXIT_FAILURE);
	}
	if (argc >= 3) {
		col = atoi(argv[2]);
	}
	line = DStringNew();
	if (argc == 4) {
		header = atoi(argv[3]);
		if (header) {
			skip_header(stdin,line,NULL,NULL);
		}
	}
	result = DStringArrayNew(col+2);
	hashtable = hash_init();
	while (!DStringGetTab(line,stdin,col,result,0,NULL)) {
		chromosome = result->data + col;
		bucket = hash_get(hashtable, (void *)chromosome, hash_DString_hash, hash_DString_compare, &new,1);
		if (new == 0) {
			/* key was already present in the hashtable */
			o = hash_getvalue(bucket);
		} else {
			buffer = DStringNew();
			DStringAppend(buffer,argv[1]);
			DStringAppendS(buffer,chromosome->string,chromosome->size);
			o = fopen64_or_die(buffer->string,"a");
			hash_setvalue(bucket,o)
			DStringSetS(buffer,chromosome->string,chromosome->size);
			hash_setkey(bucket,buffer)
			buffer = NULL;
		}
		fprintf(o,"%s\n", line->string);
	}
	bucket = hash_first(hashtable,&iter);
	while(bucket != NULL) {
		DString *ds = hash_getkey(bucket);
		DStringDestroy(ds);
		o = hash_getvalue(bucket);
		FCLOSE(o);
		bucket = hash_next(&iter);
	}
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	hash_destroy(hashtable,NULL,NULL);
	exit(EXIT_SUCCESS);
}
