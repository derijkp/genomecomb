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
#include <math.h>
#include "tools.h"
#include "debug.h"
#include "hash.h"

typedef struct Colinfo {
	char type;
	Hash_table *hashtable;
	int numnum;
	double min;
	double max;
	char *extra;
} ColInfo;

int isint(char *string, int *test) {
	char *end;
	*test = strtol(string,&end,10);
	if (end == string || *end != '\0') {
		return 0;
	} else {
		return 1;
	}
}

int isdouble(char *string, double *test) {
	char *end;
	*test = strtod(string,&end);
	if (end == string || *end != '\0') {
		return 0;
	} else {
		return 1;
	}
}

int main(int argc, char *argv[]) {
	ColInfo *colinfo;
	FILE *f1,*f2,*f3;
	DStringArray *result1=NULL;
	DString *line1 = NULL, *buffer = NULL;
	Hash_bucket *bucket;
	Hash_iter iter;
	char *cur,*prefix;
	off_t fpos;
	uint64_t count = -1, offset = 0L, progress = 20000000L;
	uint64_t data;
	double testd;
	long int temp;
	int chr1pos,start1pos,end1pos,type1pos,max1,i,isnum,new,testi,size;
	if ((argc != 10)) {
		fprintf(stderr,"Format is: bcol_indexfile file indexfile indexfilebin chrpos startpos endpos typepos colinfoprefix columns");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[4]);
	start1pos = atoi(argv[5]);
	end1pos = atoi(argv[6]);
	type1pos = atoi(argv[7]);
	prefix = argv[8];
	cur = argv[9];
	max1 = 1;
	while (*cur != '\0') {
		if (*cur == ' ' || *cur == '\n') {max1++;}
		cur++;
	}
	if (start1pos > max1) {max1 = start1pos;}
	if (end1pos > max1) {max1 = end1pos;}
	if (type1pos > max1) {max1 = type1pos;} ;
	colinfo = (ColInfo *)malloc(max1*sizeof(ColInfo));
	for(i = 0 ; i < max1 ; i++) {
		colinfo[i].type = 'i';
		colinfo[i].hashtable = hash_init_size(500);
		colinfo[i].numnum = 0;
		colinfo[i].min = INFINITY;
		colinfo[i].max = -INFINITY;
		colinfo[i].extra = 0;
	}
	f2 = fopen64_or_die(argv[2],"w");
	fprintf(f2,"# binary column\n");
	fprintf(f2,"# type wu\n");
	fprintf(f2,"begin\ttype\toffset\n");
	f3 = fopen64_or_die(argv[3],"w");
NODPRINT("poss: %d:%d-%d %d",chr1pos,start1pos,end1pos,type1pos)
	/* allocate */
	line1 = DStringNew();
	result1 = DStringArrayNew(max1+2); /* need one extra for the remainder */
	skip_header(f1,line1,NULL,NULL);
	fpos = ftello(f1);
	fprintf(f2,"%d\t%s\t%d\n",0,"wu",0);
	while (!DStringGetTab(line1,f1,max1,result1,1,NULL)) {
		count++;
		data = (uint64_t)(fpos-offset);
		fwrite(&data,8,1,f3);
/*
DPRINT("poss: %s:%s-%s",result1->data[chr1pos].string,result1->data[start1pos].string,result1->data[end1pos].string)
		{
		char *chromosome1;
		int start1,end1;
		chromosome1 = result1->data[chr1pos].string;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
		NODPRINT("%d\t%s\t%d\t%d\t%d",1,chromosome1,start1,end1)
		}
*/
		fpos = ftello(f1);
		if (fpos > progress) {
			fprintf(stderr,"filepos: %llu\n",(unsigned long long)fpos);
			progress += 20000000L;
		}
		for(i = 0 ; i < max1 ; i++) {
			if (i >= result1->size) break;
			if (isint(result1->data[i].string,&testi)) {
				isnum = 1;
			} else if (isdouble(result1->data[i].string,&testd)) {
				isnum = 2;
			} else {
				isnum = 0;
			}
			if (isnum && colinfo[i].numnum >= 256) {
				/* only put < 256 numbers in hashtable, otherwise only min and max data */
			} else if (colinfo[i].hashtable->datasize < 500) {
				bucket = hash_get(colinfo[i].hashtable, (void *)(result1->data+i), hash_DString_hash, hash_DString_compare, &new,1);
				if (new == 0) {
					/* key was already present in the hashtable */
					temp = (long int)hash_getvalue(bucket);
					temp++;
					hash_setvalue(bucket,temp);
				} else {
					DString *buffer = DStringNew();
					DStringSetS(buffer,result1->data[i].string,result1->data[i].size);
					hash_setkey(bucket,buffer);
					hash_setvalue(bucket,1);
					buffer = NULL;
				}
			}
			if (colinfo[i].type == 'i') {
				if (isnum == 1) {
					if (new) {colinfo[i].numnum++;}
					if (testi < colinfo[i].min) {colinfo[i].min = testi;}
					if (testi > colinfo[i].max) {colinfo[i].max = testi;}
				} else if (isnum == 2) {
					colinfo[i].type = 'd';
				} else if (new && colinfo[i].extra != NULL) {
					colinfo[i].type = 's';
					continue;
				} else if (new) {
					colinfo[i].extra = strdup(result1->data[i].string);
				}
			}
			if (colinfo[i].type == 'd') {
				if (isnum == 1) {
					if (new) {colinfo[i].numnum++;}
					if (testi < colinfo[i].min) {colinfo[i].min = testi;}
					if (testi > colinfo[i].max) {colinfo[i].max = testi;}
				} else if (isnum == 2) {
					if (new) {colinfo[i].numnum++;}
					if (testd < colinfo[i].min) {colinfo[i].min = testd;}
					if (testd > colinfo[i].max) {colinfo[i].max = testd;}
				} else if (new && colinfo[i].extra != NULL) {
					colinfo[i].type = 's';
					continue;
				} else if (new) {
					colinfo[i].extra = strdup(result1->data[i].string);
				}
			}
		}
	}
	if (count == -1) {count = 0;}
	fprintf(f2,"%llu\tend\t0\n",(long long int)count);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fprintf(stderr,"Writing columns\n");
	if (line1 != NULL) {DStringDestroy(line1);}
	if (result1 != NULL) {DStringArrayDestroy(result1);}
	cur = argv[9];
	buffer = DStringNew();
	for(i = 0 ; i < max1 ; i++) {
		DStringSet(buffer,prefix);
		size = 0;
		while (cur[size] != ' ' && cur[size] != '\n' &&cur[size] != '\0') {
			size++;
		}
		DStringAppendS(buffer,cur,size);
		DStringAppend(buffer,".colinfo");
		cur++;
		fprintf(stderr,"Writing %s\n",buffer->string);
		cur += size;
		f1 = fopen64_or_die(buffer->string,"w");
		fprintf(f1,"# type: %c\n",colinfo[i].type);
		if (colinfo[i].type == 'i' || colinfo[i].type == 'd') {
			fprintf(f1,"# null: %s\n",colinfo[i].extra);
		}
		if (colinfo[i].min == INFINITY) {
			fprintf(f1,"# min: NaN\n");
		} else {
			fprintf(f1,"# min: %lf\n",colinfo[i].min);
		}
		if (colinfo[i].max == -INFINITY) {
			fprintf(f1,"# max: NaN\n");
		} else {
			fprintf(f1,"# max: %lf\n",colinfo[i].max);
		}
		bucket = hash_first(colinfo[i].hashtable,&iter);
		while(bucket != NULL) {
			DString *ds = hash_getkey(bucket);
			temp = (long int)hash_getvalue(bucket);
			fprintf(f1,"%s\t%ld\n",ds->string,temp);
			DStringDestroy(ds);
			bucket = hash_next(&iter);
		}
		hash_destroy(colinfo[i].hashtable,NULL,NULL);
		if (colinfo[i].hashtable->datasize >= 500 || (colinfo[i].numnum >= 256)) {
			fprintf(f1,"...\n");
		}
		fclose(f1);
		if (colinfo[i].extra != NULL) {free(colinfo[i].extra);}
	}
	DStringDestroy(buffer);
	free(colinfo);
	exit(EXIT_SUCCESS);
}
