#ifndef TOOLS_H_LOADED
#define TOOLS_H_LOADED 1

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define DSTRING_STATICLEN 5

#define DStringArrayGet(dstringarray,pos) (dstringarray->data+pos)

typedef struct DString {
	int memsize;
	int size;
	char *string;
	char staticspace[DSTRING_STATICLEN];
} DString;

typedef struct DStringArray {
	DString *data;
	int size;
	int memsize;
	DString *datablock;
} DStringArray;

typedef struct Buffer {
	int memsize;
	int size;
	char *data;
	char *pos;
} Buffer;

/* $Format: "#define GENOMECOMB_VERSION \"$ProjectMajorVersion$.$ProjectMinorVersion$\""$ */
#define GENOMECOMB_VERSION "0.100"
#define FILEVERSION "0.10.0"
void DStringInit(DString *dstring);
DString *DStringNew();
DString *DStringNewFromChar(char *string);
DString *DStringNewFromCharS(char *string,int size);
DString *DStringNewFromInt(int i);
DString *DStringEmtpy();
DString *DStringDup(DString *dstring);
void DStringDestroy(DString *dstring);
void DStringClear(DString *dstring);
void DStringSetSize(DString *dstring, int size);
void DStringSet(DString *dstring, char *string);
void DStringSetS(DString *dstring, char *string,int size);
void DStringPrintf(DString *dstring, char *format, ...);
void DStringCopy(DString *dest, DString *src);
void DStringputs(DString *string,FILE *f);
void DStringArrayPuts(DStringArray *array,char *join,FILE *f);
void charputs(char *cur,int size,FILE *f);
int DStringCompare(DString *a, DString *b);
int DStringLocCompare(DString *a, DString *b);
char *Loc_ChrString(DString *ds);
void DStringAppend(DString *dstring, char *string);
void DStringAppendS(DString *dstring, char *string,int size);
int DStringGetLine(DString *dstring,	FILE *f1);
void SkipLine(FILE *f1);

void InitBuffer(Buffer *buffer,int size);
void DelBuffer(Buffer *buffer);
int DStringGetLine_b(DString *linePtr,	FILE *f1,Buffer *buffer);

DStringArray *DStringArrayNew(int size);
DStringArray *DStringArrayFromChar(char *string,char sep);
DStringArray *DStringArrayFromCharM(char *string,char *seps);
DStringArray *DStringArrayAppend(DStringArray *dstringarray,char *string,int size);
DStringArray *DStringArraySet(DStringArray *dstringarray,int pos,char *string,int size);
int DStringArraySearch(DStringArray *dstringarray,char *string,int size);
DStringArray *DStringArrayRange(DStringArray *dstringarray,int start, int end);
#define DStringArrayGet(dstringarray,pos) (dstringarray->data+pos)
void DStringArrayDestroy(DStringArray *dstringarray);
void check_numfieldserror(int numfields,int numfields2,DString *line,char *filename,unsigned int *linenum);
int DStringGetTab(DString *line,	FILE *f1, int max, DStringArray *result, int setzero,unsigned int *numfields);
void DStringPrintTab(FILE *f, DString *linePtr);
int DStringSplitTab(DString *line, int max, DStringArray *result, int setzero,unsigned int *numfields);
FILE *fopen64_or_die(char *filename,char *mode);
int checksort(DString *prevchromosome1,int *prevstart1,int *prevend1,DString *prevtype1,DString *prevalt1,DString *chromosome1,int start1,int end1,DString *type1,DString *alt1,char *filename,int *nextpos,int fillprev);
int checksortreg(DString *prevchromosome1,int *prevstart1,int *prevend1,DString *chromosome1,int start1,int end1,char *file);

int parse_pos(char *arg, int **rresult, int *rnum);

#define FINISHED 1000000
#define NOCHROM 100000
int get_region(FILE *f1, DString *line, int chr1pos, int start1pos, int end1pos, int max1, DString **chromosome1, int *start1, int *end1);
void skip_header(FILE *f1, DString *linePtr,unsigned int *numfields,unsigned int *pos);

int naturalcompare (char const *a, char const *b,int alen,int blen);
int loccompare (char const *a, char const *b,int alen,int blen);

#endif
