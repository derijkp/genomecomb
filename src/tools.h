#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define DSTRING_STATICLEN 5

typedef struct DString {
	int memsize;
	int size;
	char *string;
	char staticspace[DSTRING_STATICLEN];
} DString;

typedef struct DStringArray {
	DString *data;
	int size;
} DStringArray;

typedef struct Buffer {
	int memsize;
	int size;
	char *data;
	char *pos;
} Buffer;

void DStringInit(DString *dstring);
DString *DStringNew();
DString *DStringEmtpy();
void DStringDestroy(DString *dstring);
void DStringClear(DString *dstring);
void DStringSetSize(DString *dstring, int size);
void DStringSet(DString *dstring, char *string);
void DStringSetS(DString *dstring, char *string,int size);
void DStringCopy(DString *dest, DString *src);
int DStringCompare(DString *a, DString *b);
int DStringLocCompare(DString *a, DString *b);
char *Loc_ChrString(DString *ds);
void DStringAppend(DString *dstring, char *string);
void DStringAppendS(DString *dstring, char *string,int size);
ssize_t DStringGetLine(DString *dstring,	FILE *f1);

void InitBuffer(Buffer *buffer,int size);
void DelBuffer(Buffer *buffer);
int DStringGetLine_b(DString *linePtr,	FILE *f1,Buffer *buffer);

DStringArray *DStringArrayNew(int size);
void DStringArrayDestroy(DStringArray *dstringarray);
void check_numfieldserror(int numfields,int numfields2,DString *line,char *filename,unsigned int *linenum);
int DStringGetTab(DString *line,	FILE *f1, int max, DStringArray *result, int setzero,unsigned int *numfields);
FILE *fopen64_or_die(char *filename,char *mode);

int parse_pos(char *arg, int **rresult, int *rnum);

#define FINISHED 1000000
#define NOCHROM 100000
int get_region(FILE *f1, DString *line, int chr1pos, int start1pos, int end1pos, int max1, DString **chromosome1, int *start1, int *end1);
void skip_header(FILE *f1, DString *linePtr,unsigned int *numfields,unsigned int *pos);

int naturalcompare (char const *a, char const *b,int alen,int blen);
int loccompare (char const *a, char const *b,int alen,int blen);
