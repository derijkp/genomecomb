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

void DStringInit(DString *dstring);
void DStringClear(DString *dstring);
void DStringSetSize(DString *dstring, int size);
void DStringSet(DString *dstring, char *string);
void DStringSetS(DString *dstring, char *string,int size);
void DStringCopy(DString *dest, DString *src);

void DStringAppend(DString *dstring, char *string);
ssize_t DStringGetLine(DString *dstring,	FILE *f1);

int parse_pos(char *arg, int **rresult, int *rnum);
int chromosomenum(char *chromosome);

#define FINISHED 1000000

int get_region(FILE *f1, DString *line, int chr1pos, int start1pos, int end1pos, int max1, DString *chromosome1, int *nchr1, int *start1, int *end1);
