#ifndef TOOLS_LOADED

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <ctype.h>
#include "debug.h"

#define DSTRING_STATICLEN 5
#define FINISHED 1000000
#define NOCHROM 100000
#define CHROMMAX 9000
#define CHROMM 9095
#define CHROMX 9096
#define CHROMY 9097

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

void DStringInit(DString *dstring) {
	dstring->memsize = DSTRING_STATICLEN-1;
	dstring->size = 0;
	dstring->string = dstring->staticspace;
	dstring->string[0]='\0';
}

DString *DStringNew() {
	DString *dstring = malloc(sizeof(DString));
	DStringInit(dstring);
	return dstring;
}

static DString empty_dstring;

DString *DStringEmtpy() {
	DStringInit(&empty_dstring);
	return &empty_dstring;
}

void DStringClear(DString *dstring) {
	if (dstring->string != dstring->staticspace && dstring->memsize != -1) {
		free(dstring->string);
		dstring->string = dstring->staticspace;
	}
	dstring->memsize = DSTRING_STATICLEN;
	dstring->size = 0;
	dstring->staticspace[0]='\0';
}

void DStringDestroy(DString *dstring) {
	if (dstring == NULL) return;
	if (dstring == &empty_dstring) return;
	DStringClear(dstring);
	free(dstring);
}

void DStringSetSize(DString *dstring, int size) {
	int ssize = dstring->size;
	if (ssize > size) {ssize = size;}
	size++;
	if (dstring->memsize < size) {
		if (dstring->string == dstring->staticspace) {
			dstring->string = malloc(size);
			strncpy(dstring->string,dstring->staticspace,ssize);
			dstring->string[ssize] = '\0';
		} else if (dstring->memsize == -1) {
			char *temp;
			temp = malloc(size);
			strncpy(temp,dstring->string,ssize);
			dstring->string = temp;
			dstring->string[ssize] = '\0';
		} else {
			dstring->string = realloc(dstring->string,size);
		}
		dstring->memsize = size;
	}
}

void DStringSet(DString *dstring, char *string) {
	int size = strlen(string);
	DStringSetSize(dstring,size);
	strncpy(dstring->string,string,size+1);
	dstring->size = size;
}

void DStringSetS(DString *dstring, char *string,int size) {
	DStringSetSize(dstring,size);
	strncpy(dstring->string,string,size+1);
	dstring->string[size] = '\0';
	dstring->size = size;
}

void DStringCopy(DString *dest, DString *src) {
	DStringSetSize(dest,src->size);
	strncpy(dest->string,src->string,src->size+1);
	dest->size = src->size;
}

void DStringputs(DString *string,FILE *f) {
	char *cur;
	int size;
	cur = string->string;
	size = string->size;
	if (size > 0) {
		while(size--) {putc_unlocked(*cur++,f);}
	}
}

void charputs(char *cur,int size,FILE *f) {
	if (size > 0) {
		while(size--) {putc_unlocked(*cur++,stdout);}
	}
}

void DStringPrintf(DString *dstring, char *format, ...) {
	va_list args;
	int size;
	va_start(args, format);
	size = vsnprintf(NULL,0, format, args);
	va_end(args);
	DStringSetSize(dstring,dstring->size+size);
	va_start(args, format);
	vsprintf(dstring->string+dstring->size, format, args);
	dstring->size += size;
	va_end(args);
}

DString *DStringNewFromChar(char *string) {
	DString *dstring = DStringNew();
	DStringSet(dstring,string);
	return dstring;
}

DString *DStringNewFromCharS(char *string,int size) {
	DString *dstring = DStringNew();
	DStringSetS(dstring,string,size);
	return dstring;
}

DString *DStringDup(DString *dstring) {
	DString *result = DStringNew();
	DStringSetS(result,dstring->string,dstring->size);
	return result;
}

DString *DStringNewFromInt(int i) {
	DString *dstring = DStringNew();
	int size=snprintf(NULL,0,"%d",i);
	DStringSetSize(dstring,size);
	sprintf(dstring->string,"%d",i);
	return dstring;
}

/* 
    "natural" sort order: deals with mixed alphabetic and numbers
 */
#define UCHAR(c) ((unsigned char) (c))
#define NATDIGIT(c) (isdigit(UCHAR(*(c))))
#define blank(char) (char == ' ' || char == '\t')
int naturalcompare (char const *a, char const *b,int alen,int blen) {
	int diff, digitleft,digitright;
	int secondaryDiff = 0,prezero,invert,comparedigits,prevdigit=0;
	char *left=NULL,*right=NULL,*keep=NULL;
	while (alen && blank(*a)) {a++; alen--;}
	while (blen && blank(*b)) {b++; blen--;}
	left = (char *)a;
	right = (char *)b;
	/* fprintf(stdout,"---------- %s <> %s ----------\n",a,b);fflush(stdout); */
	while (1) {
		diff = *left - *right;
		/* fprintf(stdout,"diff %c <> %c: %d\n", *left, *right, diff);fflush(stdout); */
		if (!alen || *left == '\0') {
			if (!blen || *right == '\0') {return secondaryDiff;} else {break;}
		}
		if (!blen || *right == '\0') {break;}
		if (diff != 0) {	
			/* only sort on case if no other diff -> keep secondaryDiff for case diff */
			if (isupper(UCHAR(*left)) && islower(UCHAR(*right))) {
				diff = tolower(*left) - *right;
				if (diff) {
					break;
				} else if (secondaryDiff == 0) {
					secondaryDiff = -1;
				}
			} else if (isupper(UCHAR(*right)) && islower(UCHAR(*left))) {
				diff = *left - tolower(UCHAR(*right));
				if (diff) {
					break;
				} else if (secondaryDiff == 0) {
					secondaryDiff = 1;
				}
			} else {
				break;
			}
		}
		left++; alen--;
		right++; blen--;
	}
	if (diff == 0) {return secondaryDiff;}
	if (left == a) {
		prevdigit = 0;
	} else {
		/* previous is the same for left and right! */
		prevdigit = NATDIGIT(right-1);
	}
	digitright = blen && NATDIGIT(right);
	digitleft = alen && NATDIGIT(left);
	/* fprintf(stdout,"%c <> %c    digit %d <> %d    diff %d\n", *left, *right, digitleft, digitright, diff);fflush(stdout); */
	if (*right == '+') {
		if (*left == '-') {
			if (NATDIGIT(left+1) && NATDIGIT(right+1)) {return -1;} else {return 1;}
		} else if (digitleft) {
			if (prevdigit) {return 1;}
			right++; blen--;
			digitright = blen && NATDIGIT(right);
			if (!digitright) {
				return (diff<0)?-1:1;
			}
			if (secondaryDiff == 0) {secondaryDiff = -1;}
		} else {
			return (diff<0)?-1:1;
		}
	} else if (*left == '+') {
		if (*right == '-') {
			if (NATDIGIT(right+1) && NATDIGIT(left+1)) {return 1;} else {return -1;}
		} else if (digitright) {
			if (prevdigit) {return 1;}
			left++; alen--;
			digitleft = alen && NATDIGIT(left);
			if (!digitleft) {
				return (diff<0)?-1:1;
			}
			if (secondaryDiff == 0) {secondaryDiff = 1;}
		} else {
			return (diff<0)?-1:1;
		}
	} else if (*left == '-') {
		if (digitright) {
			if (!prevdigit) {
				return -1;
			} else {
				return 1;
			}
		} else {
			return ((*left - *right)<0)?-1:1;
		}
	} else if (*right == '-' && digitleft) {
		if (digitleft) {
			if (!prevdigit) {
				return 1;
			} else {
				return -1;
			}
		} else {
			return ((*left - *right)<0)?-1:1;
		}
	}
	/* fprintf(stdout,"digit %s <> %s: %d, %d\n",left,right,digitleft,digitright);fflush(stdout); */

	comparedigits = 1;
	if (!digitright) {
		if (!digitleft || (right == b) || !NATDIGIT(right-1)) {comparedigits = 0;}
	}
	if (!digitleft) {
		if ((left == a) || !NATDIGIT(left-1)) {comparedigits = 0;}
	}

	if (!comparedigits) {
		if (diff == 0) {return secondaryDiff;} else {return (diff<0)?-1:1;}
	}
	/*
	 * There are decimal numbers embedded in the two
	 * strings.  Compare them as numbers, rather than
	 * strings.  If one number has more leading zeros than
	 * the other, the number with more leading zeros sorts
	 * later, but only as a secondary choice.
	 */
	/* first take steps back (keep) to check for -0.0 */
	keep = left;
	prezero = 1;
	while (keep > a) {
		keep--;
		if (!NATDIGIT(keep)) break;
		if (*keep != '0') prezero = 0;
	}
	if (*keep == '.') {
		while (keep > a) {
			keep--;
			if (!NATDIGIT(keep)) break;
		}
		if (*keep == '-') {
			invert = -1;
		} else {
			invert = 1;
		}
		if (diff == 0) {return invert*secondaryDiff;} else {return invert*diff;}
	} else {
		if (*keep == '-') {
			invert = -1;
		} else {
			if (NATDIGIT(keep) && (*keep != '0')) {prezero = 0;}
			invert = 1;
		}
		if (prezero) {
			if (*left == '0') {
				secondaryDiff = 1;
				while (*left == '0') {
					left++; alen--;
				}
			} else if (*right == '0') {
				secondaryDiff = -1;
				while (*right == '0') {
					right++; blen--;
				}
			}
		}
		/*
		 * The code below compares the numbers in the two
		 * strings without ever converting them to integers.  It
		 * does this by first comparing the lengths of the
		 * numbers and then comparing the digit values.
		 */
		diff = 0;
		while (1) {
			if (diff == 0) {
				diff = *left - *right;
			}
			if (!blen || !NATDIGIT(right)) {
				if (alen && NATDIGIT(left)) {
					return invert;
				} else {
					/*
					 * The two numbers have the same length. See
					 * if their values are different.
					 */
	
					if (diff != 0) {
						return invert*diff;
					} else {
						return invert*secondaryDiff;
					}
				}
			} else if (!alen || !NATDIGIT(left)) {
				return -1*invert;
			}
			right++; blen--;
			left++; alen--;
		}
	}
	fprintf(stderr,"shouldnt get here\n"); exit(1);
	return 0;
}

int DStringCompare(DString *a, DString *b) {
	if (a == b) {return 0;}
	if (a == NULL) {
		return 1;
	} else if (b == NULL) {
		return -1;
	}
	return naturalcompare(a->string,b->string,a->size,b->size);
}

int loccompare (char const *a, char const *b,int alen,int blen) {
	if (alen >= 3) {
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			a += 3; alen -= 3;
			if (alen && a[0] == '-') {
				a++; alen--;
			}
		}
	}
	if (blen >= 3) {
		if ((b[0] == 'C' || b[0] == 'c') && (b[1] == 'H' || b[1] == 'h') && (b[2] == 'R' || b[2] == 'r')) {
			b += 3; blen -= 3;
			if (blen && b[0] == '-') {
				b++; blen--;
			}
		}
	}
	return naturalcompare(a,b,alen,blen);
}

int DStringLocCompare(DString *a, DString *b) {
	if (a == b) {return 0;}
	if (a == NULL) {
		return 1;
	} else if (b == NULL) {
		return -1;
	}
	return loccompare(a->string,b->string,a->size,b->size);
}

char *Loc_ChrString(DString *ds) {
	char *s;
	if (ds == NULL) {
		return "";
	} else {
		s = ds->string;
		if (ds->size > 3) {
			if ((s[0] == 'C' || s[0] == 'c') && (s[1] == 'H' || s[1] == 'h') && (s[2] == 'R' || s[2] == 'r')) {
				if (ds->size > 4 && s[3] == '-') {
					s += 4;
				} else {
					s += 3;
				}
			}
		}
		return s;
	}
}

void DStringAppend(DString *dstring, char *string) {
	int size = strlen(string);
	int nsize = dstring->size + size;
	DStringSetSize(dstring,nsize);
	strncpy(dstring->string+dstring->size,string,size+1);
	dstring->size = nsize;
}

void DStringAppendS(DString *dstring, char *string,int size) {
	int nsize = dstring->size + size;
	DStringSetSize(dstring,nsize);
	strncpy(dstring->string+dstring->size,string,size+1);
	dstring->string[nsize] = '\0';
	dstring->size = nsize;
}

int DStringGetLine(DString *linePtr,	FILE *f1) {
	register char *buf,*buf_end;
	int bsz;
	register char *cur;
	register int cnt;
	register int	c;
	if (linePtr->memsize != -1 && linePtr->string != linePtr->staticspace) {
		buf = linePtr->string;
		buf_end = linePtr->string+linePtr->memsize;
		bsz = linePtr->memsize;
	} else {
		buf=malloc(128);
		bsz=128;
		buf_end=buf+bsz;
	}
	cur=buf;
	cnt=0;
	while ((c=getc_unlocked(f1))!=EOF) {
		if (c == '\n') break;
		*cur++=c;
		cnt++;
		if (cur == buf_end) {
			buf=realloc(buf,bsz*2);
			cur=buf+bsz;
			bsz*=2;
			buf_end=buf+bsz;
		}
	}
	*cur='\0';
	linePtr->string = buf;
	linePtr->size = cnt;
	linePtr->memsize = buf_end-buf;
	if (c == EOF && cnt == 0) return -1;
	return cnt;
}

void SkipLine(FILE *f1) {
	register int c;
	while ((c=getc_unlocked(f1))!=EOF) {
		if (c == '\n') break;
	}
}

void InitBuffer(Buffer *buffer,int size) {
	buffer->memsize = size;
	buffer->data = (char *)malloc(size+1);
	buffer->pos = buffer->data;
	buffer->size = 0;
}

void DelBuffer(Buffer *buffer) {
	buffer->memsize = 0;
	free(buffer->data);
	buffer->data = NULL;
	buffer->pos = NULL;
	buffer->size = 0;
}

int DStringGetLine_b(DString *linePtr,	FILE *f1,Buffer *buffer) {
	char *cur = linePtr->string;
	int c,newdata=0;
	ssize_t size = 0;
NODPRINT("getline");
	while (1) {
		if (!buffer->size) {
			buffer->size = fread(buffer->data,sizeof(char),buffer->memsize,f1);
			buffer->pos = buffer->data;
NODPRINT("read buffer %d",buffer->size);
			if (!buffer->size) break;
		}
		newdata = 1;
		buffer->size--;
		c = *(buffer->pos)++;
		if  (c == '\n') {
			break;
		}
		if (linePtr->memsize <= size) {
			linePtr->size = size;
			DStringSetSize(linePtr,2*linePtr->memsize);
			cur = linePtr->string+size;
		}
		*cur++ = c;
		++size;
	}
	if (!newdata) {
		*cur = '\0';
		return -1;
	}
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,2*linePtr->memsize);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	return size;
}

DStringArray *DStringArrayNew(int size) {
	DStringArray *dstringarray;
	int i;
	dstringarray = (DStringArray *)malloc(1*sizeof(DStringArray));
	dstringarray->data = (DString *)malloc(size*sizeof(DString));
	dstringarray->size = 0;
	dstringarray->memsize = size;
	dstringarray->datablock = NULL;
	for (i = 0; i < size ; i++) {
		DStringInit(dstringarray->data+i);
	}
	return dstringarray;
}

DStringArray *DStringArrayAppend(DStringArray *dstringarray,char *string,int size) {
	if (dstringarray->size == dstringarray->memsize) {
		DString *temp = NULL;
		int oldmemsize = dstringarray->memsize,i;
		dstringarray->memsize *= 2;
		/* 
			need to correct string pointer in all DStrings that point to the static space
		*/
		temp = (DString *)malloc(dstringarray->memsize*sizeof(DString));
		memcpy(temp,dstringarray->data,oldmemsize*sizeof(DString));
		for (i = 0 ; i < dstringarray->size ; i++) {
			if (dstringarray->data[i].string == dstringarray->data[i].staticspace) {
				temp[i].string = temp[i].staticspace;
			}
		}
		free(dstringarray->data);
		dstringarray->data = temp;
	}
	DStringInit(dstringarray->data+dstringarray->size);
	if (size < 0) {size = strlen(string);}
	DStringSetS(dstringarray->data+dstringarray->size, string, size);
	dstringarray->size++;
	return dstringarray;
}

DStringArray *DStringArrayFromChar(char *string,char sep) {
	DStringArray *result;
	char *cur,*prev;
	int count=1;
	cur=string;
	while(*cur) {
		if (*cur == sep) count++;
		cur++;
	}
	result = DStringArrayNew(count+1);
	cur = string;
	prev = cur;
	while(*cur) {
		if (*cur == sep) {
			DStringArrayAppend(result,prev,cur-prev);
			prev = cur+1;
		}
		cur++;
	}
	DStringArrayAppend(result,prev,cur-prev);
	return result;
}

DStringArray *DStringArrayFromCharM(char *string,char *seps) {
	DStringArray *result;
	char *cur,*prev;
	result = DStringArrayNew(5);
	cur = string;
	prev = cur;
	while(*cur) {
		if (strchr((const char *)seps,*cur)) {
			DStringArrayAppend(result,prev,cur-prev);
			do {cur++;} while(strchr((const char *)seps,*cur) && *cur);
			prev = cur;
			if (!*cur) break;
		} else {
			cur++;
		}
	}
	if (cur-prev > 0) {
		DStringArrayAppend(result,prev,cur-prev);
	}
	return result;
}

DStringArray *DStringArrayRange(DStringArray *dstringarray,int start, int end) {
	DStringArray *result;
	int i;
	result = DStringArrayNew(end-start+1);
	for (i = start; i <= end ; i++) {
		DStringCopy(result->data+(i-start), dstringarray->data+i);
	}
	result->size = end-start+1;
	return result;	
}

DStringArray *DStringArraySet(DStringArray *dstringarray,int pos,char *string,int size) {
	int i;
	if (pos >= dstringarray->memsize) {
		dstringarray->memsize = pos+1;
		dstringarray->data = (DString *)realloc(dstringarray->data,dstringarray->memsize*sizeof(DString));
		i = dstringarray->size;
		while (i <= pos) {
			DStringInit(dstringarray->data+i);
			i++;
		}
		dstringarray->size=pos+1;
	}
	if (size < 0) {size = strlen(string);}
	DStringSetS(dstringarray->data+pos, string, size);
	return dstringarray;
}

int DStringArraySearch(DStringArray *dstringarray,char *string,int size) {
	DString *astring;
	int i;
	for (i = 0; i < dstringarray->size ; i++) {
		astring = DStringArrayGet(dstringarray,i);
		if (size != -1) {
			if (astring->size == size && strncmp(astring->string,string,size) == 0) {
				return i;
			}
		} else {
			if (strcmp(astring->string,string) == 0) {
				return i;
			}
		}
	}
	return -1;
}

void DStringArrayDestroy(DStringArray *dstringarray) {
	if (dstringarray == NULL) return;
	int i=0;
	for (i =0; i < dstringarray->memsize ; i++) {
		DStringClear(dstringarray->data+i);
	}
	free(dstringarray->data);
	if (dstringarray->datablock) DStringDestroy(dstringarray->datablock);
	free(dstringarray);
}

void check_numfieldserror(int numfields,int numfields2,DString *line,char *filename,unsigned int *linenum) {
	if (numfields != numfields2) {
		fprintf(stderr,"ERROR in file %s",filename);
		if (linenum != NULL) {
			fprintf(stderr," at line %d",*linenum);
		}
		fprintf(stderr,": number of columns different from the header: %d instead of %d for:\n%s\n",numfields,numfields2,line->string);
		exit(1);
	}
}

int DStringGetTab(
	DString *linePtr,	FILE *f1, int maxtab, DStringArray *result, int setzero,unsigned int *numfields
) {
	register char *cur = linePtr->string;
	register int c,newdata=0;
	register unsigned int count=0,othertab=0;
	ssize_t size = 0;
NODPRINT("maxtab=%d result->memsize=%d",maxtab,result->memsize)
	maxtab += 1;
	if (maxtab > result->memsize) {
		fprintf(stderr,"cannot DStringGetTab, size of allocated DStringArray must be >= maxtab+2\n");
		exit(1);
	}
	result->data[count].string = cur;
	result->data[count].memsize = -1;
	/* fill current pos (size) in the array for now, convert to array element sizes later */
	while (1) {
		c = getc_unlocked(f1);
		if  (c == '\t') {
			if (count < maxtab) {
				result->data[count].size = size;
				/* fprintf(stdout,"count=%d size=%d\n",count, result->data[count].size); */
				count++;
			} else {
				othertab++;
			}
		} else if (c == '\n' || c == EOF) {
			if (count < maxtab) {
				result->data[count].size = size;
			}
			break;
		}
		newdata = 1;
		if (linePtr->memsize <= (size+1)) {
			linePtr->size = size;
			DStringSetSize(linePtr,2*linePtr->memsize);
			cur = linePtr->string+size;
		}
		*cur++ = c;
		++size;
	}
	*cur = '\0';
	if (numfields != NULL) {
		*numfields = count+othertab+1;
	}
	if (count < maxtab) {
		count++;
		result->size = count;
	} else {
		result->size = maxtab;
	}
	/* add \0 to line*/
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,linePtr->memsize+2);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	if (!newdata) {
		return -1;
	}
	/* fill rest of tab separated elements in array */
	while (count <= maxtab) {
		result->data[count++].size = size;
	}
	/* make array */
	{
	register int prevsize;
	count = 0;
	prevsize = 0;
	while (count <= maxtab) {
		int pos = result->data[count].size;
		result->data[count].size = result->data[count].size-prevsize;
		result->data[count].string = linePtr->string + prevsize;
		if (setzero) {
			result->data[count].string[result->data[count].size] = '\0';
		}
		result->data[count].memsize = -1;
NODPRINT("final count=%d size=%d %s",count, result->data[count].size,result->data[count].string)
		prevsize = pos+1;
		count++;
	}
	}
	NODPRINT("\n==== result ====\n%d %s",size,linePtr->string);
	NODPRINT("==================== getline finished ====================");
	return 0;
}

void DStringPrintTab(FILE *f, DString *linePtr) {
	int i;
	for(i = 0 ; i < linePtr->size ; i++) {
		if (linePtr->string[i] == '\0') {fputc('\t',stderr);} else {fputc(linePtr->string[i],stderr);}
	}
}

int DStringSplitTab(
	DString *linePtr,	int maxtab, DStringArray *result, int setzero,unsigned int *numfields
) {
	register char *cur = linePtr->string;
	register int c;
	register unsigned int count=0,othertab=0;
	ssize_t size = 0;
NODPRINT("maxtab=%d result->memsize=%d",maxtab,result->memsize)
	maxtab += 1;
	if (maxtab > result->memsize) {
		fprintf(stderr,"cannot DStringGetTab, array memsize < maxtab+1\n");
		exit(1);
	}
	result->data[count].string = cur;
	result->data[count].memsize = -1;
	/* fill current pos (size) in the array for now, convert to array element sizes later */
	while (1) {
		c = *cur;
		if  (c == '\t') {
			if (count < maxtab) {
				result->data[count].size = size;
				/* fprintf(stdout,"count=%d size=%d\n",count, result->data[count].size); */
				count++;
			} else {
				othertab++;
			}
		} else if (size >= linePtr->size) {
			if (count < maxtab) {
				result->data[count].size = size;
			}
			break;
		}
		cur++;
		++size;
	}
	if (numfields != NULL) {
		*numfields = count+othertab+1;
	}
	if (count < maxtab) {
		count++;
		result->size = count;
	} else {
		result->size = maxtab;
	}
	/* fill rest of tab separated elements in array */
	while (count <= maxtab) {
		result->data[count++].size = size;
	}
	/* make array */
	{
	register int prevsize;
	count = 0;
	prevsize = 0;
	while (count <= maxtab) {
		int pos = result->data[count].size;
		result->data[count].size = result->data[count].size-prevsize;
		result->data[count].string = linePtr->string + prevsize;
		if (setzero) {
			result->data[count].string[result->data[count].size] = '\0';
		}
		result->data[count].memsize = -1;
NODPRINT("final count=%d size=%d %s",count, result->data[count].size,result->data[count].string)
		prevsize = pos+1;
		count++;
	}
	}
	NODPRINT("\n==== result ====\n%d %s",size,linePtr->string);
	NODPRINT("==================== getline finished ====================");
	return 0;
}

int parse_pos(char *arg, int **rresult, int *rnum) {
	char *pch;
	int *result;
	int num,memsize = 5;
	result = (int *)malloc(memsize*sizeof(int));
	pch = strtok (arg, " \t,-");
	num = 0;
	while (pch != NULL) {
		if (*pch == 'x') {
			result[num++] = -1;
		} else {
			result[num++] = atoi(pch);
		}
		if (num >= memsize) {
			memsize += memsize;
			result = (int *)realloc(result, memsize*sizeof(int));
		}
		pch = strtok (NULL, " \t,-");
	}
	*rnum = num;
	*rresult = result;
	return 0;
}

int get_region(FILE *f1, DString *linePtr, int chr1pos, int start1pos, int end1pos, int max1, DString **chromosome1, int *start1, int *end1) {
	char *linepos = NULL,*scanpos = NULL,*endpos;
	ssize_t read;
	int count;
	if (chromosome1 == NULL) {
		return 1;
	}
	if (f1 == NULL) {
		DStringDestroy(*chromosome1);
		*chromosome1 = NULL;
		return 1;
	}
	while ((read = DStringGetLine(linePtr, f1)) != -1) {
NODPRINT("%s",linePtr->string)
		linepos = linePtr->string;
		scanpos = linePtr->string;
		endpos = linePtr->string+linePtr->size;
		count = 0;
		while (1) {
			if (*linepos == '\t' || *linepos == '\0' || linepos == endpos) {
				if (count == chr1pos) {
					char *pos=scanpos;
					while (*pos != '\t' && *pos != '\0' && pos != endpos) {++pos;}
					DStringSetS(*chromosome1,scanpos,pos-scanpos);
				} else if (count == start1pos) {
					sscanf(scanpos,"%d",start1);
				} else if (count == end1pos) {
					sscanf(scanpos,"%d",end1);
				}
				if (count >= max1) break;
				count++;
				scanpos = linepos+1;
			}
			if (*linepos == '\0') break;
			if (linepos == endpos) break;
			linepos++;
		}
		if (count >= max1) break;
	}
	if (read == -1) {
		DStringDestroy(*chromosome1);
		*chromosome1 = NULL;
		return 1;
	} else {
		return 0;
	}
}

void skip_header(FILE *f1, DString *linePtr, unsigned int *numfields,unsigned int *pos) {
	ssize_t read;
	unsigned int curpos=0;
	if (pos != NULL) {curpos=*pos;}
	read = DStringGetLine(linePtr, f1);
	if (read == -1) return;
	if (strlen(linePtr->string) >= 16 && strncmp(linePtr->string,"##fileformat=VCF",16) == 0) {
		/* vcf style header */
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				DStringGetLine(linePtr, f1); curpos++;
				break;
			}
			if (linePtr->string[0] != '#' || linePtr->string[1] != '#') {
				break;
			}
			read = DStringGetLine(linePtr, f1); curpos++;
		}
	} else {
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				DStringGetLine(linePtr, f1); curpos++;
				break;
			}
			if (linePtr->string[0] != '#' && linePtr->string[0] != '>') break;
			read = DStringGetLine(linePtr, f1); curpos++;
		}
	}
	if (pos != NULL) {*pos = curpos;}
	if (numfields != NULL) {
		char *buffer;
		unsigned int count;
		buffer = linePtr->string;
		count = 1;
		while (*buffer != '\0') {
			if (*buffer == '\t') {count++;}
			buffer++;
		}
		*numfields = count;
	}
	NODPRINT("%s",linePtr->string)
}

FILE *fopen64_or_die(char *filename,char *mode) {
	FILE *f;
	f = fopen64(filename,mode);
	if (f == NULL) {
		fprintf(stderr,"Error opening file %s: %s.\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return(f);
}

int checksort(DString *prevchromosome1,int *prevstart1,int *prevend1,DString *prevtype1,DString *prevalt1,DString *chromosome1,int start1,int end1,DString *type1,DString *alt1,char *filename,int *nextpos) {
	int comp,comptype,compalt, result = 2;
 	comp = DStringLocCompare(chromosome1,prevchromosome1);
	comptype = DStringCompare(type1,prevtype1);
	compalt = DStringCompare(alt1,prevalt1);
	if (comp < 0 || (comp == 0 && 
		(start1 < *prevstart1 || (start1 == *prevstart1 && 
		(end1 < *prevend1 || (end1 == *prevend1 &&
		(comptype < 0 || (comptype == 0 && compalt < 0)
	))))))) {
		fprintf(stderr,"File (%s) is not correctly sorted (sort correctly using \"cg select -s -\")\n",filename);
		fprintf(stderr,"%s:%d-%d:%s:%s came before %s:%d-%d:%s:%s\n",prevchromosome1->string,*prevstart1,*prevend1,prevtype1->string,prevalt1->string, chromosome1->string,start1,end1,type1->string,alt1->string);
		exit(1);
	} else if (comp > 0) {
		/* prevchromosome1 = chromosome1; */
		DStringCopy(prevchromosome1,chromosome1);
		if (nextpos) {*nextpos = 0;}
	}
	if (comp == 0 && comptype == 0 && start1 == *prevstart1 && end1 == *prevend1) {
		if (compalt == 0) {
			result = 0;
		} else {
			result = 1;
		}
	}
	*prevstart1 = start1; *prevend1 = end1;
	if (comptype != 0) {DStringCopy(prevtype1,type1);}
	if (compalt != 0) {DStringCopy(prevalt1,alt1);}
	return result;
}

/* checksortreg(prevchromosome1,prevstart1,prevend1,chromosome1,start1,end1,argv[1]); */
int checksortreg(DString *prevchromosome,int *prevstart,int *prevend,DString *chromosome,int start,int end,char *file) {
	int comp;
	if (!prevchromosome || !chromosome) {
		return 0;
	}
 	comp = DStringLocCompare(chromosome,prevchromosome);
	if (comp < 0 || (comp == 0 && (start < *prevstart || (start == *prevstart && end < *prevend)))) {
		fprintf(stderr,"File (%s) is not correctly sorted (sort correctly using \"cg select -s -\")\n",file);
		fprintf(stderr,"%*.*s:%d-%d came before %*.*s:%d-%d\n",
			prevchromosome->size,prevchromosome->size,prevchromosome->string,
			*prevstart,*prevend,
			chromosome->size,chromosome->size,chromosome->string,
			start,end
		);
		exit(1);
	}
	*prevstart = start; *prevend = end;
	if (comp > 0) {
		DStringCopy(prevchromosome,chromosome);
		return 1;
	} else {
		return 0;
	}
}
#endif
#define TOOLS_LOADED 1