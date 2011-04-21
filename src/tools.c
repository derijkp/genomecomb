#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "debug.h"

#define DSTRING_STATICLEN 5
#define FINISHED 1000000
#define NOCHROM 100000
#define CHROMMAX 9000
#define CHROMM 9095
#define CHROMX 9096
#define CHROMY 9097

typedef struct DString {
	int memsize;
	int size;
	char *string;
	char staticspace[DSTRING_STATICLEN];
} DString;

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
	dstring->size = nsize;
}

int DStringGetLine(DString *linePtr,	FILE *f1) {
	char *cur = linePtr->string;
	int c;
	ssize_t size = 0;
	while (1) {
		c = getc(f1);
		if ((c == EOF)&&(size == 0)) {return -1;}
		if  (c == '\n') break;
		if (linePtr->memsize <= size) {
			linePtr->size = size;
			DStringSetSize(linePtr,2*linePtr->memsize);
			cur = linePtr->string+size;
		}
		*cur++ = c;
		++size;
	}
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,2*linePtr->memsize);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	return size;
}

DString *DStringArrayNew(int size) {
	DString *dstringarray;
	int i;
	dstringarray = (DString *)malloc((size+2)*sizeof(DString));
	for (i = 0; i < size ; i++) {
		DStringInit(dstringarray+i);
	}
	dstringarray[i].string=NULL;
	return dstringarray;
}

void DStringArrayDestroy(DString *dstringarray) {
	int i=0;
	while (1) {
		if (dstringarray[i].string == NULL) break;
		DStringClear(dstringarray+i);
		i++;
	}
	free(dstringarray);
}

int DStringGetTab(DString *line,	FILE *f1, int max, DString *result) {
	char *linepos = NULL,*prevpos = NULL;
	ssize_t read;
	int count=0;
	int c;
	while ((read = DStringGetLine(line, f1)) != -1) {
		prevpos = line->string;
		linepos = line->string;
		count = 0;
		while (1) {
			c = *linepos;
			if (c == '\0') {
				result[count].string = prevpos;
				result[count].size = linepos-prevpos;
				result[count].memsize = -1;
				count++;
				break;
			} else if (c == '\n') {
				*linepos = '\0';
				result[count].string = prevpos;
				result[count].size = linepos-prevpos;
				result[count].memsize = -1;
				count++;
				break;
			} else if (c == '\t') {
				*linepos = '\0';
				result[count].string = prevpos;
				result[count].size = linepos-prevpos;
				result[count].memsize = -1;
				prevpos = linepos;
				count++;
				if (count > max) break;
				prevpos = linepos+1;
			}
			if (c == '\0') break;
			linepos++;
		}
		if (count >= max) break;
	}
	if (count < max) {return 1;}
	if (c == '\0') {
		DStringSetS(result+count,"",0);
	} else {
		result[count].string = linepos+1;
		result[count].size = (line->string+line->size) - linepos;
		result[count].memsize = -1;
	}
	if (read == -1) {return 1;}
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
		result[num++] = atoi(pch);
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

int chromosomenum(char *chromosome) {
	int i;
	if (strncmp(chromosome,"chr",3) == 0) {
		chromosome += 3;
	} else if (strncmp(chromosome,"CHR",3) == 0) {
		chromosome += 3;
	}
	if (*chromosome == 'M') {
		return CHROMM;
	} else if (*chromosome == 'X') {
		return CHROMX;
	} else if (*chromosome == 'Y') {
		return CHROMY;
	} else {
		i = atoi(chromosome);
		if (i > CHROMMAX) {
			fprintf(stderr,"wrong chromosome %s",chromosome);
			exit(EXIT_FAILURE);
		}
		return i;
	}
}

char* num2chromosome(char *chromosome,int num) {
	if (num == CHROMM) {
		sprintf(chromosome,"M");
	} else if (num == CHROMX) {
		sprintf(chromosome,"X");
	} else if (num == CHROMY) {
		sprintf(chromosome,"Y");
	} else {
		sprintf(chromosome,"%d",num);
	}
	return chromosome;
}

int get_region(FILE *f1, DString *linePtr, int chr1pos, int start1pos, int end1pos, int max1, DString *chromosome1, int *nchr1, int *start1, int *end1) {
	char *linepos = NULL,*scanpos = NULL,*endpos;
	ssize_t read;
	int count;
	if (f1 == NULL) {
		*nchr1 = FINISHED;
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
					DStringSetS(chromosome1,scanpos,pos-scanpos);
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
		*nchr1 = FINISHED;
		return 1;
	} else {
		*nchr1 = chromosomenum(chromosome1->string);
		return 0;
	}
}

void skip_header(FILE *f1, DString *linePtr) {
	ssize_t read;
	read = DStringGetLine(linePtr, f1);
	if (read == -1) return;
	if (linePtr->string[0] == '#' && linePtr->string[1] == '#') {
		/* vcf style header */
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				DStringGetLine(linePtr, f1);
				break;
			}
			if (linePtr->string[0] != '#' || linePtr->string[1] != '#') {
				break;
			}
			read = DStringGetLine(linePtr, f1);
		}
	} else {
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				DStringGetLine(linePtr, f1);
				break;
			}
			if (linePtr->string[0] != '#' && linePtr->string[0] != '>') break;
			read = DStringGetLine(linePtr, f1);
		}
	}
	NODPRINT("%s",linePtr->string)
}
