#ifndef GZTOOLS_LOADED
#define TOOLS_LOADED 1

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <ctype.h>
#include "debug.h"
#include "tools.h"
#include "lz4tools.h"
#include "razf.h"
#include "bgzf.h"
#include "gztools.h"

GZFILE *gz_open(char *filename) {
	GZFILE *result = (GZFILE *)malloc(sizeof(GZFILE));
	int type = UNCOMPRESSED, len = strlen(filename);
	char *cur;
	result->filename = filename;
	if (len == 1 && *filename == '-') {
		type = IN;
	} else if (len >= 3) {
		cur = filename + len - 1;
		if (len >= 4 && *cur == '4') {
			if (*(--cur) == 'z' && *(--cur) == 'l' && *(--cur) == '.') {
				type = LZ4;
			}
		} else if (len >= 4 && *cur == 't') {
			if (*(--cur) == 's' && *(--cur) == 'z' && *(--cur) == '.') {
				type = ZSTD;
			}
		} else if (*cur == 'z') {
			cur--;
			if (*(cur) == 'g') {
				cur--;
				if (*(cur) == '.') {
					type = GZ;
				} else if ((len >=4 && *(cur) == 'b' && *(--cur) == '.')) {
					type = BGZ;
				}
			} else if (*(cur) == 'r' && *(--cur) == '.') {
				type = RZ;
			}
		}
	}
	if (type == IN) {
		result->fun = (FILE *)stdin;
		type = UNCOMPRESSED;
	} else if (type == UNCOMPRESSED) {
		result->fun = fopen64_or_die(filename,"r");
	} else if (type == LZ4) {
		result->lz4 = lz4_openfile(filename,0);
	} else if (type == ZSTD) {
		result->zstd = zstd_openfile(filename);
	} else if (type == GZ) {
		result->gz = gzopen(filename, "r");
	} else if (type == BGZ) {
		result->bgz = bgzf_open(filename, "r");
	} else if (type == RZ) {
		result->rz = razf_open(filename, "r");
	}
	result->type = type;
	return result;
}

void gz_seek(GZFILE *f, int64_t pos, int where) {
	int type = f->type;
	if (type == UNCOMPRESSED) {
		fseeko(f->fun, pos, where);
	} else if (type == LZ4) {
		lz4_seek(f->lz4, pos, where);
	} else if (type == ZSTD) {
		zstd_seek(f->zstd, pos, where);
	} else if (type == GZ) {
		gzseek(f->gz, pos, where);
	} else if (type == BGZ) {
		bgzf_seek(f->bgz, pos, where);
	} else if (type == RZ) {
		razf_seek(f->rz, pos, where);
	}
}

int gz_read(GZFILE *f,void *data,uint64_t size) {
	int type = f->type;
	if (type == UNCOMPRESSED) {
		uint64_t count = size,c;
		unsigned char *buffer = data;
		FILE *f1 = f->fun;
		while (count--) {
			c = getc_unlocked(f1);
			if (c == EOF) break;
			*buffer++ = c;
		}
		return(size-(count+1));
	} else if (type == LZ4) {
		return(lz4_read(f->lz4, data, size));
	} else if (type == ZSTD) {
		return(zstd_read(f->zstd, data, size));
	} else if (type == GZ) {
		return(gzread(f->gz, data, size));
	} else if (type == BGZ) {
		return(bgzf_read(f->bgz, data, size));
	} else {
		return(razf_read(f->rz, data, size));
	}
}

int gz_get(GZFILE *f) {
	int type = f->type;
	if (type == UNCOMPRESSED) {
		return(getc_unlocked(f->fun));
	} else if (type == LZ4) {
		return(lz4_get(f->lz4));
	} else if (type == ZSTD) {
		return(zstd_get(f->zstd));
	} else if (type == GZ) {
		return(gzgetc(f->gz));
	} else if (type == BGZ) {
		return(bgzf_getc(f->bgz));
	} else {
		char c;
		if (!razf_read(f->rz, (void *)(&c), 1)) {return EOF;}
		return(c);
	}
}

void gz_close(GZFILE *f) {
	int type;
	if (f == NULL) return;
	type = f->type;
	if (type == UNCOMPRESSED && f->fun != (FILE *)stdin) {
		FCLOSE(f->fun);
	} else if (type == LZ4) {
		lz4_close(f->lz4);
	} else if (type == ZSTD) {
		zstd_close(f->zstd);
	} else if (type == GZ) {
		gzclose(f->gz);
	} else if (type == BGZ) {
		bgzf_close(f->bgz);
	} else if (type == RZ) {
		razf_close(f->rz);
	}
	free(f);
}


int gz_DStringGetLine(DString *linePtr,	GZFILE *f1) {
	register char *buf,*buf_end;
	int type = f1->type;
	int bsz;
	register char *cur;
	register int cnt;
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
	if (type == UNCOMPRESSED) {
		FILE *f = f1->fun;
		register int	c;
		while ((c=getc_unlocked(f))!=EOF) {
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
		if (c == EOF && cnt == 0) {cnt = -1;}
	} else if (type == LZ4) {
		LZ4res *lz4 = f1->lz4;
		register int c;
		while ((c=lz4_get(lz4))!=EOF) {
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
		if (c == EOF && cnt == 0) {cnt = -1;}
	} else if (type == ZSTD) {
		ZSTDres *zstd = f1->zstd;
		register int c;
		while ((c=zstd_get(zstd))!=EOF) {
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
		if (c == EOF && cnt == 0) {cnt = -1;}
	} else if (type == GZ) {
		gzFile gz = f1->gz;
		register int c;
		while ((c=gzgetc(gz))!=EOF) {
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
		if (c == EOF && cnt == 0) {cnt = -1;}
	} else if (type == BGZ) {
		BGZF *bgz = f1->bgz;
		register int c;
		while ((c=bgzf_getc(bgz))!=EOF) {
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
		if (c == EOF && cnt == 0) {cnt = -1;}
	} else if (type == RZ) {
		RAZF *rz = f1->rz;
		char c;
		while (1) {
			if (!razf_read(rz, (void *)(&c), 1)) break;
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
		if (c == EOF && cnt == 0) {cnt = -1;}
	}
	if (buf_end == cur) {
		bsz += 1;
		buf=(char *)realloc(buf,bsz);
		buf_end=buf+bsz;
	}
	*cur='\0';
	linePtr->string = buf;
	linePtr->size = cnt;
	linePtr->memsize = buf_end-buf;
	return cnt;
}



/*
 * this skips the header (file position will be just after the header line),
 * and returns the header line in linePtr.
 * if numfields != NUL, it will contain the number of fields in the header
 */
void gz_skip_header(GZFILE *f1, DString *linePtr,unsigned int *numfields,unsigned int *pos) {
	ssize_t read;
	unsigned int curpos=0;
	if (pos != NULL) {curpos=*pos;}
	read = gz_DStringGetLine(linePtr, f1);
	if (read == -1) return;
	if (linePtr->size >= 16 && strncmp(linePtr->string,"##fileformat=VCF",16) == 0) {
		/* vcf style header */
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				gz_DStringGetLine(linePtr, f1); curpos++;
				break;
			}
			if (linePtr->string[0] != '#' || linePtr->string[1] != '#') {
				break;
			}
			read = gz_DStringGetLine(linePtr, f1); curpos++;
		}
	} else {
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				gz_DStringGetLine(linePtr, f1); curpos++;
				break;
			}
			if (linePtr->string[0] != '#' && linePtr->string[0] != '>') break;
			read = gz_DStringGetLine(linePtr, f1); curpos++;
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

int gz_DStringGetTab(DString *linePtr,	GZFILE *f1, int maxtab, DStringArray *result, int setzero,unsigned int *numfields) {
	register char *cur = linePtr->string;
	register int c,newdata=0;
	register unsigned int count=0,othertab=0;
	ssize_t size = 0;
	FILE *f = f1->fun;
	LZ4res *lz4 = f1->lz4;
	ZSTDres *zstd = f1->zstd;
	BGZF *bgz = f1->bgz;
	gzFile gz = f1->gz;
	RAZF *rz = f1->rz;
	int type = f1->type;
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
		if (type == UNCOMPRESSED) {
			c = getc_unlocked(f);
		} else if (type == LZ4) {
			c = lz4_get(lz4);
		} else if (type == ZSTD) {
			c = zstd_get(zstd);
		} else if (type == BGZ) {
			c = bgzf_getc(bgz);
		} else if (type == GZ) {
			c = gzgetc(gz);
		} else {
			char tc;
			if (!razf_read(rz, (void *)(&tc), 1)) {c = EOF;} else {c = (int)tc;}
		}
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
	while (count < maxtab) {
		result->data[count++].size = size;
	}
	/* make array */
	{
	register int prevsize;
	count = 0;
	prevsize = 0;
	while (count < maxtab) {
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

#endif
