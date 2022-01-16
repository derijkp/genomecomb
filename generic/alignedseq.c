/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <limits.h>
#include <inttypes.h>
#include "debug.h"
#include "tcl.h"
#include "tools.h"
#include "genomecomb.h"

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

#define AMPLICONSBUFFERSIZE 40

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract suplementary alignment (part of chimeric alignment) */
#define BAM_FSUPPL      2048

typedef struct Amplicon {
	DString *chr2;
	int start2;
	int end2;
	DString *name;
} Amplicon;

typedef struct Cigar {
	int size;
	int memsize;
	int *num;
	char *action;
} Cigar;

int cigar_refsize(Cigar *cigar) {
	int count = cigar->size;
	int *num = cigar->num;
	char *action = cigar->action;
	int result = 0;
	while(count--) {
		if (*action == 'M' || *action == 'D' || *action == 'N' || *action == '=' || *action == 'X') {
			result += *num;
		}
		action++;
		num++;
	}
	return result;
}

int cigar_clipleft(Cigar *cigar) {
	int *num = cigar->num;
	char *action = cigar->action;
	if (*action == 'H') {
		action++;
		num++;
	}
	if (*action == 'S') {
		return *num;
	} else {
		return 0;
	}
}

int cigar_clipright(Cigar *cigar) {
	int count = cigar->size;
	int *num = cigar->num+count-1;
	char *action = cigar->action+count-1;
	if (*action == 'H') {
		action--;
		num--;
	}
	if (*action == 'S') {
		return *num;
	} else {
		return 0;
	}
}

int parse_cigar(Cigar *cigar,char *string) {
	if (cigar->memsize == 0) {
		cigar->size = 0;
		cigar->memsize = 5;
		cigar->num = malloc(cigar->memsize*sizeof(int));
		cigar->action = malloc(cigar->memsize*sizeof(char));
	}
	cigar->size = 0;
	while (*string && *string != '\t') {
		if (cigar->size >= cigar->memsize) {
			cigar->memsize += 5;
			cigar->num = realloc(cigar->num,cigar->memsize*sizeof(int));
			cigar->action = realloc(cigar->action,cigar->memsize*sizeof(char));
		}
		cigar->num[cigar->size] = (int)strtol(string,&string,10);
		cigar->action[cigar->size++] = *string++;
	}
	return(cigar->size);
}

int cigar_refend(Cigar *cigar,int begin) {
	int count = cigar->size;
	int *num = cigar->num;
	char *action = cigar->action;
	while(count--) {
		if (*action == 'M' || *action == 'D' || *action == 'N' || *action == '=' || *action == 'X') {
			begin += *num;
		}
		action++;
		num++;
	}
	return begin;
}

/* gives coordinate on query, 0 based (does not know if end, so not "half-open") */
int cigar_ref2seq(Cigar *cigar,int begin, int pos,int end_ins) {
	int count = cigar->size;
	int *num = cigar->num, cur = 0, prev = 0, prevbegin = begin;
	char *action = cigar->action;
/*
fprintf(stderr,"begin=%d pos=%d\n",begin,pos);
fflush(stderr);
*/
	/* clip at start counts for position in sequence */
	if (*action == 'S' || *action == 'H') {
		cur += *num;
		prevbegin = begin;
		prev = cur;
		action++;
		num++;
		count--;
	}
	if (begin >= pos) {
		return cur;
	}
	/* begin = pos on ref, cur = pos on query */
	while(count--) {
/* fprintf(stderr,"%c %d\n",*action,*num); */
		if (*action == 'M' || *action == '=' || *action == 'X') {
			begin += *num;
			cur += *num;
		} else if (*action == 'D' || *action == 'N') {
			begin += *num;
		} else if (*action == 'I') {
			cur += *num;
		} else if (*action == 'S' || *action == 'H') {
			/* only start clip is counted, not end: alignment stops there */
			return cur;
		}
		if (begin >= pos) break;
		prevbegin = begin;
		prev = cur;
		action++;
		num++;
	}
/*
if (action > cigar->action+cigar->size-1) {
fprintf(stderr,"bigger");
fflush(stderr);
} else {
fprintf(stderr,"%c %d cur=%d prev=%d pos=%d prevbegin=%d\n",*action,*num,cur,prev,pos,prevbegin);
fflush(stderr);
}
*/
	if (action > cigar->action+cigar->size-1 || *action == 'M' || *action == '=' || *action == 'X') {
		int temp = prev+(pos-prevbegin);
/* fprintf(stderr,"temp=%d\n",temp); */
		if (temp < cur) {
			cur = temp;
		} else if (end_ins && (action < cigar->action+cigar->size-1) && *(action+1) == 'I') {
			cur += *(num+1);
/* fprintf(stderr,"cur2=%d\n",cur); */
		}
/*
		if (action > cigar->action && *(action-1) == 'I') {
			cur += *(num-1);
		}
*/
		return cur;
	} else {
		return prev;
	}
}

int genomecomb_alignedseqObjCmd(ClientData clientData,	Tcl_Interp *interp, int argc, Tcl_Obj *CONST argv[])
{
	Cigar cigar;
	char *chrstring,*strand,*cigarstring,*seqstring,*reqchrstring;
	int chrsize,begin,end,reqbegin,reqend,seqsize,reqchrsize,addinsseq=0;
	int qhardclip,qhardclip2,qbegin,qend,seqstart,preN,postN,size;
/*
fprintf(stdout,"start");
fflush(stdout);
*/
	if (argc != 9 && argc != 10) {
		Tcl_WrongNumArgs(interp, 1, argv, "begin end cigar seq reqchr reqbegin reqend ?addinsseq?");
		return TCL_ERROR;
	}
	chrstring = Tcl_GetStringFromObj(argv[1],&chrsize);
	if (Tcl_GetIntFromObj(interp, argv[2], &begin) != TCL_OK) {
		return TCL_ERROR;
	}
	if (Tcl_GetIntFromObj(interp, argv[3], &end) != TCL_OK) {
		return TCL_ERROR;
	}
	cigarstring = Tcl_GetStringFromObj(argv[4],NULL);
	seqstring = Tcl_GetStringFromObj(argv[5],&seqsize);
	reqchrstring = Tcl_GetStringFromObj(argv[6],&reqchrsize);
	if (Tcl_GetIntFromObj(interp, argv[7], &reqbegin) != TCL_OK) {
		return TCL_ERROR;
	}
	if (Tcl_GetIntFromObj(interp, argv[8], &reqend) != TCL_OK) {
		return TCL_ERROR;
	}
	if (argc > 9) {
		if (Tcl_GetIntFromObj(interp, argv[9], &addinsseq) != TCL_OK) {
			return TCL_ERROR;
		}
	}
	if (chrsize > 3 && chrstring[0] == 'c' && chrstring[1] == 'h' && chrstring[2] == 'r') {
		chrstring = chrstring + 3; chrsize += -3;
	}
	if (reqchrsize > 3 && reqchrstring[0] == 'c' && reqchrstring[1] == 'h' && reqchrstring[2] == 'r') {
		reqchrstring = reqchrstring + 3; reqchrsize += -3;
	}

	if (reqend <= begin || reqbegin >= end || chrsize != reqchrsize || strncmp(chrstring,reqchrstring,chrsize) != 0) {
		Tcl_SetObjResult(interp,Tcl_NewObj());
		return TCL_OK;
	}
	cigar.size = 0;
	cigar.memsize = 0;
	parse_cigar(&cigar,cigarstring);
	if (cigar.action[0] == 'H') {
		qhardclip = cigar.num[0];
	} else {
		qhardclip = 0;
	}
	if (cigar.action[cigar.size-1] == 'H') {
		qhardclip2 = cigar.num[cigar.size-1];
	} else {
		qhardclip2 = 0;
	}
	if (addinsseq) {
		qbegin = cigar_ref2seq(&cigar,begin,reqbegin,0);
		qend = cigar_ref2seq(&cigar,begin,reqend,1);
	} else {
		qbegin = cigar_ref2seq(&cigar,begin,reqbegin,1);
		qend = cigar_ref2seq(&cigar,begin,reqend,0);
	}
	if (qbegin < qhardclip) {
		preN = qhardclip - qbegin;
		seqstart = 0;
	} else {
		preN = 0;
		seqstart = qbegin - qhardclip; 
	}
	reqend = qend - qhardclip;
	if (reqend > seqsize) {
		postN = reqend - seqsize;
	} else {
		postN = 0;
	}
/*
fprintf(stdout,"begin: %d end: %d ",begin,end);
fprintf(stdout,"qbegin: %d qend: %d qhardclip:%d ",qbegin,qend,qhardclip);
fprintf(stdout,"preN: %d,postN: %d ",preN,postN);
fprintf(stdout,"\n",preN,postN);
fflush(stdout);
*/
	if (preN != 0 || postN != 0) {
		Tcl_SetObjResult(interp,Tcl_NewStringObj("?",1));
		return TCL_OK;
	}
	size = min(reqend,seqsize) - seqstart;
	if (size == 0) {
		Tcl_SetObjResult(interp,Tcl_NewStringObj("-",1));
	} else {
		Tcl_Obj *result=Tcl_NewObj();
		Tcl_SetStringObj(result,seqstring + seqstart,size);
		Tcl_SetObjResult(interp,result);
	}
	return TCL_OK;
}
