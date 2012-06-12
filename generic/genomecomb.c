/*	
 *	 File:    genomecomb.c
 *	 Purpose: genomecomb extension to Tcl
 *	 Author:  Copyright (c) 1995 Peter De Rijk
 *
 *	 See the file "README" for information on usage and redistribution
 *	 of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "debug.h"
#include "tcl.h"
#include "tools.h"
#include "genomecomb.h"
#define EXACT		0
#define GLOB		1
#define REGEXP		2

#define BUFFERSIZE 1000

int DStringGetLine_tab(
	Tcl_Interp *interp, DString *linePtr,	Tcl_Channel in,
	Tcl_Obj *buffer,int *keepbufferpos,
	int maxtab,DString *result
) {
	char *bufferstring,*bufferpos;
	char *cur = linePtr->string;
	int bufferlen;
	int c,newdata=0,read,prevsize;
	int count=0;
	ssize_t size = 0;
	maxtab += 1;
	NODPRINT("==================== getline ==================== \n");
	bufferstring = Tcl_GetStringFromObj(buffer,&bufferlen);
	bufferpos = bufferstring + *keepbufferpos;
	bufferlen -= *keepbufferpos;
	NODPRINT("bufferlen=%d bufferstring =%d bufferpos=%d keepbufferpos=%d\n",bufferlen,bufferstring,bufferpos,*keepbufferpos);
	result[count].string = cur;
	result[count].memsize = -1;
	while (1) {
		if (!bufferlen) {
			read=Tcl_ReadChars(in, buffer, BUFFERSIZE, 0);
			bufferstring = Tcl_GetByteArrayFromObj(buffer,&bufferlen);
			bufferpos = bufferstring;
			NODPRINT("\n==== read buffer (%d)\n%.BUFFERSIZEs\n",read, bufferpos);
			bufferlen = read;
			if (read == -1) break;
			if (!read && Tcl_Eof(in)) {
				result[count].size = cur-result[count].string;
				break;
			}
		}
		newdata = 1;
		bufferlen--;
		c = *bufferpos++;
		if (count <= maxtab) {
			if  (c == '\t') {
				result[count].size = size;
//fprintf(stdout,"count=%d size=%d\n",count, result[count].size);
				count++;
			}
		}
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
		*keepbufferpos = 0;
		*cur = '\0';
		return -1;
	}
	/* fill rest of tab separated elements in array */
	if (count < maxtab) {
		result[count].size = size;
		while (count < maxtab) {
			result[++count].size = size;
		}
	}
	/* add \0 to line*/
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,2*linePtr->memsize);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	/* make array */
	count = 0;
	prevsize = 0;
	while (count < maxtab) {
		int pos = result[count].size;
		result[count].size = result[count].size-prevsize;
		result[count].string = linePtr->string + prevsize;
		result[count].memsize = -1;
//fprintf(stdout,"final count=%d size=%d pos=%d\n",count, result[count].size,result[count].string);
		prevsize = pos+1;
		count++;
	}
	*keepbufferpos = bufferpos-bufferstring;
	NODPRINT("-after: bufferlen=%d bufferstring =%d bufferpos=%d keepbufferpos=%d\n",bufferlen,bufferstring,bufferpos,*keepbufferpos);
	NODPRINT("\n==== result ====\n%s\n",linePtr->string);
	NODPRINT("==================== getline finished ==================== \n");
	return size;
}

int
genomecomb_tsv_select_ObjCmd(Tcl_Interp *interp, Tcl_Obj *listObj, Tcl_Obj *testObj,int *result)
{
/* Dstring */
	DString *line = NULL,*array;
	ssize_t read;
	int maxtab=2;
	line = DStringNew();
	array = DStringArrayNew(maxtab+2);
	while ((read = DStringGetTab(line,stdin,maxtab,array)) != -1) {
		if (line->string[0] == '\0') break;
		fprintf(stdout,"%.*s\n",array[maxtab].size,array[maxtab].string);
	}
	if (line) {DStringDestroy(line);}
	if (array) {DStringArrayDestroy(array);}

/* readchars 2 */
/*
	DString *line = NULL,*array;
	Tcl_Channel in,out;
	Tcl_Obj *buffer;
	char *bufferstring;
	int read,keepbufferpos,maxtab=2;
	in=Tcl_GetStdChannel(TCL_STDIN);
	out=Tcl_GetStdChannel(TCL_STDOUT);
	Tcl_SetChannelOption(interp,in,"-encoding","binary");
	Tcl_SetChannelOption(interp,out,"-encoding","binary");
	buffer = Tcl_NewByteArrayObj((CONST unsigned char *)"",0);
	Tcl_IncrRefCount(buffer);
	bufferstring = Tcl_GetStringFromObj(buffer,NULL);
	keepbufferpos = 0;
NODPRINT("start bufferstring =%d keepbufferpos=%d\n",bufferstring,keepbufferpos);
	line = DStringNew();
	array = DStringArrayNew(maxtab+1);
	while ((read = DStringGetLine_tab(interp,line, in,buffer,&keepbufferpos,maxtab,array)) != -1) {
		//fprintf(stdout,"-- %s\n",line->string);
		//fprintf(stdout,"%d %s\n",array[maxtab].size,array[maxtab].string);
		fprintf(stdout,"%.*s\n",array[maxtab].size,array[maxtab].string);
	}
	if (line) {DStringDestroy(line);}
	Tcl_DecrRefCount(buffer);
*/
/* readchars 1 */
/*
	Tcl_Channel in,out;
	int read,error;
	Tcl_Obj *buffer;
	buffer = Tcl_NewObj();
	in=Tcl_GetStdChannel(TCL_STDIN);
	out=Tcl_GetStdChannel(TCL_STDOUT);
	Tcl_SetChannelOption(interp,in,"-encoding","binary");
	Tcl_SetChannelOption(interp,out,"-encoding","binary");
	while (1) {
		read=Tcl_ReadChars(in, buffer, 1000, 0);
		if (read == -1) break;
		if (read == 0 && Tcl_Eof(in)) break;
		error=Tcl_WriteChars(out, (const char *)Tcl_GetByteArrayFromObj(buffer,NULL), read);
	}
	Tcl_DecrRefCount(buffer);
*/
	return TCL_OK;
}
