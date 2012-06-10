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

int DStringGetLine_tcl(Tcl_Interp *interp, DString *linePtr,	Tcl_Channel in,Tcl_Obj *buffer,int *keepbufferpos) {
	char *bufferstring,*bufferpos;
	char *cur = linePtr->string;
	int bufferlen;
	int c,newdata=0,read;
	ssize_t size = 0;
	NODPRINT("==================== getline ==================== \n");
	bufferstring = Tcl_GetStringFromObj(buffer,&bufferlen);
	bufferpos = bufferstring + *keepbufferpos;
	bufferlen -= *keepbufferpos;
	NODPRINT("bufferlen=%d bufferstring =%d bufferpos=%d keepbufferpos=%d\n",bufferlen,bufferstring,bufferpos,*keepbufferpos);
	while (1) {
		if (!bufferlen) {
			read=Tcl_ReadChars(in, buffer, BUFFERSIZE, 0);
			bufferstring = Tcl_GetByteArrayFromObj(buffer,&bufferlen);
			bufferpos = bufferstring;
			NODPRINT("\n==== read buffer (%d)\n%.BUFFERSIZEs\n",read, bufferpos);
			bufferlen = read;
			if (read == -1) break;
			if (!read && Tcl_Eof(in)) break;
		}
		newdata = 1;
		bufferlen--;
		c = *bufferpos++;
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
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,2*linePtr->memsize);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	*keepbufferpos = bufferpos-bufferstring;
	NODPRINT("-after: bufferlen=%d bufferstring =%d bufferpos=%d keepbufferpos=%d\n",bufferlen,bufferstring,bufferpos,*keepbufferpos);
	NODPRINT("\n==== result ====\n%s\n",linePtr->string);
	NODPRINT("==================== getline finished ==================== \n");
	return size;
}

int
genomecomb_tsv_select_ObjCmd(Tcl_Interp *interp, Tcl_Obj *listObj, Tcl_Obj *testObj,int *result)
{
/* readchars 2 */
	DString *line = NULL;
	Tcl_Channel in,out;
	Tcl_Obj *buffer;
	char *bufferstring;
	int read,error,keepbufferpos;
	in=Tcl_GetStdChannel(TCL_STDIN);
	out=Tcl_GetStdChannel(TCL_STDOUT);
	Tcl_SetChannelOption(interp,in,"-encoding","binary");
	Tcl_SetChannelOption(interp,out,"-encoding","binary");
	buffer = Tcl_NewByteArrayObj((CONST unsigned char *)"",0);
	Tcl_IncrRefCount(buffer);
//	Tcl_SetByteArrayLength(buffer,1001);
	bufferstring = Tcl_GetStringFromObj(buffer,NULL);
	keepbufferpos = 0;
DPRINT("start bufferstring =%d keepbufferpos=%d\n",bufferstring,keepbufferpos);
	line = DStringNew();
	while ((read = DStringGetLine_tcl(interp,line, in,buffer,&keepbufferpos)) != -1) {
		fprintf(stdout,"%s\n",line->string);
	}
	if (line) {DStringDestroy(line);}
	Tcl_DecrRefCount(buffer);
/* Dstring with buffer */
/*
	Buffer buffer;
	DString *line = NULL;
	ssize_t read;
	InitBuffer(&buffer,1000);
	line = DStringNew();
	while ((read = DStringGetLine_b(line, stdin,&buffer)) != -1) {
		if (line->string[0] == '\0') break;
		fprintf(stdout,"%s\n",line->string);
	}
	if (line) {DStringDestroy(line);}
	DelBuffer(&buffer);
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
