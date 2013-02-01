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

/*
 * 
 * tsv_selectc_indexed queryproc varcolumns outcolumns indexfiles
 * version of tsv_selectc that used indexfiles (files containing one column) instead of stdin
 *
 * This is a rather dirty function, for speed reasons:
 * This bypasses Tcl channels completely
 * It uses non-threadsafe file access (getc_unlocked), and bytearrays.
 * All this makes it order of magnitudes faster than going the clean way through Tcl
 */

int 
genomecomb_tsv_select_indexed_ObjCmd (ClientData clientData,	Tcl_Interp *interp, int argc, Tcl_Obj *CONST argv[])
{
	/* Dstring */
	Tcl_CmdInfo cmdinfo,cmdinfo2;
	FILE **inds=NULL;
	Tcl_Obj **objv,*queryresult,*querycmd=NULL;
	int listobjc,listoutc,listindc,*cols=NULL,*outs=NULL;
	Tcl_Obj **listobjv,**listoutv,**listindv;
	DString *line = NULL,*array;
	int read;
	int maxtab=0,objc,show,i,error;
	unsigned int line_nr=0;
	if ((argc < 5)||(argc > 5)) {
		Tcl_WrongNumArgs(interp, 1, argv, "queryproc varcolumns outcolumns infiles");
		return TCL_ERROR;
	}
	if (Tcl_ListObjGetElements(interp, argv[2], &listobjc, &listobjv) != TCL_OK) {
		return TCL_ERROR;
	}
	if (Tcl_ListObjGetElements(interp, argv[3], &listoutc, &listoutv) != TCL_OK) {
		return TCL_ERROR;
	}
	if (Tcl_ListObjGetElements(interp, argv[4], &listindc, &listindv) != TCL_OK) {
		return TCL_ERROR;
	}
	cols = (int *)Tcl_Alloc(listobjc*sizeof(int));
	for (i = 0 ; i < listobjc ; i++) {
		error = Tcl_GetIntFromObj(interp,listobjv[i],cols+i);
		if (error) {return error;}
		if (cols[i] > maxtab) maxtab = cols[i];
	}
	outs = (int *)Tcl_Alloc(listoutc*sizeof(int));
	for (i = 0 ; i < listoutc ; i++) {
		error = Tcl_GetIntFromObj(interp,listoutv[i],outs+i);
		if (error) {
			outs[i] = -1;
		} else if (outs[i] > maxtab) {
			maxtab = outs[i];
		}
	}
	inds = (FILE **)Tcl_Alloc(listindc*sizeof(FILE *));
	for (i = 0 ; i < listindc ; i++) {
		inds[i] = fopen(Tcl_GetStringFromObj(listindv[i],NULL),"r");
		if (inds[i] == NULL) {
			Tcl_ResetResult(interp);
			Tcl_AppendResult(interp, "could not open file \"", Tcl_GetStringFromObj(listindv[i],NULL), "\"", (char *) NULL);
			return TCL_ERROR;
		}
	}
	line = DStringNew();
	array = DStringArrayNew(listindc);
	objv = (Tcl_Obj **)Tcl_Alloc((listobjc+2)*sizeof(Tcl_Obj *));
	querycmd = argv[1];
	objv[0]= querycmd;
	i = Tcl_GetCommandInfo(interp,Tcl_GetStringFromObj(querycmd,NULL),&cmdinfo);
	if (i==0) {return TCL_ERROR;}
	for (i = 1 ; i <= listobjc ; i++) {
		objv[i] = Tcl_NewByteArrayObj((unsigned char *)"",0);
		Tcl_IncrRefCount(objv[i]);
	}
	objv[i] = NULL;
	objc = listobjc+1;
	while (1) {
		for (i = 0 ; i < listindc ; i++) {
			read = DStringGetLine(array+i,inds[i]);
			if (read == -1) break;
		}
		if (read == -1) break;
		NODPRINT("line = %s",line->string);
		for (i = 0 ; i < listobjc ; i++) {
			if (Tcl_IsShared(objv[i+1])) {
				objv[i+1] = Tcl_NewObj();
				Tcl_IncrRefCount(objv[i+1]);
			}
			if (cols[i] != -1) {
				NODPRINT("col = %d, val=%.*s",cols[i],array[cols[i]].size,array[cols[i]].string);
				Tcl_SetByteArrayObj(objv[i+1],(unsigned char *)array[cols[i]].string,array[cols[i]].size);
			} else {
				Tcl_SetIntObj(objv[i+1],line_nr);
			}
		}
		error = cmdinfo.objProc(cmdinfo.objClientData,interp,objc,objv);
		if (error) {return error;}
		queryresult = Tcl_GetObjResult(interp);
		error = Tcl_GetBoolean(interp,Tcl_GetStringFromObj(queryresult,NULL),&show);
		if (error) {return error;}
		if (show) {
			register char *cur;
			int size,c;
			if (listoutc == 0) {
				cur = line->string; size = line->size;
				while (size--) {
					c = *cur++;
					if (c == '\0') {
						putc_unlocked('\t',stdout);
					} else {
						putc_unlocked(c,stdout);
					}
				}
				putc_unlocked('\n',stdout);
			} else {
				int out,first=1;
				for (i = 0 ; i < listoutc ; i++) {
					out=outs[i];
					if (!first) {putc_unlocked('\t',stdout);}
					first = 0;
					if (out != -1) {
						cur = array[out].string;
						size = array[out].size;
						while(size--) {putc_unlocked(*cur++,stdout);}
					} else {
						char *string;
						int len;
						string = Tcl_GetStringFromObj(listoutv[i],&len);
						if (len == 0) {
							fprintf(stdout,"%d",line_nr);
						} else {
							objv[0]=listoutv[i];
							error = Tcl_GetCommandInfo(interp,string,&cmdinfo2);
							if (error==0) {
								Tcl_ResetResult(interp);
								Tcl_AppendResult(interp, "unknown command: \"", string, "\"", (char *) NULL);
								return TCL_ERROR;
							}
							error = cmdinfo2.objProc(cmdinfo2.objClientData,interp,objc,objv);
							if (error) {return error;}
							queryresult = Tcl_GetObjResult(interp);
							cur = Tcl_GetStringFromObj(queryresult,&size);
							while(size--) {putc_unlocked(*cur++,stdout);}
						}
					}
				}
				putc_unlocked('\n',stdout);
			}
		}
		line_nr++;
	}
	for (i = 0 ; i < listindc ; i++) {
		fclose(inds[i]);
	}
	for (i = 1 ; i <= listobjc ; i++) {
		Tcl_DecrRefCount(objv[i]);
	}
	if (cols) {Tcl_Free((char *)cols);}
	if (outs) {Tcl_Free((char *)outs);}
	if (line) {DStringDestroy(line);}
	if (array) {DStringArrayDestroy(array);}
	return TCL_OK;
}
