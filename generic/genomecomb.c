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

int 
genomecomb_tsv_select_ObjCmd (ClientData clientData,	Tcl_Interp *interp, int argc, Tcl_Obj *CONST argv[])
{
	/* Dstring */
	Tcl_CmdInfo cmdinfo;
	Tcl_Obj **objv,*queryresult;
	int listobjc,listoutc,listresc,*cols=NULL,*outs=NULL;
	Tcl_Obj **listobjv,**listoutv,**listresv;
	DString *line = NULL,*array;
	ssize_t read;
	int maxtab=0,objc,show,i,error;
	if ((argc < 3)||(argc > 3)) {
		Tcl_WrongNumArgs(interp, 1, argv, "varcolumns outcolumns");
		return TCL_ERROR;
	}
	if (Tcl_ListObjGetElements(interp, argv[1], &listobjc, &listobjv) != TCL_OK) {
		return TCL_ERROR;
	}
	if (Tcl_ListObjGetElements(interp, argv[2], &listoutc, &listoutv) != TCL_OK) {
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
		if (error) {return error;}
		if (outs[i] > maxtab) maxtab = outs[i];
	}
	line = DStringNew();
	array = DStringArrayNew(maxtab+2);
	objv = (Tcl_Obj **)Tcl_Alloc((listobjc+2)*sizeof(Tcl_Obj *));
	objv[0]=Tcl_NewStringObj("tsv_selectc_query",17);
	Tcl_IncrRefCount(objv[0]);
	i = Tcl_GetCommandInfo(interp,Tcl_GetStringFromObj(objv[0],NULL),&cmdinfo);
	if (i==0) {return TCL_ERROR;}
	for (i = 1 ; i <= listobjc ; i++) {
		objv[i] = Tcl_NewByteArrayObj((unsigned char *)"",0);
		Tcl_IncrRefCount(objv[i]);
	}
	objv[i] = NULL;
	objc = listobjc+1;
	while ((read = DStringGetTab(line,stdin,maxtab,array)) != -1) {
		if (line->string[0] == '\0') break;
		NODPRINT("line = %s",line->string);
		error = cmdinfo.objProc(cmdinfo.objClientData,interp,objc,objv);
		for (i = 0 ; i < listobjc ; i++) {
			NODPRINT("col = %d, val=%.*s",cols[i],array[cols[i]].size,array[cols[i]].string);
			Tcl_SetByteArrayObj(objv[i+1],(unsigned char *)array[cols[i]].string,array[cols[i]].size);
		}
		error = cmdinfo.objProc(cmdinfo.objClientData,interp,objc,objv);
		if (error) {return error;}
		queryresult = Tcl_GetObjResult(interp);
		if (Tcl_ListObjGetElements(interp, queryresult, &listresc, &listresv) != TCL_OK) {
			return TCL_ERROR;
		}
		if (!listresc) {
			show = 0;
		} else {
			error = Tcl_GetBoolean(interp,Tcl_GetStringFromObj(listresv[0],NULL),&show);
			if (error) {return error;}
		}
		if (show) {
			if (listoutc == 0) {
				fprintf(stdout,"%.*s\n",line->size,line->string);
			} else {
				char *sep="";
				int out,pos=1;
				for (i = 0 ; i < listoutc ; i++) {
					out=outs[i];
					if (out != -1) {
						fprintf(stdout,"%s%.*s",sep,array[out].size,array[out].string);
					} else if (pos < listresc) {
						fprintf(stdout,"%s%s",sep,Tcl_GetStringFromObj(listresv[pos++],NULL));
					}
					sep="\t";
				}
				fprintf(stdout,"\n");
			}
		}
	}
	Tcl_DecrRefCount(objv[0]);
	for (i = 1 ; i <= listobjc ; i++) {
		Tcl_DecrRefCount(objv[i]);
	}
	if (cols) {Tcl_Free((char *)cols);}
	if (outs) {Tcl_Free((char *)outs);}
	if (line) {DStringDestroy(line);}
	if (array) {DStringArrayDestroy(array);}
	return TCL_OK;
}
