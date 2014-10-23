/*	
 *	 File:    genomecomb.c
 *	 Purpose: genomecomb extension to Tcl
 *	 Author:  Copyright (c) 1995 Peter De Rijk
 *
 *	 See the file "README" for information on usage and redistribution
 *	 of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <limits.h>
#include "debug.h"
#include "tcl.h"
#include "tools.h"
#include "genomecomb.h"

/*
 * 
 * tsv_selectc queryproc varcolumns outcolumns
 * this proc goes over stdin line by line interpreting it as a tab seperated file (no "", fields just cannot contain tabs)
 * for each line queryproc is called with the values in columns given by varcolumns (numbers) as parameters
 * if the result is true, output for this line is generated
 * outcolumns contains a list of columns (indicated by number) to output
 * if an element in the list is not an integer number, it will be used as a function that is called
 * with the same parameters as queryproc, and the result is put in the output column
 * If outcolumns is empty, output will consist of the entire line
 *
 * This is a rather dirty function, for speed reasons:
 * It just accesses stdin and stdout directly, bypassing Tcl channels completely
 * It uses non-threadsafe file access (getc_unlocked), and bytearrays.
 * All this makes it order of magnitudes faster than going the clean way through Tcl
 */

int 
genomecomb_tsv_select_ObjCmd (ClientData clientData,	Tcl_Interp *interp, int argc, Tcl_Obj *CONST argv[])
{
	/* Dstring */
	Tcl_CmdInfo cmdinfo,cmdinfo2;
	Tcl_Obj **objv,*queryresult,*querycmd=NULL;
	int listobjc,listoutc,*cols=NULL,*outs=NULL;
	Tcl_Obj **listobjv,**listoutv;
	DStringArray *array = NULL;
	DString *line = NULL;
	ssize_t read;
	int maxtab=0,objc,show,i,error,verbose = 0,sampling=0;
	unsigned int line_nr=0;
	uint64_t totalread = 0,next = 4294967296LL;
	if ((argc < 4)||(argc > 6)) {
		Tcl_WrongNumArgs(interp, 1, argv, "queryproc varcolumns outcolumns ?verbose? ?samplingskip?");
		return TCL_ERROR;
	}
	if (argc >= 5) {
		if (Tcl_GetIntFromObj(interp, argv[4], &verbose) != TCL_OK) {
			return TCL_ERROR;
		}
		if (verbose != 0) {
			next = abs(verbose);
		}
	}
	if (argc >= 6) {
		if (Tcl_GetIntFromObj(interp, argv[5], &sampling) != TCL_OK) {
			return TCL_ERROR;
		}
	}
	if (Tcl_ListObjGetElements(interp, argv[2], &listobjc, &listobjv) != TCL_OK) {
		return TCL_ERROR;
	}
	if (Tcl_ListObjGetElements(interp, argv[3], &listoutc, &listoutv) != TCL_OK) {
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
		} else if (outs[i] > maxtab) maxtab = outs[i];
	}
	line = DStringNew();
	array = DStringArrayNew(maxtab+2);
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
	totalread = 0;
	while ((read = DStringGetTab(line,stdin,maxtab,array,1,NULL)) != -1) {
		NODPRINT("%d,%d %s",maxtab, array->size, line->string);
		if (array->size <= maxtab) {
			fprintf(stderr, "wrong number of fields for line (%d without comments): %*.*s\n",line_nr,line->size,line->size,line->string);
			exit(1);
		}
		if (verbose < 0) {
			totalread += line->size+1;
			if (totalread >= next) {
				fprintf(stderr,"%ld\n",totalread);
				next += abs(verbose);
			}
		} else if (verbose != 0 && line_nr >= next) {
			fprintf(stderr,"%d\n",line_nr);
			next += abs(verbose);
		}
		NODPRINT("line = %s",line->string);
		for (i = 0 ; i < listobjc ; i++) {
			if (Tcl_IsShared(objv[i+1])) {
				objv[i+1] = Tcl_NewObj();
				Tcl_IncrRefCount(objv[i+1]);
			}
			if (cols[i] != -1) {
				NODPRINT("col = %d, val=%.*s",cols[i],array->data[cols[i]].size,array->data[cols[i]].string);
				Tcl_SetByteArrayObj(objv[i+1],(unsigned char *)array->data[cols[i]].string,array->data[cols[i]].size);
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
						DStringputs(array->data+out,stdout);
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
		if (sampling) {
			i = sampling;
			while(i--) {
				SkipLine(stdin);
			}
		}
	}
	for (i = 1 ; i <= listobjc ; i++) {
		Tcl_DecrRefCount(objv[i]);
	}
	if (cols) {Tcl_Free((char *)cols);}
	if (outs) {Tcl_Free((char *)outs);}
	if (line) {DStringDestroy(line);}
	if (array) {DStringArrayDestroy(array);}
	if (verbose < 0) {
		fprintf(stderr,"%ld\n",totalread);
	} else if (verbose > 0) {
		fprintf(stderr,"%d\n",line_nr);
	}
	return TCL_OK;
}

/*
 * 
 */

int 
genomecomb_loc_compare_ObjCmd (ClientData clientData,	Tcl_Interp *interp, int argc, Tcl_Obj *CONST argv[])
{
	char *loc1=NULL,*loc2;
	int loc1len,loc2len,result;
	if ((argc < 3)||(argc > 3)) {
		Tcl_WrongNumArgs(interp, 1, argv, "location1 location2");
		return TCL_ERROR;
	}
	loc1 = Tcl_GetStringFromObj(argv[1],&loc1len);
	loc2 = Tcl_GetStringFromObj(argv[2],&loc2len);
/* fprintf(stdout,"naturalcompare %s,%s,%d,%d\n",loc1,loc2,loc1len,loc2len); */
	result = loccompare(loc1,loc2,loc1len,loc2len);
	Tcl_SetObjResult(interp,Tcl_NewIntObj(result));
	return TCL_OK;
}

int 
genomecomb_nat_compare_ObjCmd (ClientData clientData,	Tcl_Interp *interp, int argc, Tcl_Obj *CONST argv[])
{
	char *loc1=NULL,*loc2;
	int loc1len,loc2len,result;
	if ((argc < 3)||(argc > 3)) {
		Tcl_WrongNumArgs(interp, 1, argv, "string1 string2");
		return TCL_ERROR;
	}
	loc1 = Tcl_GetStringFromObj(argv[1],&loc1len);
	loc2 = Tcl_GetStringFromObj(argv[2],&loc2len);
	result = naturalcompare(loc1,loc2,loc1len,loc2len);
	Tcl_SetObjResult(interp,Tcl_NewIntObj(result));
	return TCL_OK;
}
