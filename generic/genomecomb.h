/*	
 *	 File:    genomecomb.h
 *	 Purpose: genomecomb extension to Tcl
 *	 Author:  Copyright (c) 1995 Peter De Rijk
 *
 *	 See the file "README" for information on usage and redistribution
 *	 of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

/*
 * Windows needs to know which symbols to export.  Unix does not.
 * BUILD_Class should be undefined for Unix.
 */

#ifdef BUILD_Extral
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS DLLEXPORT
#endif /* BUILD_Extral */

EXTERN int Extral_Init(Tcl_Interp *interp);

