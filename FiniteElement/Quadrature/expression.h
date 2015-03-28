/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: expression.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */

#ifndef PHG_EXPRESSION_H

#if !HAVE_LIBMATHEVAL

typedef struct {
    int op;		/* type of operator, variable, or number */
    void *arg1, *arg2;	/* pointers to operands */
} EXPR;

#else	/* !HAVE_LIBMATHEVAL */
typedef void *EXPR;
#endif	/* HAVE_LIBMATHEVAL */

#ifdef __cplusplus
extern "C" {
#endif

char *phgDump3DFunction(EXPR *tree);
EXPR *phgDefine3DFunction(const char *function_expression);
double phgEvaluate3DFunction(EXPR *tree, double x, double y, double z,
	double t);
EXPR *phgDup3DFunction(EXPR *tree);
void phgFree3DFunction(EXPR *tree);

#ifdef __cplusplus
}
#endif

#define PHG_EXPRESSION_H
#endif
