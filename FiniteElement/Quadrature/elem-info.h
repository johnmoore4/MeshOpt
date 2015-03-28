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

/* $Id: elem-info.h,v 1.1 2014/04/04 04:30:00 zlb Exp $ */

#ifndef PHG_ELEM_INFO_H

typedef void (*ADD_TETRA_FUNC)
			(int v0, int v1, int v2, int v3, int bound_type[]);

/* Note: FACE_INFO[0] = number of nodes, FACE_INFO[1:] = node list.
 *
 * !!! IMPORTANT: the vertices must be listed circularly which is required by
 * 		  update_neighbours()
 */
typedef short EDGE_INFO[2];	/* list of vertices */
typedef short FACE_INFO[5];	/* 0: # of vertices, 1..: list of vertices */
typedef void (*FUNC2TET)(int verts[], ADD_TETRA_FUNC func, int bound_type[]);

typedef struct {
    const char  *name;		/* name of the element type */
    EDGE_INFO   *edge_info;	/* list of edges */
    FACE_INFO   *face_info;	/* list of faces */
    short       nvert;		/* number of vertices */
    short	nedge;		/* number of edges */
    short       nface;		/* number of faces */
    FUNC2TET    func2tet;	/* function to covert to tet */
} ELEM_INFO;

extern int phgElemInfoCount;
extern ELEM_INFO phgElemInfo[];

#ifdef __cplusplus
extern "C" {
#endif

void phgPrism2Tetra(int verts[], ADD_TETRA_FUNC func, int bound_type[]);
void phgHexa2Tetra(int verts[], ADD_TETRA_FUNC func, int bound_type[]);
void phgPyramid2Tetra(int verts[], ADD_TETRA_FUNC func, int bound_type[]);

#ifdef __cplusplus
}
#endif

#define PHG_ELEM_INFO_H
#endif
