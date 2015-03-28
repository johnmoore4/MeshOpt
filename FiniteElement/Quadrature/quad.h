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

/* $Id: quad.h,v 1.5 2014/10/21 03:00:08 zlb Exp $ */

#ifndef PHG_QUAD_H

typedef double FLOAT;
typedef int INT;
typedef int SHORT;

typedef struct QUAD_ {
    char *name;			/* name of the quadrature formulae */
    int dim;			/* dimension, 1: edge, 2: face, 3: tetra */
    int order;			/* exact for polynomials of order 'order' */
    int npoints;		/* number of points */
    FLOAT *points;		/* barycentric coordinates of quad. points */
    FLOAT *weights;		/* weights */
    SHORT id;			/* id (for use with reference count) */
} QUAD;

#define QUAD_DEFAULT	-1

/* 1D quadrature rules */
extern QUAD QUAD_1D_P1_;
#define QUAD_1D_P1 (&QUAD_1D_P1_)
extern QUAD QUAD_1D_P3_;
#define QUAD_1D_P2 (&QUAD_1D_P3_)
#define QUAD_1D_P3 (&QUAD_1D_P3_)
extern QUAD QUAD_1D_P5_;
#define QUAD_1D_P4 (&QUAD_1D_P5_)
#define QUAD_1D_P5 (&QUAD_1D_P5_)
extern QUAD QUAD_1D_P7_;
#define QUAD_1D_P6 (&QUAD_1D_P7_)
#define QUAD_1D_P7 (&QUAD_1D_P7_)
extern QUAD QUAD_1D_P9_;
#define QUAD_1D_P8 (&QUAD_1D_P9_)
#define QUAD_1D_P9 (&QUAD_1D_P9_)
extern QUAD QUAD_1D_P11_;
#define QUAD_1D_P10 (&QUAD_1D_P11_)
#define QUAD_1D_P11 (&QUAD_1D_P11_)
extern QUAD QUAD_1D_P13_;
#define QUAD_1D_P12 (&QUAD_1D_P13_)
#define QUAD_1D_P13 (&QUAD_1D_P13_)
extern QUAD QUAD_1D_P15_;
#define QUAD_1D_P14 (&QUAD_1D_P15_)
#define QUAD_1D_P15 (&QUAD_1D_P15_)
extern QUAD QUAD_1D_P17_;
#define QUAD_1D_P16 (&QUAD_1D_P17_)
#define QUAD_1D_P17 (&QUAD_1D_P17_)
extern QUAD QUAD_1D_P19_;
#define QUAD_1D_P18 (&QUAD_1D_P19_)
#define QUAD_1D_P19 (&QUAD_1D_P19_)
extern QUAD QUAD_1D_P21_;
#define QUAD_1D_P20 (&QUAD_1D_P21_)
#define QUAD_1D_P21 (&QUAD_1D_P21_)

/* 2D quadrature rules */
extern QUAD QUAD_2D_P1_;
#define QUAD_2D_P1 (&QUAD_2D_P1_)
extern QUAD QUAD_2D_P2_;
#define QUAD_2D_P2 (&QUAD_2D_P2_)
extern QUAD QUAD_2D_P3_;
#define QUAD_2D_P3 (&QUAD_2D_P3_)
extern QUAD QUAD_2D_P4_;
#define QUAD_2D_P4 (&QUAD_2D_P4_)
extern QUAD QUAD_2D_P5_;
#define QUAD_2D_P5 (&QUAD_2D_P5_)
extern QUAD QUAD_2D_P6_;
#define QUAD_2D_P6 (&QUAD_2D_P6_)
extern QUAD QUAD_2D_P7_;
#define QUAD_2D_P7 (&QUAD_2D_P7_)
extern QUAD QUAD_2D_P8_;
#define QUAD_2D_P8 (&QUAD_2D_P8_)
extern QUAD QUAD_2D_P9_;
#define QUAD_2D_P9 (&QUAD_2D_P9_)
extern QUAD QUAD_2D_P10_;
#define QUAD_2D_P10 (&QUAD_2D_P10_)
extern QUAD QUAD_2D_P11_;
#define QUAD_2D_P11 (&QUAD_2D_P11_)
extern QUAD QUAD_2D_P12_;
#define QUAD_2D_P12 (&QUAD_2D_P12_)
extern QUAD QUAD_2D_P13_;
#define QUAD_2D_P13 (&QUAD_2D_P13_)
extern QUAD QUAD_2D_P14_;
#define QUAD_2D_P14 (&QUAD_2D_P14_)
extern QUAD QUAD_2D_P15_;
#define QUAD_2D_P15 (&QUAD_2D_P15_)
extern QUAD QUAD_2D_P16_;
#define QUAD_2D_P16 (&QUAD_2D_P16_)
extern QUAD QUAD_2D_P17_;
#define QUAD_2D_P17 (&QUAD_2D_P17_)
extern QUAD QUAD_2D_P18_;
#define QUAD_2D_P18 (&QUAD_2D_P18_)
extern QUAD QUAD_2D_P19_;
#define QUAD_2D_P19 (&QUAD_2D_P19_)
extern QUAD QUAD_2D_P20_;
#define QUAD_2D_P20 (&QUAD_2D_P20_)
extern QUAD QUAD_2D_P21_;
#define QUAD_2D_P21 (&QUAD_2D_P21_)
extern QUAD QUAD_2D_P22_;
#define QUAD_2D_P22 (&QUAD_2D_P22_)
extern QUAD QUAD_2D_P23_;
#define QUAD_2D_P23 (&QUAD_2D_P23_)
extern QUAD QUAD_2D_P24_;
#define QUAD_2D_P24 (&QUAD_2D_P24_)
extern QUAD QUAD_2D_P25_;
#define QUAD_2D_P25 (&QUAD_2D_P25_)
extern QUAD QUAD_2D_P26_;
#define QUAD_2D_P26 (&QUAD_2D_P26_)
extern QUAD QUAD_2D_P27_;
#define QUAD_2D_P27 (&QUAD_2D_P27_)
extern QUAD QUAD_2D_P28_;
#define QUAD_2D_P28 (&QUAD_2D_P28_)
extern QUAD QUAD_2D_P29_;
#define QUAD_2D_P29 (&QUAD_2D_P29_)

/* 3D quadrature rules */
extern QUAD QUAD_3D_P1_;
#define QUAD_3D_P1 (&QUAD_3D_P1_)
extern QUAD QUAD_3D_P2_;
#define QUAD_3D_P2 (&QUAD_3D_P2_)
extern QUAD QUAD_3D_P3_;
#define QUAD_3D_P3 (&QUAD_3D_P3_)
extern QUAD QUAD_3D_P4_;
#define QUAD_3D_P4 (&QUAD_3D_P4_)
extern QUAD QUAD_3D_P5_;
#define QUAD_3D_P5 (&QUAD_3D_P5_)
extern QUAD QUAD_3D_P6_;
#define QUAD_3D_P6 (&QUAD_3D_P6_)
extern QUAD QUAD_3D_P7_;
#define QUAD_3D_P7 (&QUAD_3D_P7_)
extern QUAD QUAD_3D_P8_;
#define QUAD_3D_P8 (&QUAD_3D_P8_)
extern QUAD QUAD_3D_P9_;
#define QUAD_3D_P9 (&QUAD_3D_P9_)
extern QUAD QUAD_3D_P10_;
#define QUAD_3D_P10 (&QUAD_3D_P10_)
extern QUAD QUAD_3D_P11_;
#define QUAD_3D_P11 (&QUAD_3D_P11_)
extern QUAD QUAD_3D_P12_;
#define QUAD_3D_P12 (&QUAD_3D_P12_)
extern QUAD QUAD_3D_P13_;
#define QUAD_3D_P13 (&QUAD_3D_P13_)
extern QUAD QUAD_3D_P14_;
#define QUAD_3D_P14 (&QUAD_3D_P14_)


#define PHG_QUAD_H
#endif
