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

/* $Id: quad-permu.h,v 1.2 2014/11/17 08:04:35 zlb Exp $ */

#ifndef PHG_QUAD_PERMU_H

# define _F(n)	n		/* no suffix */

/*-------------------------- Permutation macros ----------------------------*/

/* 0D */
#define Perm1(a)	_F(1.0)
#define Dup1(w)		_F(w)

#define Cons1(a)	Perm1(a)

/* unsymmetric, single point orbit */
#define Dup0		Dup1
#define Dup10		Dup0
#define Perm10		Perm1

/* 1D */
#define Perm2(a)	_F(0.5),_F(0.5)
#define Dup2(w)		_F(w)
#define Perm11(a)	_F(a),_F(1.)-(_F(a)), _F(1.)-(_F(a)),_F(a)
#define Dup11(w)	_F(w),_F(w)

#define Cons2(a)	Perm2(a)
#define Cons11(a)	Perm11(a)

/* unsymmetric, single point orbit */
#define Dup20		Dup0
#define Perm20(a)	_F(a),_F(1.)-(_F(a))

/* 2D */
#define Perm3(a)	_F(1.)/_F(3.),_F(1.)/_F(3.),_F(1.)/_F(3.)
#define Dup3(w)		_F(w)
#define Perm21(a)	_F(1.)-(_F(a))-(_F(a)),_F(a),_F(a), \
			_F(a),_F(1.)-(_F(a))-(_F(a)),_F(a), \
			_F(a),_F(a),_F(1.)-(_F(a))-(_F(a))
#define Dup21(w)	_F(w),_F(w),_F(w)
#define Perm111(a,b)	_F(a),_F(b),_F(1.)-(_F(a))-(_F(b)), \
			_F(a),_F(1.)-(_F(a))-(_F(b)),_F(b), \
			_F(b),_F(a),_F(1.)-(_F(a))-(_F(b)), \
			_F(b),_F(1.)-(_F(a))-(_F(b)),_F(a), \
			_F(1.)-(_F(a))-(_F(b)),_F(a),_F(b), \
			_F(1.)-(_F(a))-(_F(b)),_F(b),_F(a)
#define Dup111(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)

#define Cons3(a)	Perm3(a)
#define Cons21(a)	Perm21(a)
#define Cons111(a,b)	Perm111(a,b)

/* unsymmetric, single point orbit */
#define Dup30		Dup0
#define Perm30(a,b)	_F(a),_F(b),_F(1.)-(_F(a))-(_F(b))

/* 3D */
#define Perm4(a)	_F(0.25),_F(0.25),_F(0.25),_F(0.25)
#define Dup4(w)		_F(w)
#define Perm31(a)	_F(1.)-_F(3.)*(_F(a)),_F(a),_F(a),_F(a), \
			_F(a),_F(1.)-_F(3.)*(_F(a)),_F(a),_F(a), \
			_F(a),_F(a),_F(1.)-_F(3.)*(_F(a)),_F(a), \
			_F(a),_F(a),_F(a),_F(1.)-_F(3.)*(_F(a))
#define Dup31(w)	_F(w),_F(w),_F(w),_F(w)
#define Perm22(a)	_F(a),_F(a),_F(.5)-(_F(a)),_F(.5)-(_F(a)), \
			_F(a),_F(.5)-(_F(a)),_F(a),_F(.5)-(_F(a)), \
			_F(a),_F(.5)-(_F(a)),_F(.5)-(_F(a)),_F(a), \
			_F(.5)-(_F(a)),_F(a),_F(.5)-(_F(a)),_F(a), \
			_F(.5)-(_F(a)),_F(a),_F(a),_F(.5)-(_F(a)), \
			_F(.5)-(_F(a)),_F(.5)-(_F(a)),_F(a),_F(a)
#define Dup22(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)
#define Perm211(a,b)	_F(a),_F(a),_F(b),_F(1.)-(_F(a))-(_F(a))-(_F(b)), \
			_F(a),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(b), \
			_F(a),_F(b),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)), \
			_F(a),_F(b),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a), \
			_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(b), \
			_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(b),_F(a), \
			_F(b),_F(a),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)), \
			_F(b),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a), \
			_F(b),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(a), \
			_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(a),_F(b), \
			_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(b),_F(a), \
			_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(b),_F(a),_F(a)
#define Dup211(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w),\
			_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)
#define Perm0111(p,a,b,c) p,a,b,c, p,a,c,b, p,b,a,c, p,b,c,a, p,c,a,b, p,c,b,a
#define Perm1111(a,b,c) \
	Perm0111(_F(a),_F(b),_F(c),_F(1.)-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(b),_F(a),_F(c),_F(1.)-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(c),_F(a),_F(b),_F(1.)-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(1.)-(_F(a))-(_F(b))-(_F(c)),_F(a),_F(b),_F(c))
#define Dup1111(w)	Dup111(w), Dup111(w), Dup111(w), Dup111(w)

#define Cons4(a)	Perm4(a)
#define Cons31(a)	Perm31(a)
#define Cons22(a)	Perm22(a)
#define Cons211(a,b)	Perm211(a,b)
#define Cons1111(a,b,c)	Perm111(a,b,c)

/* unsymmetric, single point orbit */
#define Dup40		Dup0
#define Perm40(a,b,c)	_F(a),_F(b),_F(c),_F(1.)-(_F(a))-(_F(b))-(_F(c))

#define PHG_QUAD_PERMU_H
#endif
