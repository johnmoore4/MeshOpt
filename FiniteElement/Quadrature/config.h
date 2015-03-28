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

/* # $Id: config.h.in,v 1.3 2014/04/22 05:14:36 zlb Exp $ */

#ifndef PHG_CONFIG_H
#define PHG_CONFIG_H

#define HAVE_FCLOSEALL	0
#define HAVE_REGEX_H	0
#define USE_ANONYMOUS_UNION 1

/* Note: USE_HB_FOR_DG==0 ==> use DOF_Pn for DG5-DG15,
 *	 USE_HB_FOR_DG==1 ==> use DOF_HBn for DG5-DG15 */
#define	USE_HB_FOR_DG 0

#define USE_MPI 1
#define USE_MPIIO 0
#define USE_OMP 0
#define USE_PAPI 0
#define USE_TCL 1
#define USE_TK 1
#define USE_METIS 1
#define USE_ZOLTAN 1
#define USE_PARMETIS 1
#define USE_VTK 1
#define DEBUG 0
#define ALLOW_CURVED_BOUNDARY 1
#define PHG_BIGENDIAN 1
#define GRID_TYPE_TET 0
#define GRID_TYPE_HEX 1
#define GRID_TYPE_MIX 2
#define PHG_GRID_TYPE GRID_TYPE_TET

/* linear solvers */

#define USE_HYPRE 0
#define HYPRE_VERSION_MAJOR 0
#define HYPRE_VERSION_MINOR 0
#define USE_PETSC 0
#define USE_TRILINOS 0
#define HAVE_TRILINOS_VERSION_H 0
#define USE_SUPERLU 0
#define USE_SUPERLU_PARSYMBFACT	0
#define USE_PASTIX 0
#define USE_HIPS 0
#define USE_MUMPS 0
#define USE_SPOOLES 0
#define USE_PARDISO 0
#define USE_SSPARSE 0
#define USE_LASPACK 0
#define USE_MINRES 0
#define USE_X9AMG   0
#define USE_OSKI 0

/* eigen solvers */

#define USE_ARPACK 0
#define USE_JDBSYM 0
#define USE_JDBSYM_PARALLEL 0
#define USE_BLOPEX 0
#define USE_SLEPC 0
#define USE_PRIMME 0
#define USE_GWLOBPCG 0
#define USE_TRILINOS_ANASAZI 0

/* Macros for Fortran names mangling defined by the configure script */
#define F77_FUNC(s,S)	s ## _		/* for names without underscore */
#define F77_FUNC_(s,S)	s ## __		/* for names containing underscores */
#define FC_FUNC(s,S)	s ## _		/* for names without underscore */
#define FC_FUNC_(s,S)	s ## __		/* for names containing underscores */

#define USE_BLAS 0
#define USE_LAPACK 0
#define BLAS_INT int
#define SIZEOF_BLAS_INT 0
#define LAPACK_INT BLAS_INT
#define SIZEOF_LAPACK_INT SIZEOF_BLAS_INT
#define USE_SCALAPACK 0

/* integer types */

#define dFMT "d"
#define IT_PHG IT_INT
#define SIZEOF_SHORT 0
#define SIZEOF_INT 0
#define SIZEOF_LONG 0
#define SIZEOF_LONG_LONG 0
#define TYPEOF_PHG_INT 0
#define SIZEOF_PHG_INT 0

/* float types */

#define FT_PHG FT_DOUBLE
#define SIZEOF_DOUBLE 0
#define SIZEOF_FLOAT 0
#define SIZEOF_LONG_DOUBLE 0
#define SIZEOF___FLOAT128 0

#define HAVE_LIBQUADMATH 0
#define HAVE_QUADMATH_H 0
#define HAVE_Q_SUFFIX 0
#define TYPEOF_PHG_FLOAT 0
#define SIZEOF_PHG_FLOAT 0

#define HAVE_FEENABLEEXCEPT 0
#define HAVE_LIBMATHEVAL 0
#define HAVE_LIBGMP 0
#define HAVE_LIBMPFR 0

#define HAVE_PSAPI 0

/* Note: it seems shenteng.sccas.cn doesn't like playing with creating
 * communicators. The following macro is used as a workaround.
 * When it's nonzero, only MPI_COMM_WORLD and MPI_COMM_SELF are used
 * (which means grids are distributed over either one or all processes) */
#define NO_NEW_COMMUNICATOR	0

/* Whether activate MPI_Sendrecv test (for detecting some GM communication
 * errors, only enabled by configure on LSSC2) */
#define SENDRECV_TEST 0

/* whether use native MPI_Alltoallv */
#define ALLTOALLV_HACK	0

/* whether enable workaround for realloc() bug (MVAPICH2 on DeepComp 7000) */
#define REALLOC_HACK 0

/* whether declare malloc()/free() calls as critical, needed for MVAPICH2
 * on LSSC3 (it seems MVAPICH2 1.4 makes malloc()/free() thread-unsafe) */
#define OMP_ALLOC_HACK 0

#if !defined(assert) && !defined(NDEBUG)
#if USE_MPI
#if DEBUG == 0 || defined(NDEBUG)
# define NDEBUG
# define assert(c)
#else	/* DEBUG == 0 || defined(NDEBUG) */
# define assert(c) \
    if (!(c)) phgError(1, "failed assertion at %s:%d, abort\n", \
						__FILE__, __LINE__)
#endif	/* DEBUG == 0 || defined(NDEBUG) */
#else	/* USE_MPI */
#include <assert.h>
#if DEBUG == 0
# define NDEBUG
#endif
#endif	/* USE_MPI */
#endif	/* !defined(assert) && !defined(NDEBUG) */

#define NO_DIST	1

#endif  /* !defined(PHG_CONFIG_H) */
