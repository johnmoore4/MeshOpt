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

/* $Id: phg.h,v 1.157 2014/12/09 05:10:44 zlb Exp $ */

#ifndef PHG_PHG_H

#define INT int
#define FLOAT double


#include "config.h"
#include <stdio.h>
#if USE_MPI
# include <mpi.h>
#endif	/* USE_MPI */
#if USE_OMP
# include <omp.h>
#endif	/* USE_OMP */
#include "expression.h"
#include "elem-info.h"

#ifndef __GNUC__
# define inline
#endif

#define PHG_VERSION		"0.9.2"		/**< Version string */
#define PHG_VERSION_MAJOR	0		/**< Major version number */
#define PHG_VERSION_MINOR	9		/**< Minor version number */
#define PHG_VERSION_SUBMINOR	2		/**< Subminor version number */

#ifndef Dim
# define Dim		3		/**< The Space dimensions */
#endif

/* Grid Type */
#if PHG_GRID_TYPE == GRID_TYPE_TET
/* Tetrahedra */
#define NVert		(Dim + 1)	/**< Number of vertices / element */
#define NFace		NVert		/**< Number of faces / element */
#if Dim == 3
# define NEdge		6		/**< Number of edges / element */
#  define NVertFace      3	     /**< Number of vertex on face */
#  define NEdgeFace      3	     /**< Number of edge on face */
#else
# define NEdge		NFace
#endif
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
/* Hexahedron */
# define NVert          8               /**< Number of vertices / element */
# define NFace          6		/**< Number of faces / element */
# define NVertFace      4	     /**< Number of vertex on face */
# define NEdgeFace      4	     /**< Number of edge on face */
# define NEdge         12	     /**< Number of edges / element */
# define NElem          1	     /**< Number of edges / element */
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
/* Mixed */
# define NVert          NVert2(e)         /**< Number of vertices / element */
# define NFace          NFace2(e)	  /**< Number of faces / element */
# define NVertFace      NVertFace2(e)	  /**< Number of vertex on face */
# define NEdgeFace      NEdgeFace2(e)	  /**< Number of edge on face */
# define NEdge          NEdge2(e)	  /**< Number of edges / element */
# define NElem          NElem2(e)	  /**< Number of edges / element */
#endif

/* elementary data types */

/* supported floating point types (FT_PHG is defined to be the one matching
 * FLOAT in config.h) */
#define FT_FLOAT	0
#define FT_DOUBLE	1
#define FT_LONG_DOUBLE	2
#define FT___FLOAT128	3

typedef TYPEOF_PHG_FLOAT FLOAT;	/* TYPEOF_PHG_FLOAT is defined in config.h */

/* supported integer types (IT_PHG is defined to be the one matching
 * INT in config.h) */
#define IT_SHORT	0
#define IT_INT		1
#define IT_LONG		2
#define IT_LONG_LONG	3

typedef TYPEOF_PHG_INT	INT;	/* TYPEOF_PHG_INT is defined in config.h */

typedef signed short	SHORT;		/**< short, at least 16 bits */
typedef signed char	CHAR;		/**< character */
typedef unsigned char	BYTE;		/**< byte */

/* MPI datatypesi and ops */

#if USE_MPI
# if FT_PHG == FT_LONG_DOUBLE
#  define PHG_MPI_FLOAT	MPI_LONG_DOUBLE
# elif FT_PHG == FT_DOUBLE
#  define PHG_MPI_FLOAT	MPI_DOUBLE
# elif FT_PHG == FT_FLOAT
#  define PHG_MPI_FLOAT	MPI_FLOAT
# elif FT_PHG == FT___FLOAT128
#  define PHG_MPI_FLOAT MPI_FLOAT128	/* user type defined in phg-mpi.c */
# else
#  error unexpected.
# endif

#if FT_PHG == FT___FLOAT128
# define PHG_SUM MPI_OP_SUM
# define PHG_MAX MPI_OP_MAX
#else
# define PHG_SUM MPI_SUM
# define PHG_MAX MPI_MAX
#endif

# if IT_PHG == IT_SHORT
#  define PHG_MPI_INT	MPI_SHORT
# elif IT_PHG == IT_INT
#  define PHG_MPI_INT	MPI_INT
# elif IT_PHG == IT_LONG
#  define PHG_MPI_INT	MPI_LONG
# elif IT_PHG == IT_LONG_LONG
#  define PHG_MPI_INT MPI_LONG_LONG
# else
#  error unexpected.
# endif
# define PHG_MPI_SHORT	MPI_SHORT
# define PHG_MPI_CHAR	MPI_CHAR
# define PHG_MPI_BYTE	MPI_BYTE
#endif	/* USE_MPI */

/* suffix for float constants */

#if FT_PHG == FT_FLOAT || FT_PHG == FT_DOUBLE
# define _F(n)	n		/* no suffix */
#elif FT_PHG == FT_LONG_DOUBLE
# define _F(n)	n ## L		/* 'L' suffix */
#elif FT_PHG == FT___FLOAT128
# if HAVE_Q_SUFFIX
#  define _F(n)	n ## Q		/* 'Q' suffix */
# else
#  define _F(n)	n ## L		/* fallback to 'L' suffix */
# endif
#else
# error unexpected!
#endif

/* mathematical functions for the FLOAT type */

#if HAVE_QUADMATH_H
#include <quadmath.h>
#else
# define FLT128_MAX 1.18973149535723176508575932662800702e4932Q
# define FLT128_MIN 3.36210314311209350626267781732175260e-4932Q
# define FLT128_EPSILON 1.92592994438723585305597794258492732e-34Q
#endif	/* !HAVE_QUADMATH_H */

#if !HAVE_LIBQUADMATH || !HAVE_QUADMATH_H
# define powq(a,b)	powl((long double)(a),(long double)(b))
# define sqrtq(a)	sqrtl((long double)(a))
# define fabsq(a)	fabsl((long double)(a))
# define logq(a)	logl((long double)(a))
# define expq(a)	expl((long double)(a))
# define sinq(a)	sinl((long double)(a))
# define asinq(a)	asinl((long double)(a))
# define cosq(a)	cosl((long double)(a))
# define acosq(a)	acosl((long double)(a))
# define tanq(a)	tanl((long double)(a))
# define atanq(a)	atanl((long double)(a))
# define floorq(a)	floorl((long double)(a))
# define ceilq(a)	ceill((long double)(a))
# define tgammaq(a)	tgammal((long double)(a))
#endif	/* !HAVE_LIBQUADMATH || !HAVE_QUADMATH_H */

#if FT_PHG == FT_LONG_DOUBLE
# define FLOAT_MAX	LDBL_MAX
# define FLOAT_MIN	LDBL_MIN
# define FLOAT_EPSILON	LDBL_EPSILON
# define Pow/*(a,b)*/	powl/*((FLOAT)(a),(FLOAT)(b))*/
# define Sqrt/*(a)*/	sqrtl/*((FLOAT)(a))*/
# define Fabs/*(a)*/	fabsl/*((FLOAT)(a))*/
# define Log/*(a)*/	logl/*((FLOAT)(a))*/
# define Exp/*(a)*/	expl/*((FLOAT)(a))*/
# define Sin/*(a)*/	sinl/*((FLOAT)(a))*/
# define Asin/*(a)*/	asinl/*((FLOAT)(a))*/
# define Cos/*(a)*/	cosl/*((FLOAT)(a))*/
# define Acos/*(a)*/	acosl/*((FLOAT)(a))*/
# define Tan/*(a)*/	tanl/*((FLOAT)(a))*/
# define Atan/*(a)*/	atanl/*((FLOAT)(a))*/
# define Floor/*(a)*/	floorl/*((FLOAT)(a))*/
# define Ceil/*(a)*/	ceill/*((FLOAT)(a))*/
# define Gamma/*(a)*/	tgammal/*((FLOAT)(a))*/
#elif FT_PHG == FT_DOUBLE
# define FLOAT_MAX	DBL_MAX
# define FLOAT_MIN	DBL_MIN
# define FLOAT_EPSILON	DBL_EPSILON
# define Pow/*(a,b)*/	pow/*((FLOAT)(a),(FLOAT)(b))*/
# define Sqrt/*(a)*/	sqrt/*((FLOAT)(a))*/
# define Fabs/*(a)*/	fabs/*((FLOAT)(a))*/
# define Log/*(a)*/	log/*((FLOAT)(a))*/
# define Exp/*(a)*/	exp/*((FLOAT)(a))*/
# define Sin/*(a)*/	sin/*((FLOAT)(a))*/
# define Asin/*(a)*/	asin/*((FLOAT)(a))*/
# define Cos/*(a)*/	cos/*((FLOAT)(a))*/
# define Acos/*(a)*/	acos/*((FLOAT)(a))*/
# define Tan/*(a)*/	tan/*((FLOAT)(a))*/
# define Atan/*(a)*/	atan/*((FLOAT)(a))*/
# define Floor/*(a)*/	floor/*((FLOAT)(a))*/
# define Ceil/*(a)*/	ceil/*((FLOAT)(a))*/
# define Gamma/*(a)*/	tgamma/*((FLOAT)(a))*/
#elif FT_PHG == FT_FLOAT
# define FLOAT_MAX	FLT_MAX
# define FLOAT_MIN	FLT_MIN
# define FLOAT_EPSILON	FLT_EPSILON
# define Pow/*(a,b)*/	powf/*((FLOAT)(a),(FLOAT)(b))*/
# define Sqrt/*(a)*/	sqrtf/*((FLOAT)(a))*/
# define Fabs/*(a)*/	fabsf/*((FLOAT)(a))*/
# define Log/*(a)*/	logf/*((FLOAT)(a))*/
# define Exp/*(a)*/	expf/*((FLOAT)(a))*/
# define Sin/*(a)*/	sinf/*((FLOAT)(a))*/
# define Asin/*(a)*/	asinf/*((FLOAT)(a))*/
# define Cos/*(a)*/	cosf/*((FLOAT)(a))*/
# define Acos/*(a)*/	acosf/*((FLOAT)(a))*/
# define Tan/*(a)*/	tanf/*((FLOAT)(a))*/
# define Atan/*(a)*/	atanf/*((FLOAT)(a))*/
# define Floor/*(a)*/	floorf/*((FLOAT)(a))*/
# define Ceil/*(a)*/	ceilf/*((FLOAT)(a))*/
# define Gamma/*(a)*/	tgammaf/*((FLOAT)(a))*/
#elif FT_PHG == FT___FLOAT128
# define FLOAT_MAX	FLT128_MAX
# define FLOAT_MIN	FLT128_MIN
# define FLOAT_EPSILON	FLT128_EPSILON
# define Pow/*(a,b)*/	powq/*((FLOAT)(a),(FLOAT)(b))*/
# define Sqrt/*(a)*/	sqrtq/*((FLOAT)(a))*/
# define Fabs/*(a)*/	fabsq/*((FLOAT)(a))*/
# define Log/*(a)*/	logq/*((FLOAT)(a))*/
# define Exp/*(a)*/	expq/*((FLOAT)(a))*/
# define Sin/*(a)*/	sinq/*((FLOAT)(a))*/
# define Asin/*(a)*/	asinq/*((FLOAT)(a))*/
# define Cos/*(a)*/	cosq/*((FLOAT)(a))*/
# define Acos/*(a)*/	acosq/*((FLOAT)(a))*/
# define Tan/*(a)*/	tanq/*((FLOAT)(a))*/
# define Atan/*(a)*/	atanq/*((FLOAT)(a))*/
# define Floor/*(a)*/	floorq/*((FLOAT)(a))*/
# define Ceil/*(a)*/	ceilq/*((FLOAT)(a))*/
# define Gamma/*(a)*/	tgammaq/*((FLOAT)(a))*/
#endif

#undef TRUE
#define TRUE	(1)
#undef FALSE
#define FALSE	(0)
typedef int BOOLEAN;			/**< boolean */

/* Geometric types. Also used for refinement types: EDGE means the refinement
   edge is an edge of the reference parallelepiped, FACE means the refinement
   edge is the diagonal of a face of the reference parallelepiped, and
   DIAGONAL means the refinement edge is the diagonal of the reference
   parallelepiped */
#define GTypeName(t) (					\
	((t) == VERTEX)		? "VERTEX" :		\
	((t) == EDGE)		? "EDGE" :		\
	((t) == FACE)		? "FACE" :		\
	((t) == DIAGONAL)	? "DIAGONAL" :		\
	((t) == OPPOSITE)	? "OPPOSITE" :		\
	((t) == MIXED)		? "MIXED" : "UNKNOWN")
enum {VERTEX, EDGE, FACE, DIAGONAL, OPPOSITE, MIXED};
typedef BYTE GTYPE;	/* force to BYTE to save memory */

#define VOLUME		DIAGONAL
#define P_VERTEX	OPPOSITE

/* types of faces (REMOTE means it has a remote neighbour), the type
   UNREFERENCED means the object does not belong to the leaf elements.

   Note:
   1. multiple bits may be set, e.g., DIRICHLET | REMOTE | INTERIOR
   2. REMOTE bit set ==> INTERIOR bit set
*/
#define BTypeName(t) (				\
	((t) & UNREFERENCED)	? "unused" :	\
	((t) & INTERIOR)	? (((t) & REMOTE) ? "remote" : "interior") : \
	((t) & DIRICHLET)	? "Dirichlet":	\
	((t) & NEUMANN)		? "Neumann" :	\
	((t) & BDRY_USER0)	? "user 0" :	\
	((t) & BDRY_USER1)	? "user 1" :	\
	((t) & BDRY_USER2)	? "user 2" :	\
	((t) & BDRY_USER3)	? "user 3" :	\
	((t) & BDRY_USER4)	? "user 4" :	\
	((t) & BDRY_USER5)	? "user 5" :	\
	((t) & BDRY_USER6)	? "user 6" :	\
	((t) & BDRY_USER7)	? "user 7" :	\
	((t) & BDRY_USER8)	? "user 8" :	\
	((t) & BDRY_USER9)	? "user 9" :	\
	((t) & UNDEFINED)	? "undefined" : "invalid")
enum {
    UNREFERENCED= 0,	/**< not referenced by leaf elements of the submesh */
    OWNER	= 1,	/**< assigned to current submesh */
    INTERIOR    = 2,	/**< interior of the mesh */
    REMOTE	= 4,	/**< shared with a remote process */
    DIRICHLET   = 8,	/**< belongs to a Dirichlet boundary */
    NEUMANN	= 16,	/**< belongs to a Neumann boundary */
    BDRY_USER0	= 32,	/**< belongs to a user type 0 boundary */
    BDRY_USER1	= 64,	/**< belongs to a user type 1 boundary */
    BDRY_USER2	= 128,	/**< belongs to a user type 2 boundary */
    BDRY_USER3	= 256,	/**< belongs to a user type 3 boundary */
    BDRY_USER4	= 512,	/**< belongs to a user type 4 boundary */
    BDRY_USER5	= 1024,	/**< belongs to a user type 5 boundary */
    BDRY_USER6	= 2048,	/**< belongs to a user type 6 boundary */
    BDRY_USER7	= 4096,	/**< belongs to a user type 7 boundary */
    BDRY_USER8	= 8192, /**< belongs to a user type 8 boundary */
    BDRY_USER9	= 16384,/**< belongs to a user type 9 boundary */
    UNDEFINED	= 32768	/**< belongs to a boundary but has no bdry type */
};
typedef unsigned short BTYPE;	/**< force BTYPE to ushort to save memory */

#define BDRY_SHIFT	3
#define BDRY_MASK	((BTYPE)(((-1) >> BDRY_SHIFT) << BDRY_SHIFT))
#define BDRY_UB		((BDRY_MASK >> BDRY_SHIFT) + 1)
#define IS_BDRY(b)	((b) & BDRY_MASK)

/* GRID flags */
enum {
    VERT_FLAG = 1,
    EDGE_FLAG = 2,
    FACE_FLAG = 4,
    ELEM_FLAG = 8,
    GEOM_FLAG = 16
};

typedef enum {
    GRID_TETRA	= 0,
    GRID_HEXA	= 1,
    GRID_MIXED	= 2
} GRID_ELEM_TYPE;

typedef enum {
    TETRA	= 0,
    PYRAMID	= 1,
    PRISM	= 2,
    HEXA	= 3
} ELEM_TYPE;


typedef FLOAT COORD[Dim];		/* space coordinates */

/**
 @brief The ELEMENT struct
 */
typedef struct ELEMENT_ {
    struct ELEMENT_	*children[2];		/**< Pointers to children */
    void		*neighbours[NFace];	/**< Pointers to neighbours.
						   It stores the
						   index in g->neighbours.list
						   if the neighbour is remote */
#if 0
    void		*userdata;		/**< Pointer to userdata, for
						   use by applications */
#endif
    void		*parent;		/* pointer to parent */
    INT			verts[NVert];		/**< List of vertices */
    INT			edges[NEdge];		/**< List of edges */
#if (Dim == 3)
    INT			faces[NFace];		/**< List of faces */
#endif
    INT			index;			/**< Element index */
    int			mark;			/**< Refi./coars./part. mark */
    int			region_mark;		/**< Regional mark */
#if ALLOW_CURVED_BOUNDARY
    /* Note: storing the indices vertex-wise instead of edge-wise saves
       memory, but it has problems when dealing with vertices across the
       interface of two surfaces */
    CHAR		bound_func[NEdge];	/**< Indices to curved boundary
						   edges (-1 if not curved) */
#endif
    BTYPE		bound_type[NFace];	/**< Boundary types */
    int 		ordering;		/**< Ordering of vertices:
						   - bits 0-2:   1st vertex
						   - bits 3-5:   2nd vertex
						   - bits 6-8:   3rd vertex
						   - bits 9-11:  4th vertex
						   - bits 12-14: 5th vertex
						   - bits 15-17: 6th vertex
						   - bits 18-20: 7th vertex
						   - bits 21-23: 8th vertex
						   */
    CHAR		generation;		/**< Generation (level) */
    /*union {*/
	BYTE		flag;			/**< Scratch flag */
	CHAR		hp_order;		/**< for phgHPSetup() */
    /*};*/
#if (Dim == 3)
    GTYPE		type;			/**< Refinement type */
#endif
    ELEM_TYPE           elem_type;              /**< Element type */
} ELEMENT;

/* for backward compatibility */
#if PHG_GRID_TYPE == GRID_TYPE_TET
typedef ELEMENT SIMPLEX;
#endif	/* PHG_GRID_TYPE == GRID_TYPE_TET */

#if USE_MPI
# define GlobalVertex(g,no) \
	((g) == NULL || (g)->L2Gmap_vert==NULL ? (no):(g)->L2Gmap_vert[no])
# define GlobalEdge(g,no) \
	((g) == NULL || (g)->L2Gmap_edge==NULL ? (no):(g)->L2Gmap_edge[no])
# define GlobalFace(g,no) \
	((g) == NULL || (g)->L2Gmap_face==NULL ? (no):(g)->L2Gmap_face[no])
# define GlobalElement(g,no) \
	((g) == NULL || (g)->L2Gmap_elem==NULL ? (no):(g)->L2Gmap_elem[no])
#else
# define GlobalVertex(g,no) (no)
# define GlobalEdge(g,no) (no)
# define GlobalFace(g,no) (no)
# define GlobalElement(g,no) (no)
#endif

#if USE_MPI
/* information about a remote neighbour of an element */
typedef struct {
    ELEMENT	*remote;	/* addr of the neighbour on remote process */
    ELEMENT	*local;		/* local element */
#if USE_ANONYMOUS_UNION
    union {
#endif	/* USE_ANONYMOUS_UNION */
	unsigned int path;	/* bitwise flags recording which child to
				   select when refining a boundary edge.
				   0=child 0, 1=child 1 (used by refine.c) */
	INT	peer_face;	/* local index of the shared face in the peer
				   process (set by phgUpdateBoundaryTypes and
				   used by build_L2Vmap) */
#if USE_ANONYMOUS_UNION
    };
#endif	/* USE_ANONYMOUS_UNION */
    int		rank;		/* remote process rank */
    BYTE	rface[NVertFace];/* the matching face of remote neighbour */
    BYTE	vertex;		/* the face shared with remote process */
    BYTE	op_vertex;	/* opposite vertex on remote process */
    BYTE	depth;		/* refinement count (depth) of the element */
} RNEIGHBOUR;
#endif

/**
 @brief The GRID struct
 */
#define nleaf_global	nelem_global
typedef struct GRID_ {
    FLOAT	lif;		/**< Load imbalance factor */
    FLOAT	bbox[2][Dim];	/**< BoundingBox of the mesh */
    FLOAT	volume;		/**< Volume of the mesh */
    char	*filename;	/**< name of the mesh file */
    COORD	*verts;		/**< List of vertices */
    ELEMENT	*roots;		/**< Root elements */
    ELEMENT	**elems;	/**< Array of size g->nlem, pointers to leaf
				     elements (for mesh traversal), NULL if
				     an element index is unreferenced.
				     Note for any i, i>=0 && i<nelems,
				     the following assertion holds:
					 assert(elems[i] == NULL ||
						elems[i]->index == i) */
    GRID_ELEM_TYPE  elem_type;  /**< Grid element types */  
    struct MAPPING_	*mapping;/**< Mapping from reference coord to
				      real world coord */
    struct DOF_	*geom;		/**< DOF containing the geometric data,
				     maintained by geom.c */
    struct DOF_	**dof;		/**< NULL terminated list of DOFs
				     defined on the mesh */
    struct HP_TYPE_ **hp;	/**< NULL terminated list of HP_TYPEs
				     defined on the mesh */

#if ALLOW_CURVED_BOUNDARY
    EXPR	**bdry_funcs;	/**< NULL terminated list of projection
				   functions for curved boundaries */
#endif
    BTYPE	*types_vert;	/**< Types of vertices (bit flags) */
    BTYPE	*types_edge;	/**< Types of edges (bit flags) */
    BTYPE	*types_face;	/**< Types of faces (bit flags) */
    BTYPE	*types_elem;	/**< Types of elements (bit flags) */

    /* inter-grid links */
    struct GRID_ *alien;	/**< Pointer to another grid */
    INT		*alien_map;	/**< Map of elements between the two grids */

    /* struct for handling periodic BC */
    struct PERIOD_ *period;

    /* Arrays containing ranks of, and local indices in the owner process
     * for vertices. They are non nil only when needed and are updated by
     * phgUpdateBoundaryTypes in utils.c and used by build_L2Vmap() in map.c.
     *
     * Note: an entry has value -1 iff it's uniquely owned. */
    INT		*owner_index_vert;
    int		*owner_rank_vert;

#if USE_MPI
    struct {
	/* elements with neighbours on other processes are stored here */
	int		*counts;	/* Counts of remote faces */
	int		*displs;	/* Start positions in list */
	RNEIGHBOUR	*list;		/* List of remote faces */
	INT		count;		/* Total number of remote faces */
	INT		allocated;	/* Allocated size of list */
    }		neighbours;	/**< List of remote faces */
    /* Note: don't access L2Gmaps directly, use the GlobalXxxx macros instead */
    INT		*L2Gmap_vert;	/**< Local to global map of vertices */
    INT		*L2Gmap_edge;	/**< Local to global map of edges */
    INT		*L2Gmap_face;	/**< Local to global map of faces */
    INT		*L2Gmap_elem;	/**< Local to global map of element indices */

    /* Arrays containing ranks of, and local indices in the owner process
     * for edges, faces and elements. They are non nil only when
     * needed and are updated by phgUpdateBoundaryTypes in utils.c
     * and used by build_L2Vmap() in map.c.
     *
     * Note: an entry has value -1 iff it's uniquely owned. */
    INT		*owner_index_edge;
    int		*owner_rank_edge;

    /* Note: owner_xxxx_face and owner_xxxx_elem are currently not used */
    INT		*owner_index_face;
    int		*owner_rank_face;

    INT		*owner_index_elem;
    int		*owner_rank_elem;

    MPI_Comm	comm;		/**< The mesh's communicator */
#else	/* USE_MPI */
    int		comm;
#endif	/* USE_MPI */
    INT		nleaf;		/**< Number of leaf elements in the submesh */
    INT		nvert;		/**< Number of vertices in the submesh */
    INT		nedge;		/**< Number of edges in the submesh */
#if Dim == 3
    INT		nface;		/**< Number of faces in the submesh */
#endif
    INT		nelem;		/**< Number of element indices in the submesh */

    /* counts of owned vertices/edges/faces/elements, they are updated by
     * the function phgUpdateBoundaryTypes */
    INT		nvert_owned;
    INT		nedge_owned;
    INT		nface_owned;
    INT		nelem_owned;

    /* global counters. 2D, V:E:T \approx 1:3:2. 3D, V:E:F:T \approx 1:6:10:5 */
    INT		nvert_global;	/**< Number of vertices in the global mesh */
    INT		nedge_global;	/**< Number of edges in the global mesh */
#if Dim == 3
    INT		nface_global;	/**< Number of faces in the global mesh */
#endif
    INT		nelem_global;	/**< Number of tetrahedra in the global mesh */

    INT		nroot;		/**< Number of root elements */
    INT		ntree;		/**< Total number of elements in the tree */

    int		serial_no;	/**< Number of grid changes */
    int		last_partitioner; /**< Last partitioner used */

    int		rank;		/**< Process rank (submesh no) */
    int		nprocs;		/**< Number of processes for the mesh */

    int		bc_alloc;	/**< Allocated size of bc_index and bc_rmap */
    int		bc_n;		/**< Actual size of bc_index and bc_rmap */
    int		*bc_list;	/**< sorted list of BDRY types */
    int		*bc_rmap;	/**< PHG:bc_list[i] <==> Medit:bc_rmap[i] */

    BYTE	flags;		/**< {VERT,EDGE,FACE,ELEM,GEOM}_FLAG bits */
} GRID;

#include "phg/utils.h"
#include "phg/refine.h"
#include "phg/coarsen.h"
#include "phg/phg-mpi.h"
#include "phg/grid.h"
#include "phg/albert.h"
#if USE_MPI
#if USE_PETSC && defined(NEED_PETSCKSP_H)
/* Note: "petscksp.h" must be included before "mpi-utils.h" since the former
 * defines some MPI functions as macros and the latter will redefine them */
# include <petscksp.h>
#endif	/* USE_PETSC */
# include "phg/mpi-utils.h"
# include "phg/mpi-2level.h"
#endif
#include "phg/medit.h"
#include "phg/gambit.h"
#include "phg/dof.h"
#include "phg/mapping.h"

#include "phg/distribute.h"

#include "phg/lagrange.h"
#include "phg/mass-lumping.h"
#include "phg/nedelec.h"
#include "phg/hfeb.h"
#include "phg/hierarchical-basis.h"
#include "phg/dg.h"

#include "phg/quad.h"
#include "phg/geom.h"

#include "phg/matvec.h"
#include "phg/map.h"
#include "phg/solver.h"

#include "phg/vtk.h"
#include "phg/opendx.h"
#include "phg/ensight.h"

#include "phg/option.h"
#include "phg/perf.h"
#include "phg/eigen.h"
#include "phg/dof-utils.h"
#include "phg/mark.h"

#include "phg/surface-cut.h"
#include "phg/intergrid.h"

#include "phg/hp.h"
#include "phg/periodic.h"

#include "phg/checkpoint.h"

#define PHG_PHG_H
#endif
