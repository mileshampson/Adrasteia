/*
 *  gen_dirch.h	generating Dirichlet-distributed random vectors
 *  Copyright (C) 2000	Kevin Karplus
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; 
 *  version 2.1 of the License.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  See http://www.gnu.org/copyleft/lesser.html for the license details,
 *  or write to the Free Software Foundation, Inc., 59 Temple Place, Suite
 *  330, Boston, MA 02111-1307 USA
 */


/* This is a header file for a random number generator
 * that generates probability vectors according to a Dirichlet distribution.
 *
 * Kevin Karplus
 * 6 Nov 2000
 *
 * The functions here are made to look as much like a C++ class as
 * is possible in ANSI C.
 */


#ifndef GENDIRCHH
#define GENDIRCHH

#include "gen_beta.h"


/* The sample from a Dirichlet distribution is generated by a tree of
 * beta generators, each of which splits the probability passed down
 * from above.
 *
 * The leaves of the tree are the num_dim places in the output
 * probability vector (0..num_dim-1).
 *
 * The internal nodes of the tree are indicated by indexes >=num_dim.
 * 
 */


typedef struct _gen_dirch_param
    {   int num_dim;	/* number of dimensions of probability vector */
        double * alpha;	/* 2*num_dim-1 alpha values 
			 * First the leaves, then the internal nodes.
			 */
	int root_index;
	gen_beta_param *betas;	/* num_dim-1 internal nodes */
	int *left, *right;	/* indexes of subtrees for internal nodes*/
    } gen_dirch_param;



/* store the alpha values and build a beta_tree for the generator */
void gen_dirch_initialize(gen_dirch_param *gen, int numd, const double *alp);

/* generate a probability vector according to the distribution and
 * store it in probs
 */
void gen_dirch(const gen_dirch_param *gen, double *probs);


/* free the storage internal to gen, but not gen itself */
void gen_dirch_free(gen_dirch_param *gen);

#endif
