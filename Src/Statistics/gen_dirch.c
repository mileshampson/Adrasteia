#include "gen_dirch.h"
#include <assert.h>

/*
 *  gen_dirch.c	generating Dirichlet-distributed random vectors
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


/* gen_dirch.c
 * 7 Nov 2000		Kevin Karplus
 *
 *	Generator for random probability vectors with a Dirichlet
 * 	distribution. 
 */


/* The sample from a Dirichlet distribution is generated by a tree of
 * beta generators, each of which splits the probability passed down
 * from above.
 *
 * The leaves of the tree are the num_dim places in the output
 * probability vector (0..num_dim-1).
 *
 * The internal nodes of the tree are indicated by indexes >=num_dim.
 * The alpha value for an internal node is the sum of the alpha values
 * for the children.
 * 
 */



/* A heap is used to balance the tree to have as nearly equal alpha
 * values in siblings as possible (as in Huffman coding)
 */

/* Heapify the array heap (so that alpha[parent] < alpha[children]) */
static void heapify(int *heap, const double *alpha, int heap_size)
{
    int child, parent, i;	/* indexes of elements in heap */
    int moving;
    double alpha_moving, alpha_child;
    
    for (i=heap_size/2 - 1; i>=0; i--)
    {  	moving = heap[i];	/* move this down to make heap leagal */
        alpha_moving = alpha[moving];
	parent = i;
	for(;;)
	{   child = 2*parent+1;
	    if (child >= heap_size) break;
	    alpha_child = alpha[heap[child]];
	    if (child+1 < heap_size)
	    {   if (alpha[heap[child+1]] < alpha_child)
	    	{    alpha_child = alpha[heap[++child]];
		}
	    }
	    if (alpha_moving < alpha_child) break;
	    heap[parent] = heap[child];
	    parent=child;
	}
	heap[parent]=moving;
        /* i through heap_size-1 now satisfies heap property */
    }
}

static int pop_heap(int *heap, const double *alpha, int* heap_size)
{
    int child, parent, ret;	/* indexes of elements in heap */
    int moving;
    double alpha_moving, alpha_child;

    assert(*heap_size > 0 );
    ret = heap[0];
    moving = heap[--(*heap_size)];	/* trickle this down from root */
    alpha_moving = alpha[moving];
    parent=0;
    for(;;)
    {   child = 2*parent+1;
	if (child >= *heap_size) break;
	alpha_child = alpha[heap[child]];
	if (child+1 < *heap_size)
	{   if (alpha[heap[child+1]] < alpha_child)
	    {    alpha_child = alpha[heap[++child]];
	    }
	}
	if (alpha_moving < alpha_child) break;
	heap[parent] = heap[child];
	parent=child;
	/* i through *heap_size-1 now is satisfies heap property */
    }
    heap[parent]=moving;
    return ret;
}

static int push_heap(int *heap, const double *alpha, int* heap_size, int new_el)
{
    int child, parent;	/* indexes of elements in heap */
    double alpha_new;
    
    child=(*heap_size)++;
    alpha_new = alpha[new_el];
    while(child>0)
    {   parent = (child-1)/2;
    	if (alpha[heap[parent]] < alpha_new)	break;
	heap[child] = heap[parent];
	child = parent;
    }
    heap[child] = new_el;
}




/* store the alpha values and build a beta_tree for the generator */
void gen_dirch_initialize(gen_dirch_param *gen, int numd, const double *alp)
{
   int *heap;	int heap_size;
   int leaf, internal;	/* indexes of nodes in tree */
   double *alpha;	/* shorthand pointer for gen->alpha */

   gen->num_dim = numd;

   /* allocate space for the alphas---both for the leaves and the
    * internal nodes
    */
   gen->alpha = alpha = (double*)calloc(2*numd-1, sizeof(double));
   assert(gen->alpha);

   /* allocate space for the internal nodes */
   gen->betas = (gen_beta_param*) calloc(numd-1, sizeof(gen_beta_param));
   assert(gen->betas);
   gen->left = (int*) calloc(numd-1, sizeof(int));
   assert(gen->left);
   gen->right = (int*) calloc(numd-1, sizeof(int));
   assert(gen->right);

   /* Build a "Huffman encoding" tree using a heap.
    * First, fill the heap array with the indexes of the leaves
    */
    heap = (int*)calloc(numd, sizeof(int));
    for (leaf=numd-1; leaf>=0; leaf--)
    {  alpha[leaf] = alp[leaf];
       heap[leaf] = leaf;
    }
    heapify(heap, alpha, numd);
    
    heap_size=numd;
    internal = 0;	/* number of internal nodes created so far */
    while(heap_size>=2)
    {   gen->left[internal] = pop_heap(heap, alpha, &heap_size);
    	gen->right[internal] = pop_heap(heap, alpha, &heap_size);
	gen_beta_initialize(&(gen->betas[internal]), 
		alpha[gen->left[internal]], alpha[gen->right[internal]]);
	alpha[internal+numd] =	alpha[gen->left[internal]] 
				+ alpha[gen->right[internal]];
	push_heap(heap, alpha, &heap_size, internal+numd);
	internal++;
    }
    gen->root_index = internal+numd-1;
    assert(gen->root_index == 2*numd-2);
    free(heap);
}

/* Generate the probabilities for a subtree of the whole beta-tree
 * that defines the Dirichlet distribution, with the sum of the 
 * subtree being total_prob.
 */
static gen_dirch_subtree(const gen_dirch_param *gen, double *probs,
	double total_prob, int root)
{
    int internal;	/* index in left, right, and betas arrays */
    double x;	/* beta-distributed share for left subtree */
    internal = root-gen->num_dim;
    if (internal < 0)
    {    /* leaf node */
        probs[root] = (double)total_prob;
	return;
    }
    x = gen_beta(&(gen->betas[internal]));
    gen_dirch_subtree(gen, probs, total_prob*x, gen->left[internal]);
    gen_dirch_subtree(gen, probs, total_prob*(1-x), gen->right[internal]);
}

/* generate a probability vector according to the distribution and
 * store it in probs
 */
void gen_dirch(const gen_dirch_param *gen, double *probs)
{
    double sum=0.;
    int i;
    gen_dirch_subtree(gen, probs, 1.0, gen->root_index);
    /* rescale the probs to sum to 1, to reduce round-off error */
    for (i=gen->num_dim-1; i>=0; i--)
    {    sum += probs[i];
    }
    for (i=gen->num_dim-1; i>=0; i--)
    {	probs[i] /= sum;
    }    
}

/* free the storage internal to gen, but not gen itself */
void gen_dirch_free(gen_dirch_param *gen)
{   free(gen->right);
    free(gen->left);
    free(gen->betas);
    free(gen->alpha);
}
