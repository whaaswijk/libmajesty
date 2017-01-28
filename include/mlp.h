/* C header file

File name [mlp.h]

Package [Majority Logic Package]

Synopsis [EPFL, LSI lab., Majority Logic Package (MLP)]

Description [External functions and data structures of MLP]

Author [Luca Amaru]

Copyright [Copyright (c) EPFL, Integrated Systems Laboratory]

*/

#ifndef MLP_H
#define MLP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MAXCAR 128

typedef enum{
        AND,
		OR,
		PI // primary input
}PrimitiveOp;

typedef struct LogicCircuit LogicCircuit;

typedef struct LogicGate LogicGate;

typedef struct MIG MIG;

typedef struct MAJ3 MAJ3;

typedef struct entry_t entry_t;

typedef struct hashtable_t hashtable_t;

struct LogicCircuit {
	unsigned int Nin, Nout, Ngates;
	LogicGate **in, **out, **gates;
	unsigned int *outcompl;
	unsigned int *outone;
	char **innames, **outnames;
	//here links to alternative graph representations
};

struct LogicGate { /* each logic gate in a logic circuit is assumed to be single output */
	unsigned int Nin, FanOut;
	unsigned int label; /*label of the node*/
	LogicGate **in;
	PrimitiveOp type; /*AND, OR, PI*/
	unsigned int *cmpl; /* complemented 1 /regular 0 edges at the inputs */
};


struct MIG {
	unsigned int Nin, Nout, Nnodes;
	MAJ3 **in, **out, **nodes, *one;
	unsigned int *outcompl;
	char **innames, **outnames;
};


struct MAJ3 {
	MAJ3 *in1, *in2, *in3;
	MAJ3 *aux;//support structure - just linking purposes
	MAJ3 **outEdges;
	unsigned int flag;//for various purposes
	unsigned int level;
	unsigned int fanout;
	unsigned int label;
	unsigned int value;//0 logic zero, 1 logic one, 2,3,4,5... uninitialized 
	unsigned int compl1, compl2, compl3;//complemented edges
	unsigned int PI;//primary input
	unsigned int PO;//primary output
};

/*unique table for MIG nodes uniqueness*/

struct entry_t {
	unsigned int key; /*key obatined by the hashing function*/
	MAJ3 *node; /*pointer here to the decision diagram node*/
	entry_t *next;  /*collision list*/
};
 
struct hashtable_t {
	unsigned int size, elements, collisions;
	entry_t **table; //just one big table - consider levelizing it
	entry_t *one;		
};

/* reader functions */ 

extern unsigned char getnoneofch(FILE *file);
extern unsigned int decode(FILE *file);
extern unsigned int aigerheader(FILE *file, unsigned int *M, unsigned int *I, unsigned int *L, unsigned int *O, unsigned int *A);
extern LogicCircuit *readaiger(FILE *file);

/* writer functions */

extern unsigned int writeVerilogLC(FILE *file, LogicCircuit *net);
extern unsigned int writeVerilogMIG(FILE *file, MIG *net);

/* cloning functions */

extern MIG *cloneLC(LogicCircuit *net);

/* memory functions */ 

extern void freemig(MIG *net);
extern void freecircuit(LogicCircuit *net);

/* MLP utilities */

extern void computedepth(MIG *net);
extern void recursivedepth(MAJ3 *node);
extern void resetflag(MIG *net);
extern MAJ3 *checkandret(MIG *net, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2, MAJ3 *in3, unsigned int compl3, unsigned int *finalcompl);
extern MAJ3 **getcriticalpath(MIG *net, unsigned int *Nx);
extern void deletemignode(MIG *mig, MAJ3 *tobedel);
extern void deletefanout(MAJ3* source, MAJ3 *tobedel);
extern void addfanout(MAJ3* source, MAJ3 *tobeadd);
extern void recursive_polish(MIG *mig, MAJ3 *fanzero);
extern MAJ3 *new_node(MIG *net, MAJ3 *parent, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2, MAJ3 *in3, unsigned int compl3);
extern unsigned int integrity_fanout_check(MIG *net);
extern unsigned int setrecursiveflag(MAJ3 *node);
extern MAJ3 *duplicatesubtree(MIG *net, MAJ3 *node);
extern MAJ3 *recursiveduplicate(MIG *net,MAJ3 *node);
extern unsigned int checkduplicateredundancy(MIG *mig);
extern unsigned int eliminateduplicateredundancy(MIG *mig);
extern unsigned int checkallduplicateredundancy(MIG *mig);
extern MAJ3 *new_node_nop(MIG *net, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2, MAJ3 *in3, unsigned int compl3);
extern void resetaux(MIG *net);
extern unsigned int substitute_in_flag(MIG *mig, MAJ3 *old, MAJ3 *nodenew, unsigned int compl_sub);
extern unsigned int maxdepth(MIG *mig);
extern unsigned int realmaj(MIG *mig);
extern unsigned int noduplicates(MIG *mig);
extern unsigned int replacenode(MIG *mig, MAJ3 *node, MAJ3 *nodenew, unsigned int cmpl);
extern unsigned int setrecursiveflagfanout1(MAJ3 *node);
extern MAJ3 *largefanoutcriticalpath(MIG *net, MAJ3 *node);
extern unsigned int setcriticaldepthflag(MAJ3 *node);
extern MAJ3 *searchflagwithmaxfanout(MIG *net);
extern MAJ3 *searchhighestoutputflag(MIG *net, MAJ3 *node);
extern MAJ3 *searchcriticalpathdominator(MIG *net, MAJ3 *node);
extern unsigned int searchneighborpair(MIG *mig, MAJ3 **pair1, MAJ3 **pair2);
extern void get_order(MAJ3 *node, MAJ3 **last, unsigned int *compllast, MAJ3 **second, unsigned int *complsecond, MAJ3 **first, unsigned int *complfirst);
extern unsigned int searchneighborpair_wone(MIG *mig, MAJ3 **pair1, MAJ3 **pair2);
extern unsigned int searchrelevancepair(MIG *mig, MAJ3 **pair1, MAJ3 **pair2);
extern unsigned int twooverthree(MIG *mig, MAJ3 *node1, MAJ3 *node2, MAJ3 **common1, unsigned int *compl1, MAJ3 **common2, unsigned int *compl2, MAJ3 **uncommon1, unsigned int *complun1, MAJ3 **uncommon2, unsigned int *complun2);
extern void shufflemig(MIG *mig);

/* hashing functions */

extern hashtable_t *ht_create(unsigned int size, MIG *mig);
extern unsigned int mlppair(unsigned int i, unsigned int j);
extern unsigned int hash_function(unsigned int keyin1, unsigned int compl1, unsigned int keyin2, unsigned int compl2, unsigned int keyin3, unsigned int compl3, hashtable_t *hashtable);
extern entry_t *ht_newpair(MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2,  MAJ3 *in3, unsigned int compl3, unsigned int key);
extern unsigned int is_it_the_node(entry_t *check, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2,  MAJ3 *in3, unsigned int compl3);
extern MAJ3 *find_or_create_node(MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2,  MAJ3 *in3, unsigned int compl3, hashtable_t *hashtable);
extern unsigned int MIGhashing(MIG *mig);
extern MAJ3 *buildMIG(MAJ3 *var, hashtable_t *hashtable);
extern void  freehashtable(hashtable_t *hashtable);
extern void  freeentry(entry_t *tobefreed);

unsigned int selective_depth(MIG *mig);
unsigned int depth_reduction(MIG *net, MAJ3 *node);
unsigned int irredundant(MIG *mig);
unsigned int maj_strategy(MIG *net, MAJ3 *node);
unsigned int aggressive_depth(MIG *mig, double overhead);

#ifdef __cplusplus
}
#endif

#endif
