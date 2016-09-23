/* C file

File name [mlpmemory.c]

Package [Majority Logic Package]

Synopsis [memory management]

Description [memory management]

Author [Luca Amaru]

Copyright [Copyright (c) EPFL, Integrated Systems Laboratory]

*/

#include "mlp.h"


/* function freecircuit

Synopsis [freecircuit]

Description [freecircuit]

Side effects []

*/ 

void freecircuit(LogicCircuit *net) {
	unsigned int i;
	
	for(i=0;i<net->Nin;i++) {
		free(net->in[i]->cmpl);
		free(net->in[i]);
		free(net->innames[i]);
	}
	free(net->in);
	
	for(i=0;i<net->Ngates;i++) {
		free(net->gates[i]->cmpl);
		free(net->gates[i]->in);
		free(net->gates[i]);
	}
	free(net->gates);
	free(net->out);
	free(net->outcompl);
	free(net->outone);
	for(i=0;i<net->Nout;i++) {
		free(net->outnames[i]);
	}
	free(net->innames);
	free(net->outnames);
	
	free(net);
	
}/*end freecircuit*/

/* function freemig

Synopsis [freemig]

Description [freemig]

Side effects []

*/ 
void freemig(MIG *net) {
	unsigned int i;
	//unsigned int tot, totedges;
	
	//tot=net->Nnodes+net->Nin+1;
	//totedges=0;
	
	for(i=0;i<net->Nin;i++) {
		//totedges=totedges+net->in[i]->fanout;
		if(net->in[i]->fanout>0) {
			free(net->in[i]->outEdges);
		}
		free(net->in[i]);
		free(net->innames[i]);
	}
	free(net->in);
	
	for(i=0;i<net->Nnodes;i++) {
		//totedges=totedges+net->nodes[i]->fanout;
		if(net->nodes[i]->fanout>0) {
			free(net->nodes[i]->outEdges);
		}
		free(net->nodes[i]);
	}
	free(net->nodes);
	free(net->out);
	free(net->outcompl);
	for(i=0;i<net->Nout;i++) {
		free(net->outnames[i]);
	}
	free(net->innames);
	free(net->outnames);
	//totedges=totedges+net->one->fanout;
	// 
	if(net->one->fanout>0) {
		free(net->one->outEdges);
	}
	free(net->one);
	free(net);
	
	//printf("nodes + inputs + one: %u, edges: %u\n",tot,totedges);
	
}/*end freemig*/
