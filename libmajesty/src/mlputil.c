/* C file

File name [mlputil.c]

Package [Majority Logic Package]

Synopsis [utilities to operate on a MIG]

Description [utilities to operate on a MIG]

Author [Luca Amaru]

Copyright [Copyright (c) EPFL, Integrated Systems Laboratory]

*/

#include "mlp.h"

/* function computedepth

Synopsis [computedepth]

Description [computedepth in an MIG]

Side effects []

*/ 
void computedepth(MIG *net) {
	unsigned int i;

	resetflag(net);

	for(i=0;i<net->Nout;i++) {
		recursivedepth(net->out[i]);
		//printf("out %u level %u\n",i,net->out[i]->level);
	}
	
	resetflag(net);

}/*end of computedepth*/

/* function recursivedepth

Synopsis [recursivedepth]

Description [recursivedepth computation in the MIG - from outputs to inputs]

Side effects []

*/ 
void recursivedepth(MAJ3 *node) {
	if(node->PI==1) {
		node->level=0;//PI are at 0 level
		node->flag=1;
	}
	else if(node->flag==0) {
		recursivedepth(node->in1);
		recursivedepth(node->in2);
		recursivedepth(node->in3);
		if(node->in1->level>node->in2->level) {
			if(node->in3->level>node->in1->level) {
				node->level=node->in3->level+1;
			}
			else {
				node->level=node->in1->level+1;
			}
		}
		else {
			if(node->in3->level>node->in2->level) {
				node->level=node->in3->level+1;
			}
			else {
				node->level=node->in2->level+1;
			}
		}
		node->flag=1;
	}

}/*end of recursivedepth*/

/* function resetlabel

Synopsis [resetlabel]

Description [resetlabels in the MIG data structure]

Side effects []

*/ 
void resetflag(MIG *net) {
	unsigned int i;
	
	for(i=0;i<net->Nin;i++) {
		net->in[i]->flag=0;
	}
	
	for(i=0;i<net->Nnodes;i++) {
		net->nodes[i]->flag=0;
	}
	net->one->flag = 0;
}/*end of resetlabel*/

/* function resetaux

Synopsis [resetaux]

Description [reset aux attribute in the MIG data structure]

Side effects []

*/ 
void resetaux(MIG *net) {
	unsigned int i;
	
	for(i=0;i<net->Nin;i++) {
		net->in[i]->aux=NULL;
	}
	
	for(i=0;i<net->Nnodes;i++) {
		net->nodes[i]->aux=NULL;
	}
	net->one->aux = NULL;
}/*end of resetaux*/

/* function checkandret

Synopsis [checkandret]

Description [check and return a node having given children and compl edges - return NULL if not existing]

Side effects []

*/ 
MAJ3 *checkandret(MIG *net, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2, MAJ3 *in3, unsigned int compl3, unsigned int *finalcompl) {
	unsigned int i, j, k, exit_condition;
	MAJ3 *ret, *aux, *a, *b, *c;
	unsigned int compla, complb, complc;
	ret=NULL;//initial condition
	exit_condition=0;//starting point
	(*finalcompl)=0;
	
	i=j=k=0;//initial condition
	while(exit_condition==0) {
		if(i<in1->fanout) {
			j=0;
			while((j<in2->fanout)&&(exit_condition==0)) {
				if(in1->outEdges[i]==in2->outEdges[j]) {
					k=0;
					while((k<in3->fanout)&&(exit_condition==0)) {
						if(in3->outEdges[k]==in2->outEdges[j]) {
							aux=in3->outEdges[k];
							a=in1;
							compla=compl1;
							b=in2;
							complb=compl2;
							c=in3;
							complc=compl3;		
							
							if(((aux->in1==a)&&(aux->compl1==compla))&&((aux->in2==b)&&(aux->compl2==complb))&&((aux->in3==c)&&(aux->compl3==complc))) {
								ret=aux;
								exit_condition=1;
							}
							else if(((aux->in1==a)&&(aux->compl1==compla))&&((aux->in2==c)&&(aux->compl2==complc))&&((aux->in3==b)&&(aux->compl3==complb))) {
								ret=aux;
								exit_condition=1;
							}
							else if(((aux->in1==b)&&(aux->compl1==complb))&&((aux->in2==c)&&(aux->compl2==complc))&&((aux->in3==a)&&(aux->compl3==compla))) {
								ret=aux;
								exit_condition=1;
							}
							else if(((aux->in1==b)&&(aux->compl1==complb))&&((aux->in2==a)&&(aux->compl2==compla))&&((aux->in3==c)&&(aux->compl3==complc))) {
								ret=aux;
								exit_condition=1;
							}
							else if(((aux->in1==c)&&(aux->compl1==complc))&&((aux->in2==b)&&(aux->compl2==complb))&&((aux->in3==a)&&(aux->compl3==compla))) {
								ret=aux;
								exit_condition=1;
							}
							else if(((aux->in1==c)&&(aux->compl1==complc))&&((aux->in2==a)&&(aux->compl2==compla))&&((aux->in3==b)&&(aux->compl3==complb))) {
								ret=aux;
								exit_condition=1;
							}
							else if(((aux->in1==a)&&(aux->compl1!=compla))&&((aux->in2==b)&&(aux->compl2!=complb))&&((aux->in3==c)&&(aux->compl3!=complc))) {
								ret=aux;
								exit_condition=1;
								(*finalcompl)=1;
							}
							else if(((aux->in1==a)&&(aux->compl1!=compla))&&((aux->in2==c)&&(aux->compl2!=complc))&&((aux->in3==b)&&(aux->compl3!=complb))) {
								ret=aux;
								exit_condition=1;
								(*finalcompl)=1;
							}
							else if(((aux->in1==b)&&(aux->compl1!=complb))&&((aux->in2==c)&&(aux->compl2!=complc))&&((aux->in3==a)&&(aux->compl3!=compla))) {
								ret=aux;
								exit_condition=1;
								(*finalcompl)=1;
							}
							else if(((aux->in1==b)&&(aux->compl1!=complb))&&((aux->in2==a)&&(aux->compl2!=compla))&&((aux->in3==c)&&(aux->compl3!=complc))) {
								ret=aux;
								exit_condition=1;
								(*finalcompl)=1;
							}
							else if(((aux->in1==c)&&(aux->compl1!=complc))&&((aux->in2==b)&&(aux->compl2!=complb))&&((aux->in3==a)&&(aux->compl3!=compla))) {
								ret=aux;
								exit_condition=1;
								(*finalcompl)=1;
							}
							else if(((aux->in1==c)&&(aux->compl1!=complc))&&((aux->in2==a)&&(aux->compl2!=compla))&&((aux->in3==b)&&(aux->compl3!=complb))) {
								ret=aux;
								exit_condition=1;
								(*finalcompl)=1;
							}
							
						}
						k++;
					}
				}
				j++;
			}
		}
		else {
			exit_condition=1;
		}
		i++;
	}

	return ret;
}/*end of checkandret*/

/* function getcriticalpath

Synopsis [getcriticalpath]

Description [getcriticalpath]

Side effects []

*/ 
MAJ3 **getcriticalpath(MIG *mig, unsigned int *Nx) {
	MAJ3 **ret;
	MAJ3 *last, *second, *first;
	MAJ3 *current = NULL;
	unsigned int lastcompl, secondcompl, firstcompl;
	unsigned int wc, Nelem, i;
	wc=0;

	for(i=0;i<mig->Nout;i++) {
		if(mig->out[i]->level>wc){
			wc=mig->out[i]->level;
			current=mig->out[i];
		}
	}	
	Nelem=1;//at least one
	ret=(MAJ3 **)malloc(sizeof(MAJ3 *)*Nelem);
	ret[Nelem-1]=current;
	
	while(current->PI==0) {
		get_order(current,&last,&lastcompl,&second,&secondcompl,&first,&firstcompl);
		current=last;
		ret=(MAJ3 **)realloc(ret,sizeof(MAJ3 *)*Nelem);
		ret[Nelem-1]=current;
	}
	
	(*Nx)=Nelem;
	return ret;
}/*end of getcriticalpath*/

/* function addfanout

Synopsis [addfanout]

Description [addfanout]

Side effects []

*/ 
void addfanout(MAJ3* source, MAJ3 *tobeadd) {
	
//	printf("node adding %u\n",(unsigned int)tobeadd);
	
	if(source->fanout==0) {
		source->fanout++;
		source->outEdges=(MAJ3**)malloc(sizeof(MAJ3 *));
		source->outEdges[0]=tobeadd;	
	} else {
		source->fanout++;
		source->outEdges=(MAJ3**)realloc(source->outEdges,
				sizeof(MAJ3*)*(source->fanout));
		source->outEdges[source->fanout-1]=tobeadd;	
	}
	
//	for(i=0;i<source->fanout;i++) {
//		printf("node outedge %u: %u\n",i, (unsigned int)source->outEdges[i]);
//	}

}/*end of addfanout*/

/* function deletefanout

Synopsis [deletefanout]

Description [deletefanout]

Side effects []

*/ 
void deletefanout(MAJ3* source, MAJ3 *tobedel) {
	unsigned int i, found, j;
	
		found=i=0;
		while((found==0)&&(i<source->fanout)) {
			if(source->outEdges[i]==tobedel) {
				found=1;
			}
			else {
				i++;	
			}
		}
	
		//free(array[i]->outEdges);
		//free(array[i]);
	
		if(found==0) {
			printf("no found no good\n");
			printf("source level %u, fanout %u, PI %u, label %u\n",source->level, source->fanout, source->PI, source->label);
			printf("tobedel level %u, fanout %u, PI %u, label %u\n",tobedel->level, tobedel->fanout, tobedel->PI, tobedel->label);
		}
	
		for(j=i;j<source->fanout-1;j++) {
			source->outEdges[j]=source->outEdges[j+1];
		}

		source->outEdges[source->fanout-1]=NULL;
		source->outEdges=(MAJ3 **)realloc(source->outEdges,sizeof(MAJ3 *)*(source->fanout-1));
		source->fanout--;
		
		if(source->fanout==0) {
			free(source->outEdges);
			source->outEdges = NULL;
		}

}/*end of deletefanout*/

/* function recursive_polish

Synopsis [recursive_polish]

Description [recursive_polish]

Side effects []

*/ 
void recursive_polish(MIG *mig, MAJ3 *fanzero) {
	MAJ3 *in1, *in2, *in3;
	//unsigned int i;
	
	if((fanzero->PI==0)&&(fanzero->PO==0)) {
		//printf("fanzero %u -- level %u, fanout %u, PI %u, label %u\n\n", (unsigned int) fanzero, fanzero->level, fanzero->fanout, fanzero->PI, fanzero->label);
		
		in1=fanzero->in1;
		in2=fanzero->in2;
		in3=fanzero->in3;
		
		//printf("in1\n");
		deletefanout(in1, fanzero);
		//printf("in2\n");
		
		//printf("in3\n");
				
		//printf("fanout integrity: %u\n",integrity_fanout_check(mig));	
		
		if((in1->fanout==0)&&(in1->PI==0)&&(in1->PO==0)) {
			//printf("in1 recursion\n");
			recursive_polish(mig,in1);
		}
	
		deletefanout(in2, fanzero);
	
		if((in2->fanout==0)&&(in2->PI==0)&&(in2->PO==0)) {
			//printf("in2 recursion\n");
			recursive_polish(mig,in2);
		}
	
		deletefanout(in3, fanzero);
	
		if((in3->fanout==0)&&(in3->PI==0)&&(in3->PO==0)) {
			//printf("in2 recursion\n");
			recursive_polish(mig,in3);
		}
		//printf("going back\n");
		
		deletemignode(mig, fanzero);
		
	}

}/*end of recursive_polish*/

/* function deletemignode

Synopsis [deletemignode]

Description [deletemignode]

Side effects []

*/ 
void deletemignode(MIG *mig, MAJ3 *tobedel) {
	unsigned int i, found, j;
	
	if((tobedel->PI==0)&&(tobedel->PO==0)) {
		found=i=0;
		while((found==0)&&(i<mig->Nnodes)) {
			if(mig->nodes[i]==tobedel) {
				found=1;
			}
			else {
				i++;	
			}
		}
	
		if(found==0) {
			printf("no found no good\n");
		}
	
		for(j=i;j<mig->Nnodes-1;j++) {
			mig->nodes[j]=mig->nodes[j+1];
			mig->nodes[j]->label=j;
		}
		
		mig->nodes[mig->Nnodes-1]=NULL;
		mig->nodes=(MAJ3 **)realloc(mig->nodes,sizeof(MAJ3 *)*(mig->Nnodes-1));
		mig->Nnodes--;
		
		/*for(i=0;i<mig->Nnodes;i++) {
			if(mig->nodes[i]->in1==tobedel) {
				printf("danger\n");
				printf("mignodes[i] -- level %u, fanout %u, PI %u, label %u\n", mig->nodes[i]->level, mig->nodes[i]->fanout, mig->nodes[i]->PI, mig->nodes[i]->label);
				printf("tobedel -- level %u, fanout %u, PI %u, label %u\n", tobedel->level, tobedel->fanout, tobedel->PI, tobedel->label);
			}
			if(mig->nodes[i]->in2==tobedel) {
				printf("danger\n");
				printf("mignodes[i] -- level %u, fanout %u, PI %u, label %u\n", mig->nodes[i]->level, mig->nodes[i]->fanout, mig->nodes[i]->PI, mig->nodes[i]->label);
				printf("tobedel -- level %u, fanout %u, PI %u, label %u\n", tobedel->level, tobedel->fanout, tobedel->PI, tobedel->label);
			}
			if(mig->nodes[i]->in3==tobedel) {
				printf("danger\n");
				printf("mignodes[i] -- level %u, fanout %u, PI %u, label %u\n", mig->nodes[i]->level, mig->nodes[i]->fanout, mig->nodes[i]->PI, mig->nodes[i]->label);
				printf("tobedel -- level %u, fanout %u, PI %u, label %u\n", tobedel->level, tobedel->fanout, tobedel->PI, tobedel->label);
			}
		}*/
			
		free(tobedel);
	}
		
}/*end of deletemignode*/
	
/* function new_node

Synopsis [new_node]

Description [new_node]

Side effects []

*/ 

MAJ3 *new_node(MIG *net, MAJ3 *parent, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2, MAJ3 *in3, unsigned int compl3) {
	MAJ3 *ret;
	ret=NULL;	
	
	ret=(MAJ3 *)malloc(sizeof(MAJ3));
	ret->in1=in1;
	addfanout(in1, ret);
	ret->in2=in2;
	addfanout(in2, ret);
	ret->in3=in3;
	addfanout(in3, ret);
	ret->fanout=1;
	ret->outEdges=(MAJ3 **)malloc(sizeof(MAJ3 *)*ret->fanout);
	ret->outEdges[0]=parent;
	ret->value=2;
	ret->flag=0;
	ret->compl1=compl1;
	ret->compl2=compl2;
	ret->compl3=compl3;
	ret->PI=0;
	ret->PO=0;

	ret->level = (in1->level > in2->level) ? in1->level : in2->level;
	ret->level = (ret->level > in3->level) ? ret->level : in3->level;
	++ret->level;
	
	net->nodes=(MAJ3 **)realloc(net->nodes,sizeof(MAJ3 *)*(net->Nnodes+1));
	net->nodes[net->Nnodes]=ret;
	ret->label=net->Nnodes;
	net->Nnodes++;
		
	return ret;
}/*end of new_node*/

/* function new_node_nop

Synopsis [new_node no parent]

Description [new_node]

Side effects []

*/ 

MAJ3 *new_node_nop(MIG *net, MAJ3 *in1, unsigned int compl1, MAJ3 *in2, unsigned int compl2, MAJ3 *in3, unsigned int compl3) {
	MAJ3 *ret;
	ret=NULL;	
	
	ret=(MAJ3 *)malloc(sizeof(MAJ3));
	ret->in1=in1;
	addfanout(in1, ret);
	ret->in2=in2;
	addfanout(in2, ret);
	ret->in3=in3;
	addfanout(in3, ret);
	ret->fanout=0;
	//ret->outEdges=(MAJ3 **)malloc(sizeof(MAJ3 *)*ret->fanout);
	//ret->outEdges[0]=parent;
	ret->value=2;
	ret->flag=0;
	ret->compl1=compl1;
	ret->compl2=compl2;
	ret->compl3=compl3;
	ret->PI=0;
	ret->PO=0;
	
	net->nodes=(MAJ3 **)realloc(net->nodes,sizeof(MAJ3 *)*(net->Nnodes+1));
	net->nodes[net->Nnodes]=ret;
	ret->label=net->Nnodes;
	net->Nnodes++;
		
	return ret;
}/*end of new_node_nop*/

/* function 

Synopsis [integrity_fanout_check]

Description [integrity_fanout_check]

Side effects []

*/ 

unsigned int integrity_fanout_check(MIG *mig) {
	unsigned int i, found, j, totfanout;
	unsigned int ret;
	ret=1;
	totfanout=0;
	
	for(i=0;i<mig->Nnodes;i++) {
		found=0;
		totfanout+=mig->nodes[i]->fanout;
		for(j=0;j<mig->nodes[i]->in1->fanout;j++) {
			if(mig->nodes[i]->in1->outEdges[j]==mig->nodes[i]) {
				found=1;
			}
		}
		if(found==0) {
			ret=0;
			printf("nodes->in1 %u, nodes %u\n", mig->nodes[i]->in1->label, mig->nodes[i]->label);
		}
		found=0;
		for(j=0;j<mig->nodes[i]->in2->fanout;j++) {
			if(mig->nodes[i]->in2->outEdges[j]==mig->nodes[i]) {
				found=1;
			}
		}
		if(found==0) {
			ret=0;
			printf("nodes->in2 %u, nodes %u\n", mig->nodes[i]->in2->label, mig->nodes[i]->label);
			for(j=0;j<mig->nodes[i]->in2->fanout;j++) {
				printf("node outedge %u %u\n",j, mig->nodes[i]->in2->outEdges[j]->label);
			}
		}
		found=0;
		for(j=0;j<mig->nodes[i]->in3->fanout;j++) {
			if(mig->nodes[i]->in3->outEdges[j]==mig->nodes[i]) {
				found=1;
			}
		}
		if(found==0) {
			ret=0;
			printf("nodes->in3 %u, nodes %u\n", mig->nodes[i]->in3->label, mig->nodes[i]->label);
		}
	}
	
	for(i=0;i<mig->Nin;i++) {
		totfanout+=mig->in[i]->fanout;
		if(mig->in[i]->fanout==0) {
			//printf("input %u disconnected\n",i);
		}
	}
	
	totfanout+=mig->one->fanout;
	
	if(totfanout!=mig->Nnodes*3) {
		printf("Attention: total number of edges: %u, total number of nodesx 3: %u\n",totfanout,mig->Nnodes*3);
	}
	
	return ret;
}/*end of integrity_fanout_check*/

/* function setrecursiveflag

Synopsis [setrecursiveflag]

Description [setrecursiveflag]

Side effects []

*/ 
unsigned int setrecursiveflag(MAJ3 *node) {

	unsigned int ret;
	ret=0;
		
	if((node->PI==0)&&(node->flag==0)) {
		node->flag=1;
		setrecursiveflag(node->in1);
		setrecursiveflag(node->in2);
		setrecursiveflag(node->in3);
	}
	
	ret=1;
	return ret;
}/* end setrecursiveflag */	


/* function setrecursiveflag

Synopsis [setrecursiveflagfanout1]

Description [setrecursiveflagfanout1]

Side effects []

*/ 
unsigned int setrecursiveflagfanout1(MAJ3 *node) {

	unsigned int ret;
	ret=0;
		
	if((node->PI==0)&&(node->PO==0)&&(node->flag==0)&&(node->fanout==1)) {
		node->flag=1;
		setrecursiveflagfanout1(node->in1);
		setrecursiveflagfanout1(node->in2);
		setrecursiveflagfanout1(node->in3);
	}
	
	ret=1;
	return ret;
}/* end setrecursiveflagfanout1 */	

/* function recursiveduplicate

Synopsis [recursiveduplicate]

Description [recursiveduplicate]

Side effects []

*/ 
MAJ3 *recursiveduplicate(MIG *net,MAJ3 *node) {

	MAJ3 *ret, *in1, *in2, *in3;
	ret=NULL;

	if(node->PI==1) {
		ret=node;//we do not duplicate PI
	}
	else if(node->flag==0) {
		in1=recursiveduplicate(net,node->in1);
		in2=recursiveduplicate(net,node->in2);
		in3=recursiveduplicate(net,node->in3);
		ret=new_node_nop(net, in1, node->compl1, in2, node->compl2, in3, node->compl3);
		node->aux=ret;
		node->flag=1;
	}
	else {
		if(node->aux==NULL) {
			printf("problem in recursiveduplicate\n");
		}
		else {
			ret=node->aux;
		}
	}

	return ret;
}/* end recursiveduplicate */	

/* function duplicatesubtree

Synopsis [duplicatesubtree]

Description [duplicatesubtree]

Side effects []

*/ 
MAJ3 *duplicatesubtree(MIG *net, MAJ3 *node) {

	MAJ3 *ret;
	ret=NULL;
	
	resetflag(net);
	resetaux(net);
		
	ret=recursiveduplicate(net,node);
	
	resetflag(net);
	resetaux(net);
		
		
	//here to create a dummy output to support depth computation	
	net->Nout++;	
	net->out=realloc(net->out,sizeof(MAJ3 *)*net->Nout);
	net->out[net->Nout-1]=ret;
	ret->PO=1;	
		
	computedepth(net);
		
	ret->PO=0;		
	//here dummy output is deleted	
	net->out[net->Nout-1]=NULL;
	net->Nout--;	
	net->out=realloc(net->out,sizeof(MAJ3 *)*net->Nout);	
		
	return ret;
}/* end duplicatesubtree */	


/* function checkduplicateredundancy

Synopsis [checkduplicateredundancy]

Description [checkduplicateredundancy]

Side effects []

*/ 
unsigned int checkduplicateredundancy(MIG *mig) {

	unsigned int ret, i, j;
	ret=0;
	
	for(i=0;i<mig->Nnodes;i++) {
		for(j=i+1;j<mig->Nnodes;j++) {
			if(mig->nodes[i]!=mig->nodes[j]) {
				if((mig->nodes[i]->in1==mig->nodes[j]->in1)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl1)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl3)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl3)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in3)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl3)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl1)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl1)) {
					ret++;		
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in2)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl2)) {
					ret++;			
				}
			}
		}
	}

	
	return ret;
}/* end checkduplicateredundancy */	

/* function checkallduplicateredundancy

Synopsis [checkallduplicateredundancy]

Description [checkallduplicateredundancy]

Side effects []

*/ 
unsigned int checkallduplicateredundancy(MIG *mig) {

	unsigned int ret, i, j;
	ret=0;
	
	for(i=0;i<mig->Nnodes;i++) {
		for(j=i+1;j<mig->Nnodes;j++) {
			if(mig->nodes[i]!=mig->nodes[j]) {
				if((mig->nodes[i]->in1==mig->nodes[j]->in1)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl1)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl3)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl3)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in3)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl3)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl1)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl1)) {
					ret++;		
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in2)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl2)) {
					ret++;			
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in1)&&(mig->nodes[i]->compl1!=mig->nodes[j]->compl1)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2!=mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3!=mig->nodes[j]->compl3)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1!=mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2!=mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3!=mig->nodes[j]->compl3)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1!=mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in3)&&(mig->nodes[i]->compl2!=mig->nodes[j]->compl3)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3!=mig->nodes[j]->compl1)) {
					ret++;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1!=mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2!=mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3!=mig->nodes[j]->compl1)) {
					ret++;		
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1!=mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2!=mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in2)&&(mig->nodes[i]->compl3!=mig->nodes[j]->compl2)) {
					ret++;			
				}
			}
		}
	}

	
	return ret;
}/* end checkallduplicateredundancy */	


/* function eliminateduplicateredundancy

Synopsis [eliminateduplicateredundancy]

Description [eliminateduplicateredundancy]

Side effects []

*/ 
unsigned int eliminateduplicateredundancy(MIG *mig) {

	unsigned int ret, i, j, flag, fanout;
	MAJ3 *current, *node, *nodenew;
	MAJ3 **list;
	ret=0;
	
	for(i=0;i<mig->Nnodes;i++) {
		for(j=i+1;j<mig->Nnodes;j++) {
			if(mig->nodes[i]!=mig->nodes[j]) {
				flag=0;
				if((mig->nodes[i]->in1==mig->nodes[j]->in1)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl1)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl3)) {
					flag=1;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in3)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl3)) {
					flag=1;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in2)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl2)&&(mig->nodes[i]->in2==mig->nodes[j]->in3)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl3)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl1)) {
					flag=1;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in2)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl2)&&(mig->nodes[i]->in3==mig->nodes[j]->in1)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl1)) {
					flag=1;
				}
				else if((mig->nodes[i]->in1==mig->nodes[j]->in3)&&(mig->nodes[i]->compl1==mig->nodes[j]->compl3)&&(mig->nodes[i]->in2==mig->nodes[j]->in1)&&(mig->nodes[i]->compl2==mig->nodes[j]->compl1)&&(mig->nodes[i]->in3==mig->nodes[j]->in2)&&(mig->nodes[i]->compl3==mig->nodes[j]->compl2)) {
					flag=1;
				}
				if(flag==1) {
					node=mig->nodes[j];
					nodenew=mig->nodes[i];
					list=(MAJ3 **)malloc(sizeof(MAJ3 *)*node->fanout);
					fanout=node->fanout;
					for(i=0;i<fanout;i++) {
						list[i]=node->outEdges[i];
					}
					for(i=0;i<fanout;i++) {
						current=list[i];
						if(current->in1==node) {
							current->in1=nodenew;
						} 
						else if(current->in2==node) {
							current->in2=nodenew;
						}
						else {
							current->in3=nodenew;
						}
						addfanout(nodenew, current);
						deletefanout(node, current);
					}
					free(list);
					deletefanout(node->in1, node);
					deletefanout(node->in2, node);
					deletefanout(node->in3, node);
					if(node->fanout==0) {
						deletemignode(mig,node);
					}	
					else {
						printf("problem\n");
						printf("fanout %u, level %u\n",node->fanout,node->level);
					}
				}
			}
		}
	}

	
	return ret;
}/* end eliminateduplicateredundancy */	


/* function substitute_in_flag

Synopsis [substitute_in_flag]

Description [substitute_in_flag]

Side effects [apply only on irredundant MIGs]

*/ 
unsigned int substitute_in_flag(MIG *mig, MAJ3 *old, MAJ3 *nodenew, unsigned int compl_sub) {

	unsigned int ret, i, fanout;
	MAJ3 *current;
	MAJ3 **list;
	ret=0;
	
	list=(MAJ3 **)malloc(sizeof(MAJ3 *)*old->fanout);
	
	fanout=old->fanout;
	for(i=0;i<fanout;i++) {
		list[i]=old->outEdges[i];
	}	
	
	for(i=0;i<fanout;i++) {
		current=list[i];
		//printf("outside %u/%u\n",i,fanout);
		//printf("current -- level %u, fanout %u, PI %u, label %u\n", current->level, current->fanout, current->PI, current->label);
		if((current!=NULL)&&(current->flag==1)) {
			ret=1;
			//printf("inside %u\n",i);
			if(current->in1==old) {
				if(compl_sub==1) {
					if(current->compl1==1) {
						current->compl1=0;
					}
					else {
						current->compl1=1;
					}
				}
				current->in1=nodenew;
				addfanout(nodenew,current);
				deletefanout(old, current);	
			} 
			else if(current->in2==old) {
				if(compl_sub==1) {
					if(current->compl2==1) {
						current->compl2=0;
					}
					else {
						current->compl2=1;
					}
				}
				current->in2=nodenew;
				addfanout(nodenew,current);
				deletefanout(old, current);	
			}
			else {//it is in3 - no other choice
				if(compl_sub==1) {
					if(current->compl3==1) {
						current->compl3=0;
					}
					else {
						current->compl3=1;
					}
				}
				current->in3=nodenew;
				addfanout(nodenew,current);
				deletefanout(old, current);	
			}
		}
	}
	
//	printf("totally outside\n");
	
	free(list);
		
	if((old->fanout==0)&&(old->PO==0)&&(old->PI==0)) {
		recursive_polish(mig,old);
	}	
	
	if(ret==1) {
		computedepth(mig);
	}

	return ret;
}/* end substitute_in_flag */	

/* function maxdepth

Synopsis [maxdepth]

Description [return the maximum depth of the MIG]

Side effects []

*/ 
unsigned int maxdepth(MIG *mig) {
	unsigned int ret, wc, i;
	
	wc=0;
	for(i=0;i<mig->Nout;i++) {
		if(mig->out[i]->level>wc){
			wc=mig->out[i]->level;
		}
	}
	
	ret=wc;

	return ret;
}/* end maxdepth */	

/* function realmaj

Synopsis [realmaj]

Description [return the number of non-trivial MAJ nodes]

Side effects []

*/ 
unsigned int realmaj(MIG *mig) {
	unsigned int ret, wc, i;
	
	wc=0;
	for(i=0;i<mig->Nnodes;i++) {
		if((mig->nodes[i]->in1!=mig->one)&&(mig->nodes[i]->in2!=mig->one)&&(mig->nodes[i]->in3!=mig->one)) {
			wc++;
		}
	}
	
	ret=wc;

	return ret;
}/* end realmaj */	


/* function noduplicates

Synopsis [noduplicates]

Description [eliminateduplicateredundancy]

Side effects []

*/ 
unsigned int noduplicates(MIG *mig) {

	unsigned int ret, i, j, flag, fanout, allparsed, counter, compl;
	MAJ3 *current, *node, *nodenew;
	MAJ3 **list;
	MAJ3 **array;
	ret=0;
	
	resetflag(mig);
	
	for(i=0;i<mig->Nin;i++) {
		mig->in[i]->flag=1;//inputs are never duplicated
	}
	mig->one->flag=1;//constants not duplicated
	
	allparsed=0;
	while(allparsed==0) {
		counter=0;
		array=(MAJ3 **)malloc(sizeof(MAJ3 *)*counter);
		for(i=0;i<mig->Nnodes;i++) {
			if((mig->nodes[i]->in1->flag==1)&&(mig->nodes[i]->in2->flag==1)&&(mig->nodes[i]->in3->flag==1)) {
				counter++;
				array=(MAJ3 **)realloc(array,sizeof(MAJ3 *)*counter);
				array[counter-1]=mig->nodes[i];
			}
		}
		
		for(i=0;i<counter;i++) {
			for(j=i+1;j<counter;j++) {
				if((array[j]!=NULL)&&(array[i]!=NULL)) {
					flag=0;
					compl=0;
					if((array[i]->in1==array[j]->in1)&&(array[i]->compl1==array[j]->compl1)&&(array[i]->in2==array[j]->in2)&&(array[i]->compl2==array[j]->compl2)&&(array[i]->in3==array[j]->in3)&&(array[i]->compl3==array[j]->compl3)) {
						flag=1;
					}
					else if((array[i]->in1==array[j]->in2)&&(array[i]->compl1==array[j]->compl2)&&(array[i]->in2==array[j]->in1)&&(array[i]->compl2==array[j]->compl1)&&(array[i]->in3==array[j]->in3)&&(array[i]->compl3==array[j]->compl3)) {
						flag=1;
					}
					else if((array[i]->in1==array[j]->in2)&&(array[i]->compl1==array[j]->compl2)&&(array[i]->in2==array[j]->in3)&&(array[i]->compl2==array[j]->compl3)&&(array[i]->in3==array[j]->in1)&&(array[i]->compl3==array[j]->compl1)) {
						flag=1;
					}
					else if((array[i]->in1==array[j]->in3)&&(array[i]->compl1==array[j]->compl3)&&(array[i]->in2==array[j]->in2)&&(array[i]->compl2==array[j]->compl2)&&(array[i]->in3==array[j]->in1)&&(array[i]->compl3==array[j]->compl1)) {
						flag=1;
					}
					else if((array[i]->in1==array[j]->in3)&&(array[i]->compl1==array[j]->compl3)&&(array[i]->in2==array[j]->in1)&&(array[i]->compl2==array[j]->compl1)&&(array[i]->in3==array[j]->in2)&&(array[i]->compl3==array[j]->compl2)) {
						flag=1;
					}
					else if((array[i]->in1==array[j]->in1)&&(array[i]->compl1!=array[j]->compl1)&&(array[i]->in2==array[j]->in2)&&(array[i]->compl2!=array[j]->compl2)&&(array[i]->in3==array[j]->in3)&&(array[i]->compl3!=array[j]->compl3)) {
						flag=1;
						compl=1;
					}
					else if((array[i]->in1==array[j]->in2)&&(array[i]->compl1!=array[j]->compl2)&&(array[i]->in2==array[j]->in1)&&(array[i]->compl2!=array[j]->compl1)&&(array[i]->in3==array[j]->in3)&&(array[i]->compl3!=array[j]->compl3)) {
						flag=1;
						compl=1;
					}
					else if((array[i]->in1==array[j]->in2)&&(array[i]->compl1!=array[j]->compl2)&&(array[i]->in2==array[j]->in3)&&(array[i]->compl2!=array[j]->compl3)&&(array[i]->in3==array[j]->in1)&&(array[i]->compl3!=array[j]->compl1)) {
						flag=1;
						compl=1;
					}
					else if((array[i]->in1==array[j]->in3)&&(array[i]->compl1!=array[j]->compl3)&&(array[i]->in2==array[j]->in2)&&(array[i]->compl2!=array[j]->compl2)&&(array[i]->in3==array[j]->in1)&&(array[i]->compl3!=array[j]->compl1)) {
						flag=1;
						compl=1;
					}
					else if((array[i]->in1==array[j]->in3)&&(array[i]->compl1!=array[j]->compl3)&&(array[i]->in2==array[j]->in1)&&(array[i]->compl2!=array[j]->compl1)&&(array[i]->in3==array[j]->in2)&&(array[i]->compl3!=array[j]->compl2)) {
						flag=1;
						compl=1;
					}
					if(flag==1) {
						node=array[j];
						nodenew=array[i];
						list=(MAJ3 **)malloc(sizeof(MAJ3 *)*node->fanout);
						fanout=node->fanout;
						for(i=0;i<fanout;i++) {
							list[i]=node->outEdges[i];
						}
						for(i=0;i<fanout;i++) {
							current=list[i];
							if(current->in1==node) {
								current->in1=nodenew;
								if(compl==1) {
									if(current->compl1==1) {
										current->compl1=0;
									}
									else {
										current->compl1=1;
									}
								}
							} 
							else if(current->in2==node) {
								current->in2=nodenew;
								if(compl==1) {
									if(current->compl2==1) {
										current->compl2=0;
									}
									else {
										current->compl2=1;
									}
								}
							}
							else {
								current->in3=nodenew;
								if(compl==1) {
									if(current->compl3==1) {
										current->compl3=0;
									}
									else {
										current->compl3=1;
									}
								}
							}
							addfanout(nodenew, current);
							deletefanout(node, current);
						}
						free(list);
						deletefanout(node->in1, node);
						deletefanout(node->in2, node);
						deletefanout(node->in3, node);
						if(node->fanout==0) {
							deletemignode(mig,node);
						}	
						else {
							printf("problem\n");
							printf("fanout %u, level %u\n",node->fanout,node->level);
						}
						array[j]=NULL;
					}
				}	
			}
		}
		
		for(i=0;i<counter;i++) {
			if(array[i]!=NULL) {
				array[i]->flag=1;//now it is not duplicated
			}
		}
		
		free(array);
		allparsed=1;
		for(i=0;i<mig->Nnodes;i++) {
			if(mig->nodes[i]->flag==0) {
				allparsed=0;
			}
		}
	}
	
	resetflag(mig);
	
	return ret;
}/* end noduplicates */	


/* function replacenode

Synopsis [replacenode]

Description [replace one node with the other]

Side effects [unsafe to use it without enviromental guarantees on the validity of the substitution]

*/ 
unsigned int replacenode(MIG *mig, MAJ3 *node, MAJ3 *nodenew, unsigned int compl) {
	
	unsigned int ret, i, fanout;
	MAJ3 *current;
	MAJ3 **list;
	ret=0;
	
	list=(MAJ3 **)malloc(sizeof(MAJ3 *)*node->fanout);
	fanout=node->fanout;
	for(i=0;i<fanout;i++) {
		list[i]=node->outEdges[i];
	}
	for(i=0;i<fanout;i++) {
		current=list[i];
		if(current->in1==node) {
			current->in1=nodenew;
			if(compl==1) {
				if(current->compl1==1) {
					current->compl1=0;
				}
				else {
					current->compl1=1;
				}
			}
		} 
		else if(current->in2==node) {
			current->in2=nodenew;
			if(compl==1) {
				if(current->compl2==1) {
					current->compl2=0;
				}
				else {
					current->compl2=1;
				}
			}
		}
		else {
			current->in3=nodenew;
			if(compl==1) {
				if(current->compl3==1) {
					current->compl3=0;
				}
				else {
					current->compl3=1;
				}
			}
		}
		addfanout(nodenew, current);
		deletefanout(node, current);
	}
	free(list);
	deletefanout(node->in1, node);
	deletefanout(node->in2, node);
	deletefanout(node->in3, node);
	if(node->fanout==0) {
		if(node->PO==1) {
			for(i=0;i<mig->Nout;i++) {
				if(mig->out[i]==node) {
					mig->out[i]=nodenew;
					mig->out[i]->PO=1;
				}
			}	
		}
		node->PO=0;
		deletemignode(mig,node);
	}	
	else {
		printf("problem\n");
		printf("fanout %u, level %u\n",node->fanout,node->level);
	}

	return ret;
	
}/* end replacenode */	


/* function largefanoutcriticalpath

Synopsis [largefanoutcriticalpath]

Description [return the node with large fanout on the critical path]

Side effects [set/reset the MIG flags]

*/ 
MAJ3 *largefanoutcriticalpath(MIG *net, MAJ3 *node) {
	unsigned int i, max, size;
	MAJ3 *candidate1;
	MAJ3 *current, **list;
	
	MAJ3 *ret;
	ret=NULL;
	
	resetflag(net);
	
	setrecursiveflag(node);
	
	for(i=0;i<net->Nnodes;i++) {
		if(net->nodes[i]->flag==1) {
			net->nodes[i]->in1->flag++;
			net->nodes[i]->in2->flag++;
			net->nodes[i]->in3->flag++;
		}
	}
	
	current=node;
	size=0;
	list=(MAJ3 **)malloc(sizeof(MAJ3 *)*size);
	while((current->level>0)&&(current->PI==0)) {
		size++;
		list=realloc(list,sizeof(MAJ3 *)*size);
		list[size-1]=current;
		if(current->in1->level>current->in2->level) {
			if(current->in1->level>current->in3->level) {
				current=current->in1;
			}
			else {
				current=current->in3;
			}
		}
		else {
			if(current->in2->level>current->in3->level) {
				current=current->in2;
			}
			else {
				current=current->in3;
			}
		}
	}
	
	max=1;
	candidate1=node->in1;//default
	for(i=0;i<size;i++) {
		if((list[i]->flag*list[i]->level>max)&&(list[i]->level<node->level)) {
			candidate1=list[i];
			max=list[i]->flag*list[i]->level;
		}
	}
	
	resetflag(net);
	
	free(list);
	
	ret=candidate1;
	
	return ret;
	
}/* end largefanoutcriticalpath */	


/* function setcriticaldepthflag

Synopsis [setcriticaldepthflag]

Description [setcriticaldepthflag]

Side effects []

*/ 
unsigned int setcriticaldepthflag(MAJ3 *node) {

	unsigned int ret;
	ret=0;
		
	if((node->PI==0)&&(node->flag==0)) {
		node->flag=1;
		if(node->in1->level==node->level-1) {
			setcriticaldepthflag(node->in1);
		}
		if(node->in2->level==node->level-1) {
			setcriticaldepthflag(node->in2);
		}
		if(node->in3->level==node->level-1) {
			setcriticaldepthflag(node->in3);
		}
	}
	
	ret=1;
	return ret;
}/* end setcriticaldepthflag */	

/* function searchflagwithmaxfanout

Synopsis [searchflagwithmaxfanout]

Description [searchflagwithmaxfanout]

Side effects []

*/ 
MAJ3 *searchflagwithmaxfanout(MIG *net) {
	unsigned int i, max, count, j;
	
	MAJ3 *ret;
	ret=NULL;
			
			
	max=0;
	for(i=0;i<net->Nnodes;i++) {
		if(net->nodes[i]->flag==1) {
			count=0;
			for(j=0;j<net->nodes[i]->fanout;j++) {
				if(net->nodes[i]->outEdges[j]->flag==1) {
					count++;
				}
			}
			if(count>max) {
				ret=net->nodes[i];
				max=count;
			}
		}
	}
	
	printf("fanout -- level %u, fanoutXXX %u, PI %u, label %u\n", ret->level, max, ret->PI, ret->label);
	
	return ret;
}/* end searchflagwithmaxfanout */	



/* function searchhighestoutputflag

Synopsis [searchhighestoutputflag]

Description [searchhighestoutputflag]

Side effects []

*/ 
MAJ3 *searchhighestoutputflag(MIG *net, MAJ3 *node) {
	unsigned int i, max;
	
	MAJ3 *ret;
	ret=NULL;
			
			
	max=0;
	for(i=0;i<node->fanout;i++) {
		if(node->outEdges[i]->flag==1) {
			if(node->outEdges[i]->level>max) {
				ret=node->outEdges[i];
				max=node->outEdges[i]->level;
			}
			//printf("level %u\n",node->outEdges[i]->level);
		}
	}
	
	printf("highest -- level %u, fanout %u, PI %u, label %u\n", ret->level, ret->fanout, ret->PI, ret->label);
	
	return ret;
}/* end searchflagwithmaxfanout */	


/* function searchcriticalpathdominator

Synopsis [searchcriticalpathdominator]

Description [searchcriticalpathdominator]

Side effects []

*/ 
MAJ3 *searchcriticalpathdominator(MIG *net, MAJ3 *node) {
	unsigned int i, j, exit_flag, count;
	
	MAJ3 *ret;
	ret=NULL;
	
	resetflag(net);
	
	setcriticaldepthflag(node);
	
	exit_flag=0;
	j=0;
	while((exit_flag==0)&&(j<node->level)) {
		i=0;
		count=0;
		while((exit_flag==0)&&(i<net->Nnodes)) {
			if(net->nodes[i]->flag==1) {
				if(net->nodes[i]->level==j) {
					count++;
					ret=net->nodes[i];
				}
			}
			i++;	
		}
		if(count==1) {
			exit_flag=1;
		}
		else {
			ret=NULL;
		}
		j++;
	}	
	
	resetflag(net);
	
	if((exit_flag==1)&&(ret!=NULL)) {
		//printf("dominator -- level %u, fanout %u, PI %u, label %u\n", ret->level, ret->fanout, ret->PI, ret->label);
	}
	
	return ret;
}/* end searchcriticalpathdominator */	



/* function searchneighborpair

Synopsis [searchneighborpair withing flag=1 nodes]

Description [searchneighborpair withing flag=1 nodes]

Side effects [quadratic complexity]

*/ 
unsigned int searchneighborpair(MIG *mig, MAJ3 **pair1, MAJ3 **pair2) {
	unsigned int i, j, k, m, max, count;
	unsigned int ret;
	ret=0;
	MAJ3 *nodei, *nodej;
	
	max=0;
	for(i=0;i<mig->Nin;i++) {
		mig->in[i]->flag=1;//set also the inputs
	}
	for(i=0;i<(mig->Nnodes+mig->Nin);i++) {
		for(j=i+1;j<(mig->Nnodes+mig->Nin);j++) {
			if(i<mig->Nin) {//inputs
				nodei=mig->in[i];
			}
			else {//other nodes
				nodei=mig->nodes[i-mig->Nin];
			}
			if(j<mig->Nin) {//inputs
				nodej=mig->in[j];
			}
			else {//other nodes
				nodej=mig->nodes[j-mig->Nin];
			}
			if((nodei->flag==1)&&(nodej->flag==1)) {
				count=0;
				for(k=0;k<nodei->fanout;k++) {
					for(m=0;m<nodej->fanout;m++) {
						if((nodei->outEdges[k]->flag==1)&&(nodej->outEdges[m]->flag==1)) {
							if(nodei->outEdges[k]==nodej->outEdges[m]) {
								count++;
							}	
						}
					}
				}
				if(count>max) {
					(*pair1)=nodei;
					(*pair2)=nodej;
					max=count;
				}
			}
		}
	}
	
	ret=max;
	
	return ret;
}/* end searchneighborpair */	

/* function searchneighborpair_wone

Synopsis [searchneighborpair_wone - including one node withing flag=1 nodes]

Description [searchneighborpair_wone - including one node withing flag=1 nodes]

Side effects [quadratic complexity]

*/ 
unsigned int searchneighborpair_wone(MIG *mig, MAJ3 **pair1, MAJ3 **pair2) {
	unsigned int i, j, k, m, max, count;
	unsigned int ret;
	ret=0;
	MAJ3 *nodei, *nodej;
	
	max=0;
	for(i=0;i<mig->Nin;i++) {
		mig->in[i]->flag=1;//set also the inputs
	}
	for(i=0;i<(mig->Nnodes+mig->Nin+1);i++) {
		for(j=i+1;j<(mig->Nnodes+mig->Nin);j++) {
			if(i==mig->Nnodes+mig->Nin) {//one
				nodei=mig->one;
			}
			else if(i<mig->Nin) {//inputs
				nodei=mig->in[i];
			}
			else {//other nodes
				nodei=mig->nodes[i-mig->Nin];
			}
			if(j<mig->Nin) {//inputs
				nodej=mig->in[j];
			}
			else {//other nodes
				nodej=mig->nodes[j-mig->Nin];
			}
			if((nodei->flag==1)&&(nodej->flag==1)) {
				count=0;
				for(k=0;k<nodei->fanout;k++) {
					for(m=0;m<nodej->fanout;m++) {
						if((nodei->outEdges[k]->flag==1)&&(nodej->outEdges[m]->flag==1)) {
							if(nodei->outEdges[k]==nodej->outEdges[m]) {
								count++;
							}	
						}
					}
				}
				if(count>max) {
					(*pair1)=nodei;
					(*pair2)=nodej;
					max=count;
				}
			}
		}
	}
	
	ret=max;
	
	return ret;
}/* end searchneighborpair_wone */	

/* function get_order

Synopsis [get_order]

Description [get_order]

Side effects []

*/ 
void get_order(MAJ3 *node, MAJ3 **last, unsigned int *compllast, MAJ3 **second, unsigned int *complsecond, MAJ3 **first, unsigned int *complfirst) {

	(*last)=(*second)=(*first)=NULL;
	(*compllast)=(*complsecond)=(*complfirst)=0;

	if(node->level>0) {
		if(node->in1->level>node->in2->level) {
			if(node->in3->level>node->in1->level) {
				(*last)=node->in3;
				(*compllast)=node->compl3;
				(*second)=node->in1;
				(*complsecond)=node->compl1;
				(*first)=node->in2;
				(*complfirst)=node->compl2;
			}
			else {
				(*last)=node->in1;
				(*compllast)=node->compl1;
				if(node->in3->level>node->in2->level) {
					(*second)=node->in3;
					(*complsecond)=node->compl3;
					(*first)=node->in2;
					(*complfirst)=node->compl2;
				}
				else {
					(*first)=node->in3;
					(*complfirst)=node->compl3;
					(*second)=node->in2;
					(*complsecond)=node->compl2;
				}
			}
		}
		else {
			if(node->in3->level>node->in2->level) {
				(*last)=node->in3;
				(*compllast)=node->compl3;
				(*second)=node->in2;
				(*complsecond)=node->compl2;
				(*first)=node->in1;
				(*complfirst)=node->compl1;
			}
			else {
				(*last)=node->in2;
				(*compllast)=node->compl2;
				if(node->in3->level>node->in1->level) {
					(*second)=node->in3;
					(*complsecond)=node->compl3;
					(*first)=node->in1;
					(*complfirst)=node->compl1;
				}
				else {
					(*first)=node->in3;
					(*complfirst)=node->compl3;
					(*second)=node->in1;
					(*complsecond)=node->compl1;
				}
			}
		}	
	}
}/*end of get_order*/


/* function searchrelevancepair

Synopsis [searchrelevancepair withing flag=1 nodes]

Description [searchrelevancepair withing flag=1 nodes]

Side effects [quadratic complexity]

*/ 
unsigned int searchrelevancepair(MIG *mig, MAJ3 **pair1, MAJ3 **pair2) {
	unsigned int i, k, max, count, size;
	unsigned int ret;
	ret=0;
	MAJ3 *nodei, **list, *current;
	MAJ3 *last, *second, *first;
	unsigned int lastcompl, secondcompl, firstcompl;
	
	(*pair2)=(*pair1)=mig->in[0];
	
	max=0;
	for(i=0;i<mig->Nnodes;i++) {
		nodei=mig->nodes[i];
		count=0;
		for(k=0;k<nodei->fanout;k++) {
			if(nodei->outEdges[k]->flag==1) {
					count++;	
			}
		}
		if(count>max) {
			(*pair1)=nodei;
			max=count;
		}
	}
		
	size=0;
	list=(MAJ3 **)malloc(sizeof(MAJ3 *)*size);
	for(k=0;k<(*pair1)->fanout;k++) {
		if(((*pair1)->outEdges[k]->flag==1)&&((*pair1)->outEdges[k]->PI==0)) {
			size++;
			list=realloc(list,sizeof(MAJ3 *)*size);
			list[size-1]=(*pair1)->outEdges[k];
		}
	}	
	
	max=0;
	for(i=0;i<size;i++) {
		current=list[i];
		get_order(current,&last,&lastcompl,&second,&secondcompl,&first,&firstcompl);
		if((second==(*pair1))&&(current->level>max)) {
			max=current->level;
			(*pair2)=current;
		}
	}
	
	ret=max;
	
	free(list);
	
	return ret;
}/* end searchrelevancepair */	


/* function twooverthree

Synopsis [twooverthree check if two nodes have 2 over three equal nodes (with equal compl edges)]

Description [twooverthree check if two nodes have 2 over three equal nodes (with equal compl edges) -- return 1 if ok, 0 otherwise]

Side effects [-]

*/ 
unsigned int twooverthree(MIG *mig, MAJ3 *node1, MAJ3 *node2, MAJ3 **common1, unsigned int *compl1, MAJ3 **common2, unsigned int *compl2, MAJ3 **uncommon1, unsigned int *complun1, MAJ3 **uncommon2, unsigned int *complun2) {
	unsigned int ret;
	ret=0;
	
	(*common1)=(*common2)=(*uncommon1)=(*uncommon2)=NULL;
	(*compl1)=(*compl2)=(*complun1)=(*complun2)=0;
	
	if((node1->PI==1)||(node2->PI==1)) {
		return 0;
	}
	
	if((node1->in1==node2->in1)&&(node1->compl1==node2->compl1)) {
		if((node1->in2==node2->in2)&&(node1->compl2==node2->compl2)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in2;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in2==node2->in3)&&(node1->compl2==node2->compl3)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in2;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in3==node2->in2)&&(node1->compl3==node2->compl2)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in3;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in3==node2->in3)&&(node1->compl3==node2->compl3)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in3;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
	}
	else if((node1->in2==node2->in2)&&(node1->compl2==node2->compl2)) {
		if((node1->in1==node2->in1)&&(node1->compl1==node2->compl1)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in1;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in3==node2->in3)&&(node1->compl3==node2->compl3)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in3;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in3==node2->in1)&&(node1->compl3==node2->compl1)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in3;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in1==node2->in3)&&(node1->compl1==node2->compl3)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in1;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
	}
	else if((node1->in3==node2->in3)&&(node1->compl3==node2->compl3)) {
		if((node1->in1==node2->in1)&&(node1->compl1==node2->compl1)) {
			ret=1;
			(*common1)=node1->in3;
			(*common2)=node1->in1;
			(*compl1)=node1->compl3;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in2==node2->in2)&&(node1->compl2==node2->compl2)) {
			ret=1;
			(*common1)=node1->in3;
			(*common2)=node1->in2;
			(*compl1)=node1->compl3;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in1==node2->in2)&&(node1->compl1==node2->compl2)) {
			ret=1;
			(*common1)=node1->in3;
			(*common2)=node1->in1;
			(*compl1)=node1->compl3;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in2==node2->in1)&&(node1->compl1==node2->compl1)) {
			ret=1;
			(*common1)=node1->in3;
			(*common2)=node1->in2;
			(*compl1)=node1->compl3;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
	}
	else if((node1->in1==node2->in2)&&(node1->compl1==node2->compl2)) {
		if((node1->in2==node2->in1)&&(node1->compl2==node2->compl1)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in2;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in2==node2->in3)&&(node1->compl2==node2->compl3)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in2;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in3==node2->in3)&&(node1->compl3==node2->compl3)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in3;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in3==node2->in1)&&(node1->compl3==node2->compl1)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in3;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
	}
	else if((node1->in1==node2->in3)&&(node1->compl1==node2->compl3)) {
		if((node1->in2==node2->in2)&&(node1->compl2==node2->compl2)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in2;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in2==node2->in1)&&(node1->compl2==node2->compl1)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in2;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl2;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in3==node2->in1)&&(node1->compl3==node2->compl1)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in3;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in3==node2->in2)&&(node1->compl3==node2->compl2)) {
			ret=1;
			(*common1)=node1->in1;
			(*common2)=node1->in3;
			(*compl1)=node1->compl1;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in2;
			(*complun1)=node1->compl2;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
	}
	else if((node1->in2==node2->in1)&&(node1->compl2==node2->compl1)) {
		if((node1->in1==node2->in2)&&(node1->compl1==node2->compl2)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in1;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in1==node2->in3)&&(node1->compl1==node2->compl3)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in1;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in3==node2->in2)&&(node1->compl3==node2->compl2)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in3;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in3;
			(*complun2)=node2->compl3;
		}
		else if((node1->in3==node2->in3)&&(node1->compl3==node2->compl3)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in3;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
	}
	else if((node1->in2==node2->in3)&&(node1->compl2==node2->compl3)) {
		if((node1->in1==node2->in1)&&(node1->compl1==node2->compl1)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in1;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in1==node2->in2)&&(node1->compl1==node2->compl2)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in1;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl1;
			(*uncommon1)=node1->in3;
			(*complun1)=node1->compl3;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
		else if((node1->in3==node2->in1)&&(node1->compl3==node2->compl1)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in3;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in2;
			(*complun2)=node2->compl2;
		}
		else if((node1->in3==node2->in2)&&(node1->compl3==node2->compl2)) {
			ret=1;
			(*common1)=node1->in2;
			(*common2)=node1->in3;
			(*compl1)=node1->compl2;
			(*compl2)=node1->compl3;
			(*uncommon1)=node1->in1;
			(*complun1)=node1->compl1;
			(*uncommon2)=node2->in1;
			(*complun2)=node2->compl1;
		}
	}
	
	
	return ret;
}/* end searchrelevancepair */	

/* function shufflemig

Synopsis [shufflemig]

Description [shufflemig elements]

Side effects [-]

*/ 
void shufflemig(MIG *mig) {
	MAJ3 *t;
	unsigned int i, j;
		
    for(i=0;i<mig->Nnodes-1;i++) {
      j = i + rand() / (RAND_MAX / (mig->Nnodes - i) + 1);
      t = mig->nodes[j];
      mig->nodes[j] = mig->nodes[i];
      mig->nodes[i] = t;
    }
	
    for(i=0;i<mig->Nnodes;i++) {
		mig->nodes[i]->label=i;
    }
	
	
}/* end shufflemig */	




