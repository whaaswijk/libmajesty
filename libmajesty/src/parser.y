%{
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cassert>
#include "mlp.h"

using namespace std;
extern "C" int yylex();
extern "C" int yyparse();
extern "C" FILE *yyin;
extern int yylineno;
extern "C" void parse_verilog(FILE *file, MIG **m);
extern "C" unsigned int char_to_int(char *cstr);

typedef unsigned int uint;

MIG *mig = NULL;
MIG **omig = NULL;
static MAJ3* one = NULL;
static map<string,MAJ3*> inmap;
static vector<pair<string,MAJ3*>> innodes;
static vector<string> outnodes;
static map<string,MAJ3*> nodemap;
static vector<MAJ3*> nodes;
static vector<MAJ3*> out;
static vector<bool> outcompl;

void handle_input(string*);
void handle_output(string*);
void handle_wire(string*);
pair<bool, MAJ3*> *handle_node(bool, string*);
void lit_assignment(string *s, int lit); 
void unary_assignment(string *s, pair<bool,MAJ3*>*);
void binary_assignment(string*, pair<bool,MAJ3*>*, char op, pair<bool,MAJ3*>*);
void maj_assignment(string*, pair<bool,MAJ3*>*, pair<bool,MAJ3*>*, pair<bool, MAJ3*>*);
void parse_wrap(void);
void init_mig(void);
void set_outcompl(int, bool, MAJ3 *);

void yyerror(const char *s);
%}

%union {
unsigned int ival;
std::string *sval;
std::pair<bool,struct MAJ3*> *pval;
}


%token <ival> MODULE
%token <ival> INPUT
%token <ival> OUTPUT
%token <ival> WIRE
%token <ival> ASSIGN
%token <ival> ENDMODULE
%token <ival> NUM
%token <ival> LIT
%token <sval> NAME
%defines
%start circuit

%type <pval> node

%%

circuit: MODULE NAME '(' tnames ')' ';' inputs outputs wires assignments mig_init ENDMODULE { parse_wrap(); }
;

tnames: /* empty */
| tnames NAME { delete $2; }
;


inputs: INPUT ilist ';'
;

ilist: /* empty */
| ilist NAME { handle_input($2); }
;

outputs: OUTPUT olist ';'
;

olist: /* empty */
| olist NAME { handle_output($2); }
;

wires: WIRE wnames ';'

wnames: /* empty */
| wnames NAME { handle_wire($2); }
;

mig_init: /* empty */ { init_mig(); }
;

assignments: /* empty */ 
| assignments assignment
;

assignment: ASSIGN NAME '=' NUM ';' /* one = 1; */ { delete $2; }
| ASSIGN NAME '=' LIT ';' { lit_assignment($2, $4); }
| ASSIGN NAME '=' node ';' { unary_assignment($2, $4); }
| ASSIGN NAME '=' node '&' node ';' { binary_assignment($2, $4, '&', $6); }
| ASSIGN NAME '=' node '|' node ';' { binary_assignment($2, $4, '|', $6); }
| ASSIGN NAME '=' '(' node '&' node ')' '|' '(' node '&' node ')' '|' '(' node '&' node ')' ';' { maj_assignment($2, $5, $7, $19); delete $11; delete $13; delete $17; }
;

node: NAME { $$ = handle_node(false, $1); }
| '~' NAME { $$ = handle_node(true, $2); }
;


%%

void yyerror(char const *e) {
cout << "Error around line " << yylineno << ": " << e << endl;
}

int output_idx(string name) {
	auto res = -1;
	for (auto i = 0u; i < outnodes.size(); i++) {
		if (outnodes[i] == name) {
			res = i;
			break;
		}
	}	
	return res;
}

void handle_input(string* s) {
	MAJ3 *n = (MAJ3*)calloc(1, sizeof(MAJ3));
	n->PI = 1;
	n->in1 = n->in2 = n->in3 = n->aux = NULL;
	n->compl1 = n->compl2 = n->compl3 = 0;
	n->outEdges = NULL;
	n->fanout = n->flag = n->level = 0;
	n->label = innodes.size();
	n->value = 2;

	auto p = pair<string,MAJ3*>(*s,n);
	innodes.push_back(p);
	inmap[*s] = n;
	delete s;
}

void handle_output(string* s) {
	auto it = inmap.find(*s);
	if (it != inmap.end()) {
		inmap[*s]->PO = 1;
	}
	outnodes.push_back(*s);
	delete s;
}

void handle_wire(string *name) {
	if (*name == "one") {
		delete name;
		return;
	}

	int idx = nodes.size();

	MAJ3 *n = (MAJ3*)calloc(1, sizeof(MAJ3));
	n->in1 = n->in2 = n->in3 = n->aux = NULL;
	n->compl1 = n->compl2 = n->compl3 = 0;
	n->outEdges = NULL;
	n->fanout = n->flag = n->PO = n->PI = n->level = 0;
	n->label = idx;
	n->value = 2;

	nodes.push_back(n);
	nodemap[*name] = n;
	delete name;
}

pair<bool, MAJ3*> *handle_node(bool cmpl, string* s) {
	if(*s == "one") {
		return new pair<bool,MAJ3*>(cmpl, one);
	}
	string name = *s;
	delete s;


	auto it = inmap.find(name);
	if (it != inmap.end()) {
		return new pair<bool,MAJ3*>(cmpl, it->second);
	}

	return new pair<bool,MAJ3*>(cmpl, nodemap[name]);
}

void lit_assignment(string *s, int lit) {
	assert(lit == 0 || lit == 1);
	int idx = -1;
	for (auto i = 0u; i < outnodes.size(); i++) {
		if (outnodes[i] == *s) {
			idx = i;
			break;
		}
	}
	delete s;
	if (idx == -1) {
		yyerror("internal node cannot be unary assignment");
		exit(1);
	}
	one->PO = 1;
	set_outcompl(idx, lit == 0 ? true : false, one);
}

void unary_assignment(string *s, pair<bool,MAJ3*> *p) {
	int idx = -1;
	for (auto i = 0u; i < outnodes.size(); i++) {
		if (outnodes[i] == *s) {
			idx = i;
			break;
		}
	}
	delete s;
	if (idx == -1) {
		yyerror("internal node cannot be unary assignment");
		delete p;
		exit(1);
	}

	p->second->PO = 1;
	set_outcompl(idx, p->first, p->second);
	delete p;
}

void binary_assignment(string *s, pair<bool,MAJ3*> *n1, char op, pair<bool,MAJ3*> *n2) {
	string name = *s;
	delete s;

	MAJ3* n = nodemap[name];
	if (n == nullptr) {
		auto node_idx = nodes.size();
		auto out_idx = output_idx(name);
		assert(out_idx != -1);
		n = (MAJ3*)calloc(1, sizeof(MAJ3));
		n->in1 = n->in2 = n->in3 = n->aux = NULL;
		n->compl1 = n->compl2 = n->compl3 = 0;
		n->outEdges = NULL;
		n->fanout = n->flag = n->PI = n->level = 0;
		n->label = node_idx;
		n->value = 2;
		n->PO = 1;
		nodes.push_back(n);
		nodemap[name] = n;
		set_outcompl(out_idx, false, n);
	}

	addfanout(n1->second, n);
	addfanout(n2->second, n);
	n->in1 = n1->second;
	n->compl1 = (unsigned int)n1->first;
	n->in2 = n2->second;
	n->compl2 = (unsigned int)n2->first;
	addfanout(one, n);
	if (op == '&') {
		n->in3 = one;
		n->compl3 = 1;
	} else {
		n->in3 = one;
		n->compl3 = 0;
	}

	delete n1;
	delete n2;
}

void maj_assignment(string *s, pair<bool,MAJ3*> *n1, pair<bool,MAJ3*> *n2, pair<bool, MAJ3*> *n3) {
	string name = *s;
	delete s;

	MAJ3 *n = nodemap[name];
	addfanout(n1->second, n);
	addfanout(n2->second, n);
	addfanout(n3->second, n);
	n->in1 = n1->second;
	n->compl1 = (unsigned int)n1->first;
	n->in2 = n2->second;
	n->compl2 = (unsigned int)n2->first;
	n->in3 = n3->second;
	n->compl3 = (unsigned int)n3->first;

	delete n1;
	delete n2;
	delete n3;
}

void init_mig(void) {
	mig = (MIG*)calloc(1, sizeof(MIG));
	mig->Nin = innodes.size();
	mig->Nout = outnodes.size();
	mig->Nnodes = nodes.size();
	mig->innames = (char**)malloc(sizeof(char*)*mig->Nin);
	mig->outnames = (char**)malloc(sizeof(char*)*mig->Nout);
	mig->outcompl = (uint*)malloc(sizeof(uint)*mig->Nout);

	mig->in = (MAJ3**)malloc(sizeof(MAJ3*)*mig->Nin);
	uint i = 0;
	for (auto it = innodes.begin(); it != innodes.end(); ++it) {
		mig->in[i] = it->second;
		mig->innames[i] = (char*)malloc(sizeof(char)*(it->first.length()+1));
		strncpy(mig->innames[i], it->first.c_str(), (it->first.length()+1));
		++i;
	}

	i = 0;
	mig->out = (MAJ3**)malloc(sizeof(MAJ3*)*mig->Nout);
	for (auto it = outnodes.begin(); it != outnodes.end(); ++it) {
		mig->outnames[i] = (char*)malloc(sizeof(char)*(it->length()+1));
		strncpy(mig->outnames[i], it->c_str(), (it->length()+1));
		mig->out[i] = out[i];
		mig->outcompl[i] = outcompl[i];
		++i;
	}

	mig->nodes = (MAJ3**)malloc(sizeof(MAJ3*)*mig->Nnodes);
	i = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it) {
		mig->nodes[i++] = *it;
	}
	mig->one = one;
}

void set_outcompl(int idx, bool c, MAJ3 *n) {
	if (out.size() == 0) {
		out.reserve(outnodes.size());
		outcompl.reserve(outnodes.size());
	}
	out[idx] = n;
	outcompl[idx] = c ? true : false;
}

void parse_wrap(void) {
	// Clear parser information that has been collected.
	inmap.clear();
	vector<pair<string,MAJ3*>>().swap(innodes);
	vector<string>().swap(outnodes);
	vector<MAJ3*>().swap(nodes);
	nodemap.clear();
	out.clear();
	outcompl.clear();
	one = NULL;

	// Set the output MIG to the MIG that has been parsed from the input file.
	*omig = mig;
}

void parse_verilog(FILE *file, MIG **m) {
	assert(one == NULL);
	one = (MAJ3*)malloc(sizeof(MAJ3));
	one->in1=one->in2=one->in3=NULL;
	one->outEdges=NULL;//to be updated
	one->value=1;
	one->fanout=one->label=0;
	one->aux=NULL;
	one->flag=0;
	one->compl1=one->compl2=one->compl3=0;
	one->PI=1;
	one->PO = false;

	omig = m;
	yyin = file;
	yyparse();
}
