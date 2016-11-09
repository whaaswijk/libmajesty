%{
#include <iostream>
#include "parser.hpp"

extern "C" int yylex(void);
extern "C" void scan_string(const char*); 
extern "C" void delete_buffer(void);
extern "C" void init_parser(void);
extern "C" void destroy_parser(void);
extern void yy_flush_buffer(YY_BUFFER_STATE buffer);

using namespace std;

#define handle_literal() { return *yytext; } 

int handle_num(void);
int handle_lit(void);
int handle_name(void);

%}

%option yylineno

whitespace [ \t\n]
operator [&|()~;=\[\]]
comment \/\/.*\n
num [0-9]+
name (\\)?[a-zA-Z][a-zA-Z0-9]*(\[[0-9]+\])?
lit 1'b[0-1]
comma ,

%%

module { return MODULE; }
input { return INPUT; }
output { return OUTPUT; }
wire { return WIRE; }
assign { return ASSIGN; }
endmodule { return ENDMODULE; }
{num} { return handle_num(); }
{lit} { return handle_lit(); }
{name} { return handle_name(); }
{whitespace} ;
{comment} ;
{operator} { return *yytext; }
{comma} ;

%%

unsigned int char_to_int(char *cstr) {
	auto str = string(cstr);
	return stoi(str);
}

int handle_num(void) {
	yylval.ival = char_to_int(yytext);
	return NUM;
}

int handle_lit(void) {
	yylval.ival = char_to_int(yytext+3);
	return LIT;
}

int handle_name(void) {
	// Strip out escape characters in Verilog varnames
	if (yytext[0] == '\\') {
		yylval.sval = new string(yytext+1);
	} else {
		yylval.sval = new string(yytext);
	}
	return NAME;
}

int yywrap(void) {
  return 1;
}

static YY_BUFFER_STATE buf = NULL;
void scan_string(const char* str) {
    buf = yy_scan_string(str);
}

void delete_buffer() {
    yy_delete_buffer(buf);
}

void init_parser() {
}

void destroy_parser() {
    yy_flush_buffer(YY_CURRENT_BUFFER);
    //yy_delete_buffer(YY_CURRENT_BUFFER);
}
