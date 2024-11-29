%option noyywrap
%option nounput

%{
	#include <iostream>
	#include <string>

	#include "../src/geom3d_object.h"
	#include "../src/point.h"
	#include "../src/line.h"
	#include "geom3d_grammar.tab.hpp"
%}

%option yylineno

%%
"make_point"				{ return MAKE_POINT; }
"make_point_in_plane"		{ return MAKE_POINT_IN_PLANE; }
"make_midpoint"				{ return MAKE_MIDPOINT; }
"make_point_ratio"			{ return MAKE_POINT_RATIO; }
"make_point_on_line"		{ return MAKE_POINT_ON_LINE; }
"translate"					{ return TRANSLATE; }
"make_point_projection"		{ return MAKE_POINT_PROJECTION; }
"make_line_intersection"	{ return MAKE_LINE_INTERSECTION; }
"make_line_plane_intersection" { return MAKE_LINE_PLANE_INTERSECTION; }
"make_foot_on_plane"		{ return MAKE_FOOT_ON_PLANE; }
"make_line"					{ return MAKE_LINE; }
"make_line_orthogonal_on_plane" { return MAKE_LINE_ORTHOGONAL_ON_PLANE; }
"make_cube"					{ return MAKE_CUBE; }
"make_tetrahedron"			{ return MAKE_TETRAHEDRON; }
"make_regular_tetrahedron"	{ return MAKE_REGULAR_TETRAHEDRON; }
"make_regular_hexagon"		{ return MAKE_REGULAR_HEXAGON; }
"make_unit_sphere"			{ return MAKE_UNIT_SPHERE; }
"make_sphere"				{ return MAKE_SPHERE; }
"make_pyramid"				{ return MAKE_PYRAMID; }
"distance"					{ return DISTANCE; }
"square_distance"			{ return SQUARE_DISTANCE; }
"make_number"				{ return MAKE_NUMBER; }
"make_plane"				{ return MAKE_PLANE; } 
"make_plane_orthogonal_on_plane" { return MAKE_PLANE_ORTHOGONAL_ON_PLANE; }
"make_plane_orthogonal_on_line"  { return MAKE_PLANE_ORTHOGONAL_ON_LINE; }
"make_tangent_plane"			 { return MAKE_TANGENT_PLANE; }
"make_parallelogram"		{ return MAKE_PARALLELOGRAM; }
"prove"                     { return PROVE; }
"point_on_line"             { return POINT_ON_LINE; }
"equal_points"				{ return EQUAL_POINTS; }
"midpoint"					{ return MIDPOINT; }
"congruent"					{ return CONGRUENT; }
"segments_in_ratio"			{ return SEGMENTS_IN_RATIO; }
"point_segment_ratio"		{ return POINT_SEGMENT_RATIO; }
"line_intersection"			{ return LINE_INTERSECTION; }
"orthogonal_lines"			{ return ORTHOGONAL_LINES; }
"parallel_lines"			{ return PARALLEL_LINES; }
"not_skew"					{ return NOT_SKEW; }
"point_in_plane"			{ return POINT_IN_PLANE; }
"point_on_sphere"			{ return POINT_ON_SPHERE; }
"parallel_planes"			{ return PARALLEL_PLANES; }
"orthogonal_planes"			{ return ORTHOGONAL_PLANES; }
"parallel_line_plane"		{ return PARALLEL_LINE_PLANE; }
"orthogonal_line_plane"		{ return ORTHOGONAL_LINE_PLANE; }
"line_in_plane"				{ return LINE_IN_PLANE; }
"line_plane_intersection"	{ return LINE_PLANE_INTERSECTION; }
"equal_angles"				{ return EQUAL_ANGLES; }
"angle"						{ return ANGLE;}
"equal_numbers"				{ return EQUAL_NUMBERS; }
"collinear"					{ return COLLINEAR; }
[+-]?[0-9]+(\.[0-9]*)?		{ yylval.number = atoi(yytext); return number_token; }
([a-zA-Z_][a-zA-Z0-9]*)*	{ yylval.name = new std::string(yytext); return name_token; }
"="							{ return *yytext; }
"{"							{ return *yytext; }
"}"							{ return *yytext; }
"*"							{ return *yytext; }
"/"							{ return *yytext; }
"+"							{ return *yytext; }
\n							{ return *yytext; }
[ \t]                       { }
"%".*                       { /* This is comment. DO NOTHING */ }
.							{ std::cerr << "Error, unknown ' " << *yytext << "'"<< std::endl;
                              exit(EXIT_FAILURE); 
                            }



%%

/*
int main(int argc, char **argv)
{
	if (argc > 1) {
		yyin = fopen(argv[1], "r");
		if (yyin == NULL) {
			perror(argv[1]);
			return 1;
		}
	}

	yylex();
	

	return 0;
}
*/
