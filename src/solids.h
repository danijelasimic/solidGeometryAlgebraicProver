#ifndef SOLIDS_H
#define SOLIDS_H

#include <vector>

#include "line.h"
#include "point.h"
#include "polinom.h"

void make_regular_tetrahedron(Point& C, Point& D, std::vector<Polinom>& construction_polinoms);
void make_tetrahedron(Point& C, Point& D, std::vector<Polinom>& construction_polinoms, 
                        std::vector<std::shared_ptr<Variable>>& variables);

/* Helper functions for making solids base. */
void make_regular_triangle(Point& C, std::vector<Polinom>& construction_polinoms);
void make_regular_hexagon(Point& A3, Point& A4, Point& A5, Point& A6, std::vector<Polinom>& construction_polinoms);
void make_regular_octagon(Point& A2, Point& A3, Point& A4, Point& A5, Point& A6, Point& A7,
                  std::vector<Polinom>& construction_polinoms);
                  
/* Make top of the pyramid. */                  
void make_top_vertex(Point& Top, int number_base, std::vector<Polinom>& construction_polinoms);

void make_parallelogram(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& construction_polinoms);

#endif //SOLIDS_H