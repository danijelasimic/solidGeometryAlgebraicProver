#ifndef STATEMENT_H
#define STATEMENT_H

#include <vector>

#include "line.h"
#include "point.h"
#include "polinom.h"
#include "plane.h"
#include "number.h"
#include "sphere.h"

void point_on_line(Point& A, Line& l, std::vector<Polinom>& statement_polinoms);
void point_on_line(Point& A, Point& P, Point& Q, std::vector<Polinom>& statement_polinoms);

void equal_points(Point& A, Point& P, std::vector<Polinom>& statement_polinoms);
void midpoint(Point& M, Point& A, Point& B, std::vector<Polinom>& statement_polinoms);
void point_segment_ratio(Point& M, Point& A, Point& B, int m, int n, std::vector<Polinom>& statement_polinoms);
void point_segment_ratio(Point& M, Point& A, Point& B, Number& m, Number& n, std::vector<Polinom>& polinoms);

void make_foot_on_plane(Point& D, Point& A, Plane& pi, std::vector<Polinom>& polinoms);
void make_point_projection(Point& D, Point& A, Plane& pi, std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& polinoms);

void congruent(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& statement_polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms);
void segments_in_ratio(Point& A, Point& B, Point& C, Point& D, int m, int n, std::vector<Polinom>& statement_polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms);

void line_intersection(Point& A, Line& l1, Line& l2, std::vector<Polinom>& polinoms);
void line_intersection(Point& A, Point& P, Point& Q, Point& R, Point& S, std::vector<Polinom>& polinoms);

void orthogonal_lines(Line& l1, Line& l2, std::vector<Polinom>& polinoms);
void orthogonal_lines(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& polinoms);

void not_skew(Line& l1, Line& l2, std::vector<Polinom>& polinoms);
void not_skew(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& polinoms);

void parallel_lines(Line& l1, Line& l2, std::vector<Polinom>& polinoms);
void parallel_lines(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& polinoms);

void make_plane_trough_points(Plane& pi, Point& A, Point& B, Point& C, std::vector<Polinom>& construction_polinoms);
void make_vector_of_plane_trough_points(Plane& pi, Point& A, Point& B, Point& C, std::vector<Polinom>& construction_polinoms);

void point_in_plane(Point& A, Plane& pi, std::vector<Polinom>& polinoms);
void point_in_plane(Point& A, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms);

void parallel_planes(Plane& pi, Plane& teta, std::vector<Polinom>& polinoms);
void parallel_planes(Point& A, Point& B, Point& C, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms);
void parallel_planes(Plane& pi, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms);

void orthogonal_planes(Plane& pi, Plane& teta, std::vector<Polinom>& polinoms);
void orthogonal_planes(Point& A, Point& B, Point& C, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms);
void orthogonal_planes(Plane& pi, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms);


void parallel_line_plane(Line& l, Plane& pi, std::vector<Polinom>& polinoms);
void parallel_line_plane(Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms);
void parallel_line_plane(Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms,
							std::vector<std::shared_ptr<Variable>>& variables,
							std::vector<Polinom>& construction_polinoms);
void parallel_line_plane(Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms,
							std::vector<std::shared_ptr<Variable>>& variables,
							std::vector<Polinom>& construction_polinoms);

void orthogonal_line_plane(Line& l, Plane& pi, std::vector<Polinom>& polinoms);
void orthogonal_line_plane(Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms);
void orthogonal_line_plane(Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms);
void orthogonal_line_plane(Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms);

void line_in_plane(Line& l, Plane& pi, std::vector<Polinom>& polinoms);
void line_in_plane(Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms);
void line_in_plane(Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms);
void line_in_plane(Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms);

void translate(Point& M, Point& A, Number& x_distance, Number& y_distance, Number& z_distance, std::vector<Polinom>& construction_polinoms);

void point_on_sphere(Point& A, Sphere& s, std::vector<Polinom>& polinoms);

void line_plane_intersection(Point& M, Line& l, Plane& pi, std::vector<Polinom>& polinoms);
void line_plane_intersection(Point& M, Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms);
void line_plane_intersection(Point& M, Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms);
void line_plane_intersection(Point& M, Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms);

void make_plane_orthogonal_on_plane(Plane& pi, Plane& tetha, Line& l, std::vector<Polinom>& construction_polinoms);
void make_plane_orthogonal_on_plane(Plane& pi, Plane& tetha, Point& A, Point& B, std::vector<Polinom>& construction_polinoms);
void make_plane_orthogonal_on_plane(Plane& pi, Point& A, Line& l, std::vector<Polinom>& construction_polinoms);
void make_plane_orthogonal_on_plane(Plane& pi, Sphere& s, Point& A, std::vector<Polinom>& construction_polinoms);

void equal_angles(Point& A, Point& B, Point& C, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms);

void make_line_orthogonal_on_plane(Line& l, Point& A, Point&  P, Point& Q, Point& R, std::vector<Polinom>& construction_polinoms);

void collinear(Point& A, Point& B, Point& C, std::vector<Polinom>& polinoms);

void equal_numbers(Number& n1, Number& n2, std::vector<Polinom>& polinoms);
void multiply(Number& n1, int broj, Number& n2, std::vector<Polinom>& polinoms);
void multiply(Number& n1, Number& n2, Number& n3, std::vector<Polinom>& polinoms);
void division(Number& n1, Number& n2, Number& n3, std::vector<Polinom>& polinoms);
void addition(Number& n1, Number& n2, Number& n3, std::vector<Polinom>& polinoms);

Polinom distance2(Point& A, Point& B, std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms);



#endif //STATEMENT_H