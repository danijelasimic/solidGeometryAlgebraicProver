#ifndef NDGCONDITIONS_H
#define NDGCONDITIONS_H

#include <vector>
#include <string>

#include "line.h"
#include "point.h"
#include "polinom.h"
#include "plane.h"
#include "number.h"
#include "sphere.h"

class NDGCondition 
{
public:
    /* Carefful, doesn't do anything! */
	NDGCondition() = default;

    NDGCondition(const std::string& ndgText, const Polinom& ndgPolinoms)
    : _ndgText {ndgText}, _ndgPolinoms {std::vector<Polinom>{ndgPolinoms}}
    {}


    NDGCondition(const std::string& ndgText, std::initializer_list<Polinom> ndgPolinoms)
    : _ndgText {ndgText}, _ndgPolinoms {ndgPolinoms}
    {}

    NDGCondition(const std::string& ndgText, std::vector<Polinom>& ndgPolinoms)
    : _ndgText {ndgText}, _ndgPolinoms {ndgPolinoms}
    {}


    ~NDGCondition() = default;

    void print(std::ofstream& out) const;

    //efficient, the caller cannot change text
    //but no new copy is made
    const std::string& getNdgText() const { return _ndgText; }

    const std::vector<Polinom>& getNdgPolinoms() const 
    {
        return _ndgPolinoms;
    }

private:
    std::string _ndgText;
    std::vector<Polinom> _ndgPolinoms;
};

/* Different point construction NDG */
void insert_ndg_points_not_equal(Point& A, Point& B, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_number_not_zero(Number& m, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_line_vector_not_zero(Line& l, std::vector<NDGCondition>& ndg_polinoms);

void insert_ndg_lines_not_parallel(Line& l, Line& g, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_lines_not_skew(Line& l, Line& g, std::vector<NDGCondition>& ndg_polinoms);

void insert_ndg_lines_not_parallel(Point& A, Point& B, Point& C, Point& D, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_lines_not_skew(Point& A, Point& B, Point& C, Point& D, std::vector<NDGCondition>& ndg_polinoms);

void insert_ndg_plane_vector_not_zero(Plane& pi, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_plane_and_line_not_parallel(Line& l, Plane& pi, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_plane_and_line_not_parallel(Point& A, Point& B, Plane& pi, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_plane_and_line_not_parallel(Line& l, Point& A, Point& B, Point& C, std::vector<NDGCondition>& ndg_polinoms);
void insert_ndg_plane_and_line_not_parallel(Point& P, Point& Q,
                                Point& A, Point& B, Point& C, std::vector<NDGCondition>& ndg_polinoms);

void insert_ndg_point_not_on_line(Point& A, Line& l, std::vector<NDGCondition>& ndg_polinoms);  
void insert_ndg_lines_with_same_vector_not_equal(Line& l1, Line& l2, std::vector<NDGCondition>& ndg_polinoms);   

//helper function for insert_ndg_lines_with_same_vector_not_equal
std::vector<Polinom> point_not_on_line(Point& A, Polinom& x_vec, Polinom& y_vec, Polinom& z_vec, Point linePoint);         
void insert_ndg_lines_with_same_vector_not_equal(Point& A, Point& B, Point& C, Point& D, std::vector<NDGCondition>& ndg_polinoms); 

void insert_ndg_point_line_distance_not_zero(Number& n, Point& M, Point& A, Point& B, std::vector<NDGCondition>& ndg_polinoms); 

#endif //NDGCONDITIONS_H