#include "point.h"


void Point::print() const
{
	std::cout << name() << "(" << "x_" << coordinate_x->variable_index() << ", "
	                           << "x_" << coordinate_y->variable_index() << ", "
	                           << "x_" << coordinate_z->variable_index()
	                    << ")"
	                    << std::endl;                  
}

std::shared_ptr<Variable> Point::coordinate(int i) const
{
	switch(i)
	{
		case 1: return coordinate_x;
		case 2: return coordinate_y;
		case 3: return coordinate_z;
		default: return nullptr;
	}
}

std::shared_ptr<Variable> Point::x() const
{
	return coordinate_x;
}

std::shared_ptr<Variable> Point::y() const
{
	return coordinate_y;
}

std::shared_ptr<Variable> Point::z() const
{
	return coordinate_z;
}


void Point::set_point_coordinate_value(int coordinate, int value)
{
	switch(coordinate)
	{
		case 1: coordinate_x->set_value(value); break;
		case 2: coordinate_y->set_value(value); break;
		case 3: coordinate_z->set_value(value); break;
		default: return; //TODO throw error
	}
}

void Point::set_all_coordinates_values(int val1, int val2, int val3)
{
	coordinate_x->set_value(val1);
	coordinate_y->set_value(val2);
	coordinate_z->set_value(val3);
}


void Point::make_midpoint(Point& A, Point& B, std::vector<Polinom>& construction_polinoms)
{
	//2mx - ax -bx
	Polinom p1{Term{2, x()}, Term{-1, A.x()}, Term{-1, B.x()}};
	construction_polinoms.push_back(p1);

	//2my - ay -by
	Polinom p2{Term{2, y()}, Term{-1, A.y()}, Term{-1, B.y()}};
	construction_polinoms.push_back(p2);

	//2mx - az -bz
	Polinom p3{Term{2, z()}, Term{-1, A.z()}, Term{-1, B.z()}};
	construction_polinoms.push_back(p3);
}


Polinom Point::distance(Point& B)
{
	Polinom p{Term{1, x(), 2}, Term{-2, {x(), B.x()}}, Term{1, B.x(), 2},
			  Term{1, y(), 2}, Term{-2, {y(), B.y()}}, Term{1, B.y(), 2},
			  Term{1, z(), 2}, Term{-2, {z(), B.z()}}, Term{1, B.z(), 2}};

	return p;
}

void Point::translate(Point& A, int x_distance, int y_distance, int z_distance, std::vector<Polinom>& construction_polinoms)
{
	Polinom p1{Term{1, x()}, Term{-1, A.x()}, Term{-x_distance}};
	Polinom p2{Term{1, y()}, Term{-1, A.y()}, Term{-y_distance}};
	Polinom p3{Term{1, z()}, Term{-1, A.z()}, Term{-z_distance}};

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3});
}

Polinom scalar_product_of_vectors(const Point& A, const Point& B, const Point& C, const Point& D)
{

	Polinom rez = scalar_product_of_vectors((Polinom)B.x() - A.x(), (Polinom)B.y() - A.y(), (Polinom)B.z() - A.z(),
											(Polinom)D.x() - C.x(), (Polinom)D.y() - C.y(), (Polinom)D.z() - C.z());

	return rez;											
}

std::vector<Polinom> vector_product_of_vectors(const Point& A, const Point& B, const Point& C, const Point& D)
{
	return vector_product_of_vectors((Polinom)B.x() - A.x(), (Polinom)B.y() - A.y(), (Polinom)B.z() - A.z(),
											(Polinom)D.x() - C.x(), (Polinom)D.y() - C.y(), (Polinom)D.z() - C.z());
}






















