#include "line.h"

void Line::print() const
{
	std::cout << name() << "(" << "x_" << _vector_x->variable_index() << ", "
                            	<< "x_" << _vector_y->variable_index() << ", "
                            	<< "x_" << _vector_z->variable_index() 
                    	<< ")"
                    	<< std::endl;                  

}

void Line::make_line_trough_points(Point& A, Point& B, std::vector<Polinom>& construction_polinoms)
{
	/* Any of the two given points belong to the line, but we choose A to be in equation. */
	/* We just copy A coordinates to the line coordinates. */
	_point_x = A.x();
	_point_y = A.y();
	_point_z = A.z();

	//_vx - a._x + b._x
	Polinom p1{Term{1, _vector_x}, Term{-1, A.x()}, Term{1, B.x()}};
	construction_polinoms.push_back(p1);

	//_vy - a._y + b._y
	Polinom p2{Term{1, _vector_y}, Term{-1, A.y()}, Term{1, B.y()}};
	construction_polinoms.push_back(p2);

	//_vz - a._z + b._z
	Polinom p3{Term{1, _vector_z}, Term{-1, A.z()}, Term{1, B.z()}};
	construction_polinoms.push_back(p3);	
}

void Line::set_zeros_for_coordinates()
{
	_vector_x->set_value(0);
	_vector_y->set_value(0);

	_point_x->set_value(0);
	_point_y->set_value(0);
	_point_z->set_value(0);
}

void Line::set_point_coordinates(const std::shared_ptr<Variable>& px, const std::shared_ptr<Variable>& py,
								const std::shared_ptr<Variable>& pz)
{
	_point_x = px;
	_point_y = py;
	_point_z = pz;
}								