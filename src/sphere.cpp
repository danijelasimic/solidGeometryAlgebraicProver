#include "sphere.h"

void Sphere::print() const
{
	std::cout << name() << ":\n" << "Center(" 
								 << "x_" << _centerX->variable_index() << ", "
                            	 << "x_" << _centerY->variable_index() << ", "
                            	 << "x_" << _centerZ->variable_index() 
                    			 << ")" 
                    		     << std::endl                
								 << "Radius = ";
	
	if (_radius->has_value()) 
		 std::cout << _radius->value() << std::endl;
	else
	{
		std::cout << "x_" << _radius->variable_index() << std::endl;
	}
}

void Sphere::set_center_and_radius(Point& center, Number& num)
{
	_centerX = center.x();
	_centerY = center.y();
	_centerZ = center.z();

	_radius = num.variable();
}

/* This is not the best design decision, but it is simple.
   Center coordinates should be the same as for coordinates given by variables
   _centerX, _centerY, _centerZ.
   This is done while calling contruction of spehere.
*/
void Sphere::make_sphere_through_four_points(Point& A, Point& B, Point& C, Point& D, Point& Center, std::vector<Polinom>& construction_polinoms)
{
	Polinom square_radius = Polinom{Term{-1, _radius, 2}};

	Polinom dist1 = A.distance(Center) + square_radius;
	Polinom dist2 = B.distance(Center) + square_radius;
	Polinom dist3 = C.distance(Center) + square_radius;
	Polinom dist4 = D.distance(Center) + square_radius;

	construction_polinoms.insert(construction_polinoms.end(), {dist1, dist2, dist3, dist4});
}