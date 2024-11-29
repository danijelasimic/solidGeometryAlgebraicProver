#include "number.h"

void Number::print() const 
{
    std::cout << name() << "= "; 
	
	if (_variable->has_value()) 
		 std::cout << _variable->value() << std::endl;
	else
	{
		std::cout << "x_" << _variable->variable_index() << std::endl;
	}
}

void Number::make_square_distance(Point& A, Point& B, std::vector<Polinom>& construction_polinoms)
{
	Polinom res = A.distance(B) - variable();

	construction_polinoms.push_back(res);
}

void Number::make_distance(Point& A, Point& B, std::vector<Polinom>& construction_polinoms)
{
	Polinom res = A.distance(B) - (Polinom)variable() * variable();

	construction_polinoms.push_back(res);
}

/* distance from point M to line determined with AB. */
/* d = |AB x MA|/|AB| */
void Number::make_distance_line(Point& M, Point& A, Point& B, std::vector<Polinom>& construction_polinoms,
	std::vector<std::shared_ptr<Variable>>& variables)
{
	auto distAB = std::make_shared<Variable>();
	auto dist2  = std::make_shared<Variable>();
	
	auto vx  = std::make_shared<Variable>();
	auto vy  = std::make_shared<Variable>();
	auto vz  = std::make_shared<Variable>();


	variables.insert(variables.end(), {distAB, dist2, vx, vy, vz});

	auto vecP = vector_product_of_vectors(A, B, M, A);

	Polinom p1 = distAB - A.distance(B);
	Polinom p2 = dist2 - (Polinom)vx*vx - (Polinom)vy*vy - (Polinom)vz*vz;
	Polinom p3 = vx - vecP[0];
	Polinom p4 = vy - vecP[1];
	Polinom p5 = vz - vecP[2];
	Polinom p6 = (Polinom)variable()*variable()*distAB - dist2;
	
	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3, p4, p5, p6});
}