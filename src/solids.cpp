#include "solids.h"

void make_regular_tetrahedron(Point& C, Point& D, std::vector<Polinom>& construction_polinoms)
{
	C.set_point_coordinate_value(3);

	// 2cx - 1
	Polinom p1{Term{2, C.x()}, Term{-1}};
	construction_polinoms.push_back(p1);

	// 4cy^2 - 3
	Polinom p2{Term{4, C.y(), 2}, Term{-3}};
	construction_polinoms.push_back(p2);

	//3dy - cy
	Polinom p3{Term{3, D.y()}, Term{-1, C.y()}};
	construction_polinoms.push_back(p3);

	//dx - cx - same first coordinate
	Polinom p4{Term{1, D.x()}, Term{-1, C.x()}};
	construction_polinoms.push_back(p4);

	//3dz^2 - 2
	Polinom p5{Term{3, D.z(), 2}, Term{-2}};
	construction_polinoms.push_back(p5);
}

void make_tetrahedron(Point& C, Point& D, std::vector<Polinom>& construction_polinoms, 
					  std::vector<std::shared_ptr<Variable>>& variables)
{
	C.set_point_coordinate_value(3);

	/* Adding condition that third D coordinate is not equal zero. */
	auto zz = std::make_shared<Variable>();
	variables.push_back(zz);

	Polinom p = (Polinom)D.z() * zz - 1;
	construction_polinoms.push_back(p);
}

// Documentation for explanation of equations
void make_regular_hexagon(Point& A3, Point& A4, Point& A5, Point& A6, std::vector<Polinom>& construction_polinoms)
{
	A3.set_point_coordinate_value(3);
	A4.set_point_coordinate_value(1, 1);
	A4.set_point_coordinate_value(3);
	A5.set_point_coordinate_value(1);
	A5.set_point_coordinate_value(3);
	A6.set_point_coordinate_value(6);

	
	Polinom p1{Term{2, A3.x()}, Term{-3}};
	Polinom p2{Term{4, A3.y(), 2}, Term{-3}};
	Polinom p3{Term{1, A4.y()}, Term{-3}};
	Polinom p4{Term{2, A6.x()}, Term{-1}};

	// Some coordinates are same, but I use polinomials to express that.
	// Can this be improved?
	Polinom p5{Term{1, A5.y()}, Term{-1, A4.y()}};
	Polinom p6{Term{1, A6.y()}, Term{-1, A3.y()}};

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3, p4, p5, p6});
}

void make_regular_triangle(Point& C, std::vector<Polinom>& construction_polinoms)
{
	C.set_point_coordinate_value(3);

	// 2cx - 1
	Polinom p1{Term{2, C.x()}, Term{-1}};
	construction_polinoms.push_back(p1);

	// 4cy^2 - 3
	Polinom p2{Term{4, C.y(), 2}, Term{-3}};
	construction_polinoms.push_back(p2);
}

void make_regular_octagon(Point& A2, Point& A3, Point& A4, Point& A5, Point& A6, Point& A7,
                  std::vector<Polinom>& construction_polinoms)
{

}				  
                  
      
void make_top_vertex(Point& Top, int number_base, std::vector<Polinom>& construction_polinoms)
{
	switch (number_base) {
		case 3: {
			Polinom p1{Term{2, Top.x()}, Term{-1}};
			construction_polinoms.push_back(p1);

			Polinom p2{Term{12, Top.y(), 2}, Term{-1}};
			construction_polinoms.push_back(p2);
			break;
		}
		case 4: {
			Polinom p1{Term{2, Top.x()}, Term{-1}};
			construction_polinoms.push_back(p1);

			Polinom p2{Term{2, Top.y()}, Term{-1}};
			construction_polinoms.push_back(p2);
			break;
		}
		case 6: {
			Polinom p1{Term{2, Top.x()}, Term{-1}};
			construction_polinoms.push_back(p1);

			Polinom p2{Term{4, Top.y(), 2}, Term{-3}};
			construction_polinoms.push_back(p2);
			break;
		}
		case 8: {
			break;
		}
		default: return; // do ovog slucaja ne bi smelo da dodje
	}
}

void make_parallelogram(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& construction_polinoms)
{
	Polinom p1 = (Polinom)C.x() - D.x() - B.x() + A.x();
	Polinom p2 = (Polinom)C.y() - D.y() - B.y() + A.y();
	Polinom p3 = (Polinom)C.z() - D.z() - B.z() + A.z();

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3});
}