#include "plane.h"

void Plane::print() const
{
	std::cout << name() << "(" << "x_" << _vector_x->variable_index() << ", "
                            	<< "x_" << _vector_y->variable_index() << ", "
                            	<< "x_" << _vector_z->variable_index() << ", "
                                << "x_" << _vector_d->variable_index() 
                    	<< ")"
                    	<< std::endl;                  

}

void Plane::make_plane_trough_points(Point& A, Point& B, Point& C, std::vector<Polinom>& construction_polinoms)
{
	Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cax = Polinom{Term{1, C.x()}, Term{-1, A.x()}};
	Polinom vec_cay = Polinom{Term{1, C.y()}, Term{-1, A.y()}};
	Polinom vec_caz = Polinom{Term{1, C.z()}, Term{-1, A.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cax, vec_cay, vec_caz);	

    Polinom vecx = rez[0] - Polinom{Term{1, _vector_x}};
    Polinom vecy = rez[1] - Polinom{Term{1, _vector_y}};
    Polinom vecz = rez[2] - Polinom{Term{1, _vector_z}};

    construction_polinoms.insert(construction_polinoms.end(), {vecx, vecy, vecz});

    // Point A on plane
    Polinom incidentA = Polinom{Term{1, {_vector_x, A.x()}}, Term{1, {_vector_y, A.y()}},
                                Term{1, {_vector_z, A.z()}}, Term{1, _vector_d}};

    construction_polinoms.push_back(incidentA);                                
}
