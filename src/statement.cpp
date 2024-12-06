#include "statement.h"

void point_on_line(Point& A, Line& l, std::vector<Polinom>& statement_polinoms)
{
	//A (a1, a2, a3)
	//l (v1, v2, v3) i point_determining_line(p1, p2, p3)
	//(a1 - p1)*v2 - (a2 - p2)*v1 = 0 --> ax*vy - px*vy - ax*vx + py*vx
	//(a1 - p1)*v3 - (a3 - p3)*v1 = 0 --> ax*vz - px*vz - az*vx + pz*vx
	//(a2 - p2)*v3 - (a3 - p3)*v2 = 0 --> ay*vz - py*vz - az*vy + pz*vy

	Polinom p1{Term{1, {A.x(), l.y_vec()}}, 
	           Term{-1, {l.x_p(), l.y_vec()}},
	           Term{-1, {A.y(), l.x_vec()}},
	           Term{1, {l.y_p(), l.x_vec()}}};
	statement_polinoms.push_back(p1);

	Polinom p2{Term{1, {A.x(), l.z_vec()}}, Term{-1, {l.x_p(), l.z_vec()}},
			   Term{-1, {A.z(), l.x_vec()}}, Term{1, {l.z_p(), l.x_vec()}}};
	statement_polinoms.push_back(p2);

	Polinom p3{Term{1, {A.y(), l.z_vec()}}, Term{-1, {l.y_p(), l.z_vec()}},
				Term{-1, {A.z(), l.y_vec()}}, Term{1, {l.z_p(), l.y_vec()}}};
	statement_polinoms.push_back(p3);
}


void point_on_line(Point& A, Point& P, Point& Q, std::vector<Polinom>& statement_polinoms)
{
	//ay - py
	Polinom aypy{Term{1, A.y()}, Term{-1, P.y()}};
	//qz - pz
	Polinom qzpz{Term{1, Q.z()}, Term{-1, P.z()}};
	//az - pz
	Polinom azpz{Term{1, A.z()}, Term{-1, P.z()}};
	//qy - py
	Polinom qypy{Term{1, Q.y()}, Term{-1, P.y()}};
	//qx - px
	Polinom qxpx{Term{1, Q.x()}, Term{-1, P.x()}};
	//ax - px
	Polinom axpx{Term{1, A.x()}, Term{-1, P.x()}};

	Polinom statement1 = aypy * qzpz - azpz * qypy;
	Polinom statement2 = azpz * qxpx - axpx * qzpz;
	Polinom statement3 = axpx * qypy - aypy * qxpx;


	statement_polinoms.push_back(statement1);
	statement_polinoms.push_back(statement2);
	statement_polinoms.push_back(statement3);
}

void equal_points(Point& A, Point& P, std::vector<Polinom>& statement_polinoms)
{
	//ax - px
	Polinom axpx{Term{1, A.x()}, Term{-1, P.x()}};
	//ay - py
	Polinom aypy{Term{1, A.y()}, Term{-1, P.y()}};
	//az - pz
	Polinom azpz{Term{1, A.z()}, Term{-1, P.z()}};

	statement_polinoms.push_back(axpx);
	statement_polinoms.push_back(aypy);
	statement_polinoms.push_back(azpz);
}


void midpoint(Point& M, Point& A, Point& B, std::vector<Polinom>& statement_polinoms)
{
	//2mx - ax -bx
	Polinom p1{Term{2, M.x()}, Term{-1, A.x()}, Term{-1, B.x()}};
	statement_polinoms.push_back(p1);

	//2my - ay -by
	Polinom p2{Term{2, M.y()}, Term{-1, A.y()}, Term{-1, B.y()}};
	statement_polinoms.push_back(p2);

	//2mz - az -bz
	Polinom p3{Term{2, M.z()}, Term{-1, A.z()}, Term{-1, B.z()}};
	statement_polinoms.push_back(p3);
}

void point_segment_ratio(Point& M, Point& A, Point& B, int m, int n, std::vector<Polinom>& statement_polinoms)
{
	// U ovoj varijanti je M izmedju A i B, tj. |MA|/|MB| = m/n
	if (m > 0) {
		//(m+n)mx - n*ax - m*bx
		Polinom p1{Term{m + n, M.x()}, Term{-n, A.x()}, Term{-m, B.x()}};
		statement_polinoms.push_back(p1);

		//(m+n)my - n*ay - m*by
		Polinom p2{Term{m + n, M.y()}, Term{-n, A.y()}, Term{-m, B.y()}};
		statement_polinoms.push_back(p2);

		//(m+n)mz - n*az - m*bz
		Polinom p3{Term{m + n, M.z()}, Term{-n, A.z()}, Term{-m, B.z()}};
		statement_polinoms.push_back(p3);	
	} else {
		// u ovoj varijanti je A izmedju M i B, tj. |AM|/|AB| = m/(n - m) (i dalje vazi |MA|/|MB| = m/n)
		m = -m;
		n = n - m;
		//(m+n)mx - n*ax - m*bx
		Polinom p1{Term{m + n, A.x()}, Term{-n, M.x()}, Term{-m, B.x()}};
		statement_polinoms.push_back(p1);

		//(m+n)my - n*ay - m*by
		Polinom p2{Term{m + n, A.y()}, Term{-n, M.y()}, Term{-m, B.y()}};
		statement_polinoms.push_back(p2);

		//(m+n)mz - n*az - m*bz
		Polinom p3{Term{m + n, A.z()}, Term{-n, M.z()}, Term{-m, B.z()}};
		statement_polinoms.push_back(p3);
	}
}

void point_segment_ratio(Point& M, Point& A, Point& B, Number& m, Number& n, std::vector<Polinom>& polinoms)
{
	//Polinom p1{Term{m + n, M.x()}, Term{-n, A.x()}, Term{-m, B.x()}};
	Polinom p1 = (Polinom)m.variable()*M.x() + (Polinom)n.variable()*M.x() - (Polinom)n.variable()*A.x() - (Polinom)m.variable()*B.x();
	polinoms.push_back(p1);

	//(m+n)my - n*ay - m*by
	//Polinom p2{Term{m + n, M.y()}, Term{-n, A.y()}, Term{-m, B.y()}};
	Polinom p2 = (Polinom)m.variable()*M.y() + (Polinom)n.variable()*M.y() - (Polinom)n.variable()*A.y() - (Polinom)m.variable()*B.y();
	polinoms.push_back(p2);

	//(m+n)mz - n*az - m*bz
	//Polinom p3{Term{m + n, M.z()}, Term{-n, A.z()}, Term{-m, B.z()}};
	Polinom p3 = (Polinom)m.variable()*M.z() + (Polinom)n.variable()*M.z() - (Polinom)n.variable()*A.z() - (Polinom)m.variable()*B.z();
	polinoms.push_back(p3);
}

Polinom distance2(Point& A, Point& B, std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms)
{
	auto x = std::make_shared<Variable>();
    auto y = std::make_shared<Variable>();
	auto z = std::make_shared<Variable>();

	variables.insert(variables.end(), {x, y, z});

	Polinom p1 = x - (Polinom)A.x() - B.x();
	Polinom p2 = y - (Polinom)A.y() - B.x();
	Polinom p3 = z - (Polinom)A.x() - B.x();

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3});

	Polinom rez = (Polinom)x*x + (Polinom)y*y + (Polinom)z*z;

	return rez;
}

void congruent(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& statement_polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms)
{
	// (a x − b x )^2 + (a y − b y )^2 + (a z − b z )^2 − (c x − d x )^2 − (c y − d y )^2 − (c z − d z )^2 = 0
	//evry line is: ax^2 - 2*ax*bx + bx^2
	// Polinom p{
	// 	Term{1, A.x(), 2}, Term{-2, {A.x(), B.x()}}, Term{1, B.x(), 2},
	// 	Term{1, A.y(), 2}, Term{-2, {A.y(), B.y()}}, Term{1, B.y(), 2},
	// 	Term{1, A.z(), 2}, Term{-2, {A.z(), B.z()}}, Term{1, B.z(), 2},
	// 	Term{-1, C.x(), 2}, Term{2, {C.x(), D.x()}}, Term{-1, D.x(), 2},
	// 	Term{-1, C.y(), 2}, Term{2, {C.y(), D.y()}}, Term{-1, D.y(), 2},
	// 	Term{-1, C.z(), 2}, Term{2, {C.z(), D.z()}}, Term{-1, D.z(), 2}
	// };

	// statement_polinoms.push_back(p);
	auto distAB = std::make_shared<Variable>();
	auto distCD = std::make_shared<Variable>();

    variables.insert(variables.end(), {distAB, distCD});

	//Polinom p1{distAB - A.distance(B)};
	//Polinom p2{distCD - C.distance(D)};
	Polinom p1{distAB - distance2(A, B, variables, construction_polinoms)};
	Polinom p2{distCD - distance2(C, D, variables, construction_polinoms)};

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2});

	Polinom rez = (Polinom)distAB - (Polinom)distCD;

	statement_polinoms.push_back(rez);
}

void segments_in_ratio(Point& A, Point& B, Point& C, Point& D, int m, int n, std::vector<Polinom>& statement_polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms)
{
	// n^2*(a x − b x )^2 + n^2*(a y − b y )^2 + n^2*(a z − b z )^2 − m^2*(c x − d x )^2 − m^2*(c y − d y )^2 − m^2*(c z − d z )^2 = 0
	//every line is: n^2*ax^2 - n^2*2*ax*bx + n^2*bx^2
	/*
	Polinom p{
		Term{n*n, A.x(), 2}, Term{-2*n*n, {A.x(), B.x()}}, Term{n*n, B.x(), 2},
		Term{n*n, A.y(), 2}, Term{-2*n*n, {A.y(), B.y()}}, Term{n*n, B.y(), 2},
		Term{n*n, A.z(), 2}, Term{-2*n*n, {A.z(), B.z()}}, Term{n*n, B.z(), 2},
		Term{-1*m*m, C.x(), 2}, Term{2*m*m, {C.x(), D.x()}}, Term{-1*m*m, D.x(), 2},
		Term{-1*m*m, C.y(), 2}, Term{2*m*m, {C.y(), D.y()}}, Term{-1*m*m, D.y(), 2},
		Term{-1*m*m, C.z(), 2}, Term{2*m*m, {C.z(), D.z()}}, Term{-1*m*m, D.z(), 2}
	};
	statement_polinoms.push_back(p);
	*/
	auto distAB = std::make_shared<Variable>();
	auto distCD = std::make_shared<Variable>();

    variables.insert(variables.end(), {distAB, distCD});

	Polinom p1{distAB - A.distance(B)};
	Polinom p2{distCD - C.distance(D)};

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2});

	Polinom rez = (Polinom)distAB * (Polinom)n * (Polinom)n - (Polinom)distCD * (Polinom)m * (Polinom)m;

	statement_polinoms.push_back(rez);
}

void line_intersection(Point& A, Line& l1, Line& l2, std::vector<Polinom>& polinoms)
{
	point_on_line(A, l1, polinoms);
	point_on_line(A, l2, polinoms);
}

void line_intersection(Point& A, Point& P, Point& Q, Point& R, Point& S, std::vector<Polinom>& polinoms)
{
	point_on_line(A, P, Q, polinoms);
	point_on_line(A, R, S, polinoms);
}

void orthogonal_lines(Line& l1, Line& l2, std::vector<Polinom>& polinoms)
{
	Polinom orth = Polinom{Term{1, {l1.x_vec(), l2.x_vec()}},
	                     	Term{1, {l1.y_vec(), l2.y_vec()}},
							Term{1, {l1.z_vec(), l2.z_vec()}}};
	polinoms.push_back(orth);

	//not_skew(l1, l2, polinoms);
}

void orthogonal_lines(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& polinoms)
{
	Polinom xx = Polinom{Term{1, {B.x(), D.x()}}, Term{-1, {B.x(), C.x()}},
						Term{-1, {A.x(), D.x()}}, Term{1, {A.x(), C.x()}}};

	Polinom yy = Polinom{Term{1, {B.y(), D.y()}}, Term{-1, {B.y(), C.y()}},
						Term{-1, {A.y(), D.y()}}, Term{1, {A.y(), C.y()}}};

	Polinom zz = Polinom{Term{1, {B.z(), D.z()}}, Term{-1, {B.z(), C.z()}},
						Term{-1, {A.z(), D.z()}}, Term{1, {A.z(), C.z()}}};												

	Polinom orth = xx + yy + zz;
	polinoms.push_back(orth);

	//not_skew(A, B, C, D, polinoms);
}

void not_skew(Line& l1, Line& l2, std::vector<Polinom>& polinoms)
{
	Polinom vec_abx = Polinom{Term{1, l1.x_vec()}};
	Polinom vec_aby = Polinom{Term{1, l1.y_vec()}};
	Polinom vec_abz = Polinom{Term{1, l1.z_vec()}};

	Polinom vec_acx = Polinom{Term{1, l1.x_p()}, Term{-1, l2.x_p()}};
	Polinom vec_acy = Polinom{Term{1, l1.y_p()}, Term{-1, l2.y_p()}};
	Polinom vec_acz = Polinom{Term{1, l1.z_p()}, Term{-1, l2.z_p()}};

	Polinom vec_cdx = Polinom{Term{1, l2.x_vec()}};
	Polinom vec_cdy = Polinom{Term{1, l2.y_vec()}};
	Polinom vec_cdz = Polinom{Term{1, l2.z_vec()}};

	// mix product of AC AB CD
	polinoms.push_back(
		mix_product_of_vectors(vec_acx, vec_acy, vec_acz,
								vec_abx, vec_aby, vec_abz,
								vec_cdx, vec_cdy, vec_cdz));
}

void not_skew(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& polinoms)
{
	Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_acx = Polinom{Term{1, C.x()}, Term{-1, A.x()}};
	Polinom vec_acy = Polinom{Term{1, C.y()}, Term{-1, A.y()}};
	Polinom vec_acz = Polinom{Term{1, C.z()}, Term{-1, A.z()}};

	Polinom vec_cdx = Polinom{Term{1, D.x()}, Term{-1, C.x()}};
	Polinom vec_cdy = Polinom{Term{1, D.y()}, Term{-1, C.y()}};
	Polinom vec_cdz = Polinom{Term{1, D.z()}, Term{-1, C.z()}};

	// mix product of AC AB CD
	polinoms.push_back(
		mix_product_of_vectors(vec_acx, vec_acy, vec_acz,
								vec_abx, vec_aby, vec_abz,
								vec_cdx, vec_cdy, vec_cdz));
}


void parallel_lines(Line& l1, Line& l2, std::vector<Polinom>& polinoms)
{
	Polinom vec_abx = Polinom{Term{1, l1.x_vec()}};
	Polinom vec_aby = Polinom{Term{1, l1.y_vec()}};
	Polinom vec_abz = Polinom{Term{1, l1.z_vec()}};

	Polinom vec_cdx = Polinom{Term{1, l2.x_vec()}};
	Polinom vec_cdy = Polinom{Term{1, l2.y_vec()}};
	Polinom vec_cdz = Polinom{Term{1, l2.z_vec()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cdx, vec_cdy, vec_cdz);
	polinoms.insert(polinoms.end(), rez.begin(), rez.end());					
}

void parallel_lines(Point& A, Point& B, Point& C, Point& D, std::vector<Polinom>& polinoms)
{
	Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cdx = Polinom{Term{1, D.x()}, Term{-1, C.x()}};
	Polinom vec_cdy = Polinom{Term{1, D.y()}, Term{-1, C.y()}};
	Polinom vec_cdz = Polinom{Term{1, D.z()}, Term{-1, C.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cdx, vec_cdy, vec_cdz);
	polinoms.insert(polinoms.end(), rez.begin(), rez.end());														
}

void make_plane_trough_points(Plane& pi, Point& A, Point& B, Point& C, std::vector<Polinom>& construction_polinoms)
{
	Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cax = Polinom{Term{1, C.x()}, Term{-1, A.x()}};
	Polinom vec_cay = Polinom{Term{1, C.y()}, Term{-1, A.y()}};
	Polinom vec_caz = Polinom{Term{1, C.z()}, Term{-1, A.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cax, vec_cay, vec_caz);	

    Polinom vecx = rez[0] - Polinom{Term{1, pi.x_vec()}};
    Polinom vecy = rez[1] - Polinom{Term{1, pi.y_vec()}};
    Polinom vecz = rez[2] - Polinom{Term{1, pi.z_vec()}};

    construction_polinoms.insert(construction_polinoms.end(), {vecx, vecy, vecz});

    // Point A on plane
    Polinom incidentA = Polinom{Term{1, {pi.x_vec(), A.x()}}, Term{1, {pi.y_vec(), A.y()}},
                                Term{1, {pi.z_vec(), A.z()}}, Term{1, pi.d_vec()}};

    construction_polinoms.push_back(incidentA);
}

void make_vector_of_plane_trough_points(Plane& pi, Point& A, Point& B, Point& C, std::vector<Polinom>& construction_polinoms)
{
	Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cax = Polinom{Term{1, C.x()}, Term{-1, A.x()}};
	Polinom vec_cay = Polinom{Term{1, C.y()}, Term{-1, A.y()}};
	Polinom vec_caz = Polinom{Term{1, C.z()}, Term{-1, A.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cax, vec_cay, vec_caz);

    Polinom vecx = rez[0] - Polinom{Term{1, pi.x_vec()}};
    Polinom vecy = rez[1] - Polinom{Term{1, pi.y_vec()}};
    Polinom vecz = rez[2] - Polinom{Term{1, pi.z_vec()}};

    construction_polinoms.insert(construction_polinoms.end(), {vecx, vecy, vecz});
}

void point_in_plane(Point& A, Plane& pi, std::vector<Polinom>& polinoms)
{
	Polinom incidentA = Polinom{Term{1, {pi.x_vec(), A.x()}}, Term{1, {pi.y_vec(), A.y()}},
                                Term{1, {pi.z_vec(), A.z()}}, Term{1, pi.d_vec()}};

    polinoms.push_back(incidentA);
}

void point_in_plane(Point& A, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms)
{
	Polinom vec_pax = Polinom{Term{1, A.x()}, Term{-1, P.x()}};
	Polinom vec_pay = Polinom{Term{1, A.y()}, Term{-1, P.y()}};
	Polinom vec_paz = Polinom{Term{1, A.z()}, Term{-1, P.z()}};

	Polinom vec_pqx = Polinom{Term{1, Q.x()}, Term{-1, P.x()}};
	Polinom vec_pqy = Polinom{Term{1, Q.y()}, Term{-1, P.y()}};
	Polinom vec_pqz = Polinom{Term{1, Q.z()}, Term{-1, P.z()}};

	Polinom vec_prx = Polinom{Term{1, R.x()}, Term{-1, P.x()}};
	Polinom vec_pry = Polinom{Term{1, R.y()}, Term{-1, P.y()}};
	Polinom vec_prz = Polinom{Term{1, R.z()}, Term{-1, P.z()}};

	Polinom incident = mix_product_of_vectors(vec_pax, vec_pay, vec_paz, vec_pqx, vec_pqy,
												vec_pqz, vec_prx, vec_pry, vec_prz);

	polinoms.push_back(incident);												
}


void parallel_planes(Plane& pi, Plane& teta, std::vector<Polinom>& polinoms)
{
	Polinom vec_pix = Polinom{Term{1, pi.x_vec()}};
	Polinom vec_piy = Polinom{Term{1, pi.y_vec()}};
	Polinom vec_piz = Polinom{Term{1, pi.z_vec()}};

	Polinom vec_tetax = Polinom{Term{1, teta.x_vec()}};
	Polinom vec_tetay = Polinom{Term{1, teta.y_vec()}};
	Polinom vec_tetaz = Polinom{Term{1, teta.z_vec()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_pix, vec_piy, vec_piz,
															vec_tetax, vec_tetay, vec_tetaz);

	polinoms.insert(polinoms.end(), {rez[0], rez[1], rez[2]});							
}

/* In order to simplify, first make plaines, and then check if statement holds.
*/
void parallel_planes(Point& A, Point& B, Point& C, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane pi(name, variables);
	Plane teta(name, variables);

	make_vector_of_plane_trough_points(pi, A, B, C, construction_polinoms);
	make_vector_of_plane_trough_points(teta, P, Q, R, construction_polinoms);

	parallel_planes(pi, teta, polinoms);
}

void parallel_planes(Plane& pi, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane teta(name, variables);
	make_vector_of_plane_trough_points(teta, P, Q, R, construction_polinoms);

	parallel_planes(pi, teta, polinoms);
}						

void orthogonal_planes(Plane& pi, Plane& teta, std::vector<Polinom>& polinoms)
{
	Polinom rez = Polinom{
		Term{1, {pi.x_vec(), teta.x_vec()}},
		Term{1, {pi.y_vec(), teta.y_vec()}},
		Term{1, {pi.z_vec(), teta.z_vec()}}};
	polinoms.push_back(rez);
}

void orthogonal_planes(Point& A, Point& B, Point& C, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane pi(name, variables);
	Plane teta(name, variables);

	make_vector_of_plane_trough_points(pi, A, B, C, construction_polinoms);
	make_vector_of_plane_trough_points(teta, P, Q, R, construction_polinoms);

	orthogonal_planes(pi, teta, polinoms);
}		

void orthogonal_planes(Plane& pi, Point& P, Point& Q, Point& R, 
						std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
						std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane teta(name, variables);
	make_vector_of_plane_trough_points(teta, P, Q, R, construction_polinoms);

	orthogonal_planes(pi, teta, polinoms);
}

void parallel_line_plane(Line& l, Plane& pi, std::vector<Polinom>& polinoms)
{
	Polinom rez = Polinom{
		Term{1, {pi.x_vec(), l.x_vec()}},
		Term{1, {pi.y_vec(), l.y_vec()}},
		Term{1, {pi.z_vec(), l.z_vec()}}};
	polinoms.push_back(rez);
}

void parallel_line_plane(Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms)
{
	Polinom rez = Polinom{
		Term{1, {pi.x_vec(), B.x()}}, Term{-1, {pi.x_vec(), A.x()}},
		Term{1, {pi.y_vec(), B.y()}}, Term{-1, {pi.y_vec(), A.y()}},
		Term{1, {pi.z_vec(), B.z()}}, Term{-1, {pi.z_vec(), A.z()}}};
	polinoms.push_back(rez);
}

void parallel_line_plane(Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms,
							std::vector<std::shared_ptr<Variable>>& variables,
							std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane pi(name, variables);
	make_vector_of_plane_trough_points(pi, P, Q, R, construction_polinoms);

	parallel_line_plane(l, pi, polinoms);
}

void parallel_line_plane(Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms, std::vector<std::shared_ptr<Variable>>& variables,
							std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane pi(name, variables);
	make_vector_of_plane_trough_points(pi, P, Q, R, construction_polinoms);

	parallel_line_plane(A, B, pi, polinoms);
}

void orthogonal_line_plane(Line& l, Plane& pi, std::vector<Polinom>& polinoms)
{
	Polinom vec_pix = Polinom{Term{1, pi.x_vec()}};
	Polinom vec_piy = Polinom{Term{1, pi.y_vec()}};
	Polinom vec_piz = Polinom{Term{1, pi.z_vec()}};

	Polinom vec_lx = Polinom{Term{1, l.x_vec()}};
	Polinom vec_ly = Polinom{Term{1, l.y_vec()}};
	Polinom vec_lz = Polinom{Term{1, l.z_vec()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_pix, vec_piy, vec_piz,
															vec_lx, vec_ly, vec_lz);

	polinoms.insert(polinoms.end(), {rez[0], rez[1], rez[2]});
}

void orthogonal_line_plane(Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms)
{
	Polinom vec_pix = Polinom{Term{1, pi.x_vec()}};
	Polinom vec_piy = Polinom{Term{1, pi.y_vec()}};
	Polinom vec_piz = Polinom{Term{1, pi.z_vec()}};

	Polinom vec_lx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_ly = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_lz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_pix, vec_piy, vec_piz,
															vec_lx, vec_ly, vec_lz);

	polinoms.insert(polinoms.end(), {rez[0], rez[1], rez[2]});
}

void orthogonal_line_plane(Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms)
{
	Polinom vec_pqx = Polinom{Term{1, Q.x()}, Term{-1, P.x()}};
	Polinom vec_pqy = Polinom{Term{1, Q.y()}, Term{-1, P.y()}};
	Polinom vec_pqz = Polinom{Term{1, Q.z()}, Term{-1, P.z()}};

	Polinom vec_prx = Polinom{Term{1, R.x()}, Term{-1, P.x()}};
	Polinom vec_pry = Polinom{Term{1, R.y()}, Term{-1, P.y()}};
	Polinom vec_prz = Polinom{Term{1, R.z()}, Term{-1, P.z()}};

	Polinom vec_lx = Polinom{Term{1, l.x_vec()}};
	Polinom vec_ly = Polinom{Term{1, l.y_vec()}};
	Polinom vec_lz = Polinom{Term{1, l.z_vec()}};

	Polinom rez1 = scalar_product_of_vectors(vec_pqx, vec_pqy, vec_pqz,
												vec_lx, vec_ly, vec_lz);

	Polinom rez2 = scalar_product_of_vectors(vec_prx, vec_pry, vec_prz,
												vec_lx, vec_ly, vec_lz);

	polinoms.insert(polinoms.end(), {rez1, rez2});
}

void orthogonal_line_plane(Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms)
{
	Polinom vec_pqx = Polinom{Term{1, Q.x()}, Term{-1, P.x()}};
	Polinom vec_pqy = Polinom{Term{1, Q.y()}, Term{-1, P.y()}};
	Polinom vec_pqz = Polinom{Term{1, Q.z()}, Term{-1, P.z()}};

	Polinom vec_prx = Polinom{Term{1, R.x()}, Term{-1, P.x()}};
	Polinom vec_pry = Polinom{Term{1, R.y()}, Term{-1, P.y()}};
	Polinom vec_prz = Polinom{Term{1, R.z()}, Term{-1, P.z()}};

	Polinom vec_lx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_ly = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_lz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom rez1 = scalar_product_of_vectors(vec_pqx, vec_pqy, vec_pqz,
												vec_lx, vec_ly, vec_lz);

	Polinom rez2 = scalar_product_of_vectors(vec_prx, vec_pry, vec_prz,
												vec_lx, vec_ly, vec_lz);

	polinoms.insert(polinoms.end(), {rez1, rez2});
}

void line_in_plane(Line& l, Plane& pi, std::vector<Polinom>& polinoms)
{
	Point point_on_lineL = l.point();
	point_in_plane(point_on_lineL, pi, polinoms);

	Polinom p2 = scalar_product_of_vectors(l.x_vec(), l.y_vec(), l.z_vec(),
											pi.x_vec(), pi.y_vec(), pi.z_vec());
	polinoms.push_back(p2);
}

void line_in_plane(Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms)
{
	point_in_plane(A, pi, polinoms);
	point_in_plane(B, pi, polinoms);
}

void line_in_plane(Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms,
					std::vector<std::shared_ptr<Variable>>& variables,
					std::vector<Polinom>& construction_polinoms)
{
	std::string name="noname";
	Plane teta(name, variables);
	make_vector_of_plane_trough_points(teta, P, Q, R, construction_polinoms);

	line_in_plane(l, teta, polinoms);
}

void line_in_plane(Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms)
{
	point_in_plane(A, P, Q, R, polinoms);
	point_in_plane(B, P, Q, R, polinoms);
}							


void translate(Point& M, Point& A, Number& x_distance, Number& y_distance, Number& z_distance, std::vector<Polinom>& construction_polinoms)
{
	Polinom p1{Term{1, M.x()}, Term{-1, A.x()}, Term{-1, x_distance.variable()}};
	Polinom p2{Term{1, M.y()}, Term{-1, A.y()}, Term{-1, y_distance.variable()}};
	Polinom p3{Term{1, M.z()}, Term{-1, A.z()}, Term{-1, z_distance.variable()}};

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3});
}

void point_on_sphere(Point& A, Sphere& s, std::vector<Polinom>& polinoms)
{
	Point O{s.centerX(), s.centerY(), s.centerZ()};	

	Polinom p1{Term{-1, s.radius(), 2}};
	Polinom res = A.distance(O) + p1;

	polinoms.push_back(res);
}

void line_plane_intersection(Point& M, Line& l, Plane& pi, std::vector<Polinom>& polinoms)
{
	point_in_plane(M, pi, polinoms);
	point_on_line(M, l, polinoms);
}

void line_plane_intersection(Point& M, Point& A, Point& B, Plane& pi, std::vector<Polinom>& polinoms)
{
	point_on_line(M, A, B, polinoms);
	point_in_plane(M, pi, polinoms);
}

void line_plane_intersection(Point& M, Line& l, Point& P, Point& Q, Point& R, std::vector<Polinom>& polinoms)
{
	point_on_line(M, l, polinoms);
	point_in_plane(M, P, Q, R, polinoms);
}	

void line_plane_intersection(Point& M, Point& A, Point& B, Point& P, Point& Q, Point& R,
                            std::vector<Polinom>& polinoms)
{
	point_on_line(M, A, B, polinoms);
	point_in_plane(M, P, Q, R, polinoms);
}	


void make_plane_orthogonal_on_plane(Plane& pi, Plane& tetha, Line& l, std::vector<Polinom>& construction_polinoms)
{
	orthogonal_planes(pi, tetha, construction_polinoms);
	line_in_plane(l, pi, construction_polinoms);
}

void make_plane_orthogonal_on_plane(Plane& pi, Plane& tetha, Point& A, Point& B, std::vector<Polinom>& construction_polinoms)
{
	orthogonal_planes(pi, tetha, construction_polinoms);
	point_in_plane(A, pi, construction_polinoms);
	point_in_plane(B, pi, construction_polinoms);
}

void equal_angles(Point& A, Point& O1, Point& B, Point& C, Point& O2, Point& D, std::vector<Polinom>& polinoms,
					std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& construction_polinoms)
{
	/*
	scalar 1 = AO1 · BO1
	scalar 2 = CO2 · DO2
	dist CO2 = |CO2| -- ovo je kvadratno rastojanje
	dist DO2 = |DO2|
	dist AO1 = |AO1|
	dist BO1 = |BO1|
	poly = scalar 1 ∗ scalar 1 ∗ dist CO2 ∗ dist DO2 − scalar 2 ∗ scalar 2 ∗ dist AO1 ∗ dist BO1	
	*/
	auto scalar1 = std::make_shared<Variable>();
    auto scalar2 = std::make_shared<Variable>();
	auto distAO1 = std::make_shared<Variable>();
	auto distBO1 = std::make_shared<Variable>();
	auto distCO2 = std::make_shared<Variable>();
	auto distDO2 = std::make_shared<Variable>();

    variables.insert(variables.end(), {scalar1, scalar2, distAO1, distBO1, distCO2, distDO2});

	//Polinom p1{distAO1 - A.distance(O1)};
	//Polinom p2{distBO1 - B.distance(O1)};
	//Polinom p3{distCO2 - C.distance(O2)};
	//Polinom p4{distDO2 - D.distance(O2)};
	Polinom p1{distAO1 - distance2(A, O1, variables, construction_polinoms)};
	Polinom p2{distBO1 - distance2(B, O1, variables, construction_polinoms)};
	Polinom p3{distCO2 - distance2(C, O2, variables, construction_polinoms)};
	Polinom p4{distDO2 - distance2(D, O2, variables, construction_polinoms)};
	Polinom p5{scalar1 - scalar_product_of_vectors(A, O1, B, O1)};
	Polinom p6{scalar2 - scalar_product_of_vectors(C, O2, D, O2)};

	construction_polinoms.insert(construction_polinoms.end(), {p1, p2, p3, p4, p5, p6});

	Polinom rez = (Polinom)scalar1 * scalar1 * distCO2 * distDO2 - (Polinom)scalar2 * scalar2 * distAO1 * distBO1;

	polinoms.push_back(rez);
}

void make_line_orthogonal_on_plane(Line& l, Point& A, Point&  P, Point& Q, Point& R, std::vector<Polinom>& construction_polinoms)
{
	l.set_point_coordinates(A.x(), A.y(), A.z());

	std::vector<Polinom> rez = vector_product_of_vectors(P, Q, P, R);	

    Polinom vecx = rez[0] - l.x_vec();
    Polinom vecy = rez[1] - l.y_vec();
    Polinom vecz = rez[2] - l.z_vec();

    construction_polinoms.insert(construction_polinoms.end(), {vecx, vecy, vecz});
}

void equal_numbers(Number& n1, Number& n2, std::vector<Polinom>& polinoms)
{
	Polinom p = n1.variable() - (Polinom)n2.variable();

	polinoms.push_back(p);
}

void multiply(Number& n1, int broj, Number& n2, std::vector<Polinom>& polinoms)
{
	Polinom p = n1.variable() - (Polinom)broj*n2.variable();

	polinoms.push_back(p);
}

void multiply(Number& n1, Number& n2, Number& n3, std::vector<Polinom>& polinoms)
{
	Polinom p = n1.variable() - (Polinom)n3.variable()*n2.variable();

	polinoms.push_back(p);
}

void division(Number& n1, Number& n2, Number& n3, std::vector<Polinom>& polinoms)
{
	Polinom p = (Polinom)n1.variable()*n2.variable() - n3.variable();

	polinoms.push_back(p);
}

void addition(Number& n1, Number& n2, Number& n3, std::vector<Polinom>& polinoms)
{
	Polinom p = (Polinom)n1.variable() - n2.variable() - n3.variable();

	polinoms.push_back(p);
}

void make_foot_on_plane(Point& D, Point& A, Plane& pi, std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& polinoms)
{
	point_in_plane(D, pi, polinoms);

	auto k = std::make_shared<Variable>();
	variables.insert(variables.end(), k);

	Polinom p1 = (Polinom)A.x() - D.x() - (Polinom)k*pi.x_vec();
	Polinom p2 = (Polinom)A.y() - D.y() - (Polinom)k*pi.y_vec();
	Polinom p3 = (Polinom)A.z() - D.z() - (Polinom)k*pi.z_vec();

	auto square = std::make_shared<Variable>();
	variables.insert(variables.end(), square);

	Polinom p5 = (Polinom)square - (Polinom)A.x()*A.x() + (Polinom)A.y()*A.y() + (Polinom)A.z()*A.z();

	Polinom p4 = (Polinom)k*square - (Polinom)pi.x_vec()*A.x() - (Polinom)pi.y_vec()*A.y() 
				 - (Polinom)pi.z_vec()*A.z() - pi.d_vec();

	polinoms.insert(polinoms.end(), {p1, p2, p3, p4, square});
}

void make_point_projection(Point& AProject, Point& A, Plane& pi, std::vector<std::shared_ptr<Variable>>& variables, std::vector<Polinom>& polinoms)
{
	// Polinom pi_magnitude = scalar_product_of_vectors(pi.x_vec(), pi.y_vec(), pi.z_vec(),
	// 													pi.x_vec(), pi.y_vec(), pi.z_vec());

	// Polinom scalarAPn = scalar_product_of_vectors(A.x() - (Polinom)P.x(), A.y() - (Polinom)P.y(), A.z() - (Polinom)P.z(),
	// 													pi.x_vec(), pi.y_vec(), pi.z_vec());

    // Polinom p1 = AProject.x()*pi_magnitude - A.x()*pi_magnitude + scalarAPn*pi.x_vec();
	// Polinom p2 = AProject.y()*pi_magnitude - A.y()*pi_magnitude + scalarAPn*pi.y_vec();
	// Polinom p3 = AProject.z()*pi_magnitude - A.z()*pi_magnitude + scalarAPn*pi.z_vec();

	// polinoms.insert(polinoms.end(), {p1, p2, p3});	

	Polinom rez = scalar_product_of_vectors(pi.x_vec(), pi.y_vec(), pi.z_vec(),
													pi.x_vec(), pi.y_vec(), pi.z_vec());

	auto pi_magnitude = std::make_shared<Variable>();
    auto k = std::make_shared<Variable>();

	variables.insert(variables.end(), {pi_magnitude, k});													

	Polinom p1 = pi_magnitude - rez;													
	Polinom p2 = AProject.x() - (Polinom)A.x() + (Polinom)k*pi.x_vec();
	Polinom p3 = AProject.y() - (Polinom)A.y() + (Polinom)k*pi.y_vec();
	Polinom p4 = AProject.z() - (Polinom)A.z() + (Polinom)k*pi.z_vec();
	Polinom p5 = (Polinom)k * pi_magnitude - (Polinom)A.x()*pi.x_vec() - 
					(Polinom)A.y()*pi.y_vec() - (Polinom)A.z()*pi.z_vec() - pi.d_vec();

	polinoms.insert(polinoms.end(), {p1, p2, p3, p4, p5});											
}

void collinear(Point& A, Point& B, Point& C, std::vector<Polinom>& polinoms)
{
	//we can use this -- the vectors AB and AC must be parallel
	parallel_lines(A, B, A, C, polinoms);
}

void make_plane_orthogonal_on_plane(Plane& pi, Point& A, Line& l, std::vector<Polinom>& construction_polinoms)
{
	point_in_plane(A, pi, construction_polinoms);
	orthogonal_line_plane(l, pi, construction_polinoms);
}

//pi_x = ax - ox
//pi_y = ay - oy
//pi_z = az - oz
//d = -(ax^2 + ay^2 + az^2) + (ax*ox + ay*oy + az*oz)
void make_plane_orthogonal_on_plane(Plane& pi, Sphere& s, Point& A, std::vector<Polinom>& polinoms)
{
	Polinom p1 = (Polinom)A.x() - s.centerX() - pi.x_vec();
	Polinom p2 = (Polinom)A.y() - s.centerY() - pi.y_vec();
	Polinom p3 = (Polinom)A.z() - s.centerZ() - pi.z_vec();	

	polinoms.insert(polinoms.end(), {p1, p2, p3});

	Polinom p4 = (Polinom)pi.d_vec() + (Polinom)A.x()*A.x() + (Polinom)A.y()*A.y() + (Polinom)A.z()*A.z() -
				 (Polinom)A.x()*s.centerX() - (Polinom)A.y()*s.centerY() - (Polinom)A.z()*s.centerZ();

	polinoms.push_back(p4);
}