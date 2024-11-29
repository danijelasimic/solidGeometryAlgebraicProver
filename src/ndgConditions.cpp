#include "ndgConditions.h"

void NDGCondition::print(std::ofstream& out) const
{
    out << "\n\n\n\\item " << _ndgText << "\n\n";
    for (const auto& p : _ndgPolinoms) {
        p.print_latex(out);
        out << " $\\neq 0$\\\\\n";
    }
}

void insert_ndg_points_not_equal(Point& A, Point& B, std::vector<NDGCondition>& ndg_polinoms)
{
    // (a x − b x )^2 + (a y − b y )^2 + (a z − b z )^2
	//evry line is: ax^2 - 2*ax*bx + bx^2
	Polinom p{
		Term{1, A.x(), 2}, Term{-2, {A.x(), B.x()}}, Term{1, B.x(), 2},
		Term{1, A.y(), 2}, Term{-2, {A.y(), B.y()}}, Term{1, B.y(), 2},
		Term{1, A.z(), 2}, Term{-2, {A.z(), B.z()}}, Term{1, B.z(), 2}
	};

    ndg_polinoms.push_back(NDGCondition{"Point " + A.name() + " is not equal to point " + B.name() + ".", {p}});
}

void insert_ndg_number_not_zero(Number& m,  std::vector<NDGCondition>& ndg_polinoms)
{
    Polinom p{Term{-1, m.variable()}};

    ndg_polinoms.push_back(NDGCondition{"Number " + m.name() + " is not zero.", {p}});
}

void insert_ndg_line_vector_not_zero(Line& l, std::vector<NDGCondition>& ndg_polinoms)
{
    Polinom p = (Polinom)l.x_vec()*l.x_vec() + (Polinom)l.y_vec()*l.y_vec() + (Polinom)l.z_vec()*l.z_vec();

    ndg_polinoms.push_back(NDGCondition{"Line " + l.name() + " is well-defined. Line vector is not zero.", {p}});
}

void insert_ndg_lines_not_parallel(Line& l1, Line& l2, std::vector<NDGCondition>& ndg_polinoms)
{
    std::string text = "Lines " + l1.name() + " and " + l2.name() + " are not parallel.";
    Polinom vec_abx = Polinom{Term{1, l1.x_vec()}};
	Polinom vec_aby = Polinom{Term{1, l1.y_vec()}};
	Polinom vec_abz = Polinom{Term{1, l1.z_vec()}};

	Polinom vec_cdx = Polinom{Term{1, l2.x_vec()}};
	Polinom vec_cdy = Polinom{Term{1, l2.y_vec()}};
	Polinom vec_cdz = Polinom{Term{1, l2.z_vec()}};

	std::vector<Polinom> p = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cdx, vec_cdy, vec_cdz);

    ndg_polinoms.push_back(NDGCondition{text, p});
}

void insert_ndg_lines_not_skew(Line& l1, Line& l2, std::vector<NDGCondition>& ndg_polinoms)
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
	ndg_polinoms.push_back(NDGCondition{"Lines " + l1.name() + " and " + l2.name() + " are not skew.",
		mix_product_of_vectors(vec_acx, vec_acy, vec_acz,
								vec_abx, vec_aby, vec_abz,
								vec_cdx, vec_cdy, vec_cdz)});

}

void insert_ndg_lines_not_parallel(Point& A, Point& B, Point& C, Point& D, std::vector<NDGCondition>& ndg_polinoms)
{
    Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cdx = Polinom{Term{1, D.x()}, Term{-1, C.x()}};
	Polinom vec_cdy = Polinom{Term{1, D.y()}, Term{-1, C.y()}};
	Polinom vec_cdz = Polinom{Term{1, D.z()}, Term{-1, C.z()}};

    std::string text = "Line one, determined with points " + A.name() + " and " + B.name() + ", " +
                       "is not paralle with line two, detemined with points " + C.name() + 
                       " and " + D.name() + ".";
    std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cdx, vec_cdy, vec_cdz);

	ndg_polinoms.push_back(NDGCondition{text, rez});
}

void insert_ndg_lines_not_skew(Point& A, Point& B, Point& C, Point& D, std::vector<NDGCondition>& ndg_polinoms)
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

    std::string text = "Line one, determined with points " + A.name() + " and " + B.name() + ", " +
                       "is not skew with line two, detemined with points " + C.name() + 
                       " and " + D.name() + ".";

	// mix product of AC AB CD
	ndg_polinoms.push_back(NDGCondition{text,
		mix_product_of_vectors(vec_acx, vec_acy, vec_acz,
								vec_abx, vec_aby, vec_abz,
								vec_cdx, vec_cdy, vec_cdz)});
}

void insert_ndg_plane_vector_not_zero(Plane& pi, std::vector<NDGCondition>& ndg_polinoms)
{
    Polinom p = (Polinom)pi.x_vec()*pi.x_vec() + (Polinom)pi.y_vec()*pi.y_vec() + (Polinom)pi.z_vec()*pi.z_vec();

    ndg_polinoms.push_back(NDGCondition{"Plane " + pi.name() + " is well-defined. Plane vector is not zero.", {p}});
}

void insert_ndg_plane_and_line_not_parallel(Line& l, Plane& pi, std::vector<NDGCondition>& ndg_polinoms)
{
    	Polinom rez = Polinom{
		Term{1, {pi.x_vec(), l.x_vec()}},
		Term{1, {pi.y_vec(), l.y_vec()}},
		Term{1, {pi.z_vec(), l.z_vec()}}};
	
        ndg_polinoms.push_back(NDGCondition{
            "Line " + l.name() + " is not parallel with plane " + pi.name() + ".", 
            rez
        });
}

void insert_ndg_plane_and_line_not_parallel(Point& A, Point& B, Plane& pi, std::vector<NDGCondition>& ndg_polinoms)
{
    Polinom rez = Polinom{
		Term{1, {pi.x_vec(), B.x()}}, Term{-1, {pi.x_vec(), A.x()}},
		Term{1, {pi.y_vec(), B.y()}}, Term{-1, {pi.y_vec(), A.y()}},
		Term{1, {pi.z_vec(), B.z()}}, Term{-1, {pi.z_vec(), A.z()}}};

    ndg_polinoms.push_back(NDGCondition{
            "Line, determined with points " + A.name() + " and " + B.name() +
            " is not parallel with plane " + pi.name() + ".", 
            rez
        });
}

void insert_ndg_plane_and_line_not_parallel(Line& l, Point& A, Point& B, Point& C, std::vector<NDGCondition>& ndg_polinoms)
{
    //make plane that have points A, B, C
    Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cax = Polinom{Term{1, C.x()}, Term{-1, A.x()}};
	Polinom vec_cay = Polinom{Term{1, C.y()}, Term{-1, A.y()}};
	Polinom vec_caz = Polinom{Term{1, C.z()}, Term{-1, A.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cax, vec_cay, vec_caz);

    Polinom poly = rez[0]*l.x_vec() + rez[1]*l.y_vec() + rez[2]*l.z_vec();
	
    ndg_polinoms.push_back(NDGCondition{
        "Line " + l.name() + " is not parallel with plane that is detemined with points " +
        A.name() + ", " + B.name() + " and " + C.name() + ".", 
        poly
    });
}

void insert_ndg_plane_and_line_not_parallel(Point& P, Point& Q,
                                Point& A, Point& B, Point& C, std::vector<NDGCondition>& ndg_polinoms)
{
    //make plane that have points A, B, C
    Polinom vec_abx = Polinom{Term{1, B.x()}, Term{-1, A.x()}};
	Polinom vec_aby = Polinom{Term{1, B.y()}, Term{-1, A.y()}};
	Polinom vec_abz = Polinom{Term{1, B.z()}, Term{-1, A.z()}};

	Polinom vec_cax = Polinom{Term{1, C.x()}, Term{-1, A.x()}};
	Polinom vec_cay = Polinom{Term{1, C.y()}, Term{-1, A.y()}};
	Polinom vec_caz = Polinom{Term{1, C.z()}, Term{-1, A.z()}};

	std::vector<Polinom> rez = vector_product_of_vectors(vec_abx, vec_aby, vec_abz,
									vec_cax, vec_cay, vec_caz);

    Polinom poly = rez[0]*(P.x() - (Polinom)Q.x()) + 
                  rez[1]*(P.y() - (Polinom)Q.y()) + 
                  rez[2]*(P.z() - (Polinom)Q.z());                                    

	
    ndg_polinoms.push_back(NDGCondition{
        "Line determined with points " + P.name() + " and " + Q.name() +
        " is not parallel with plane that is detemined with points " +
        A.name() + ", " + B.name() + " and " + C.name() + ".", 
        poly
    });
}         

void insert_ndg_point_not_on_line(Point& A, Line& l, std::vector<NDGCondition>& ndg_polinoms)
{
	std::vector<Polinom> poly;

	Polinom p1{Term{1, {A.x(), l.y_vec()}}, 
	           Term{-1, {l.x_p(), l.y_vec()}},
	           Term{-1, {A.y(), l.x_vec()}},
	           Term{1, {l.y_p(), l.x_vec()}}};
	poly.push_back(p1);

	Polinom p2{Term{1, {A.x(), l.z_vec()}}, Term{-1, {l.x_p(), l.z_vec()}},
			   Term{-1, {A.z(), l.x_vec()}}, Term{1, {l.z_p(), l.x_vec()}}};
	poly.push_back(p2);

	Polinom p3{Term{1, {A.y(), l.z_vec()}}, Term{-1, {l.y_p(), l.z_vec()}},
				Term{-1, {A.z(), l.y_vec()}}, Term{1, {l.z_p(), l.y_vec()}}};
	poly.push_back(p3);

	ndg_polinoms.push_back(NDGCondition{
        "Point " + A.name() + 
        " is not  on the line " +
        l.name() + ".", 
        poly
    });
}