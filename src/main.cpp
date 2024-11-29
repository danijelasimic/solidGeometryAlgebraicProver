#include "main.h"

void print_wolfram(int problem_num, std::vector<Polinom>& construction_polinoms, std::vector<Polinom>& statement_polinoms, 
										  std::vector<std::shared_ptr<Variable>>& variables, bool ndg)
{

	unsigned int i;

	//check variables that are zero
	//NOT IMPLEMENTED
	// for(Polinom& poly : construction_polinoms) {
	// 	poly.check_and_set_zero_vaiable();
	// }

	// std::string dat_output_name = "../output";
	// std::vector<std::string> results;
	// boost::split(results, dat_name, [](char c){return c == '.';});
	// for (int j = 0; j < (int)results.size() - 1; j++)
	// 	dat_output_name += results[j];

	// std::cout << dat_output_name << "\n";
	std::string dat_name = "../output/" + std::to_string(problem_num) + "_input.wls";
	std::cout << "Created input file: " << dat_name << std::endl;

	std::ofstream out {dat_name};

	/* Print construction polynomials. */
	bool first = true;
	out << "poly = { ";
	for (i = 0; i<construction_polinoms.size(); i++) {
		if (!construction_polinoms[i].is_zero()) {
			if (!first) {
				out << ",\n";
			} else
				first = false;

			construction_polinoms[i].print_wolfram(out);
		}
	}
	out << "}\n\n";

	/* Print variables. */
	first = true;
	out << "variables = {";
	for(i = 0; i < variables.size(); i++) {
		if (!variables[i]->has_value()) {
			if (first) first = false;
			else
				out << ", ";
			variables[i]->print_idx(out);
		}
	}
	out << "}\n\n";

	out << "gb = GroebnerBasis[poly, variables]\n";
	out << "Print[\"Groebner Basis:\"]\n";
	/* Grebnerova baza ne sme da bude 1. */
	out << "Print[gb]\n\n\n";

	out << "Print[\"Polynomial reduce:\"]\n";
	/* Print statement polynomials, call Reduce. */
	i = 0;
	for(auto statement : statement_polinoms) {
		i++;
		out << "statement" << i << " = ";
		if (statement.is_zero()) out << "0";
		else statement.print_wolfram(out);
		out << std::endl;
		out << "Print[PolynomialReduce[statement" << i << ", gb, variables][[2]]]\n\n";
	}
	out.close();

	//ndg
	if (ndg) {
		std::string dat_name_ndg = "../output/" + std::to_string(problem_num) + "_ndg_input.wls";
		std::cout << "Created input file: " << dat_name_ndg << std::endl;
		std::ofstream out_ndg {dat_name_ndg};

		/* Print construction polynomials. */
		bool first = true;
		out_ndg << "poly = { ";
		for (i = 0; i<construction_polinoms.size(); i++) {
			if (!construction_polinoms[i].is_zero()) {
				if (!first) {
					out_ndg << ",\n";
				} else
					first = false;

				construction_polinoms[i].print_wolfram(out_ndg);
			}
		}
		out_ndg << "}\n\n";

		/* Print variables. */
		first = true;
		out_ndg << "variables = {";
		for(auto var : variables) {
			if (!var->has_value()) {
				if (first) first = false;
				else
					out_ndg << ", ";
				var->print_idx(out_ndg);
			}
		}
		if (!first) out_ndg << ", ";
		out_ndg << "zz}\n\n";

		/* Print statement polynomials, call Groebner Basis. */
		i = 0;
		for(auto statement : statement_polinoms) {
				i++;
				out_ndg << "statement" << i << " = ";
				if (statement.is_zero()) out_ndg << "0";
				else statement.print_wolfram(out_ndg);
				out_ndg << std::endl;
				out_ndg << "Print[GroebnerBasis[Insert[poly, statement" << i << "* zz - 1,1], variables]]\n\n";
		}

		out_ndg.close();	
	}
}

void print_latex(int problem_num, std::vector<Polinom>& construction_polinoms, std::vector<Polinom>& statement_polinoms, 
				 bool ndg, std::vector<NDGCondition>& ndg_polinoms)
{
	// std::string dat_output_name;
	// std::vector<std::string> results;

	// boost::split(results, dat_name, [](char c){return c == '.';});

	std::ofstream out {"../output/" + std::to_string(problem_num) + "_polynomials.tex"};

	out << "\\documentclass{article}\n";
	out << "\\begin{document}\n\n";

	out << "\\section{Construction polynomials}\n\n";

	for (auto p : construction_polinoms) {
		if (!p.is_zero()) {
			p.print_latex(out);
			out << " = 0\\\\\n";
		}
	}

	out << "\n\n\n\\section{Statement polynomials}\n\n";

	for (auto p : statement_polinoms) {
			p.print_latex(out);
			out << " = 0\\\\\n";
	}

	if (ndg)
	{
		out << "\n\n\n\\section{Non degeneracy conditons}\n\n";

		out << "\n\n\n\\begin{itemize}\n\n";	

		for(const auto& ndgP : ndg_polinoms)
		{
			ndgP.print(out);			
		}	

		out << "\n\n\n\\end{itemize}\n\n";
	}	

	out << "\\end{document}";

	out.close();
}

void simplify_polynomials(std::unordered_map<std::string, Point*>& table_points, 
							std::unordered_map<std::string, Line*>& table_lines, bool solid_exists)
{
	// If solid is constructed than special coordinates are alredy chosen,
	// so nothing else should be done
	if (solid_exists)
		return;

	int point_zeros = 0;

	//todo set zeros plane

	//find three independent points and set their coordinates to
	//(0, 0, 0)
	//(0, 0, x6)
	//(0, x8, x9)
	//todo: can x6 be 1?
	for (auto& p : table_points) {
		// if the point is not dependent, its free, and its coordinates are free, than its
		// coordinates can be simplified by leting them be zero
		if ((p.second)->dependent() == false) {
			// x coordinate set 0
			(p.second)->set_point_coordinate_value(1);
			// y coordinate set 0 (secon and first point)
			if (point_zeros <= 1)
				(p.second)->set_point_coordinate_value(2);
			// z coordinate set 0 (only first point)
			if (point_zeros == 0)
				(p.second)->set_point_coordinate_value(3);

			point_zeros++;

			if (point_zeros == 3)
				return;
		}
	}

	/*
	Provera postavljanja nula.
	for (auto& p : table_points) {
		std::cout << (p.second)->name() << "\n";
		std::cout << "x" << (p.second)->x()->variable_index() << " = " << (p.second)->x()->has_value() << (p.second)->x()->value() << std::endl;
		std::cout << "x" << (p.second)->y()->variable_index() << " = " << (p.second)->y()->has_value() << (p.second)->y()->value() << std::endl;
		std::cout << "x" << (p.second)->z()->variable_index() << " = " << (p.second)->z()->has_value() << (p.second)->z()->value() << std::endl;
		std::cout << "\n\n\n";
	}
	*/

	if (point_zeros == 3)
				return;

	// if non point is set to zero
	// then one independent line is selected and 
	// its vector is set to (0, 0, x3)
	// and its point is set to (0, 0, 0)
	// todo: can x3 be 1?
	if (point_zeros == 0) {
		for(auto& l : table_lines) {
			if ((l.second)->dependent() == false)
			{
				(l.second)->set_zeros_for_coordinates();
				return;
			}
		}
	}

}
