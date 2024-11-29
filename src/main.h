#ifndef MAIN_H
#define MAIN_H

#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <map>
#include <unordered_map>

#include "../src/polinom.h"
#include "../src/variable.h"
#include "../src/line.h"
#include "../src/point.h" 
#include "../src/ndgConditions.h"


void print_wolfram(int problem_num, std::vector<Polinom>& construction_polinoms, std::vector<Polinom>& statement_polinoms, 
										  std::vector<std::shared_ptr<Variable>>& variables, bool ndg);

void print_latex(int problem_num, std::vector<Polinom>& construction_polinoms, std::vector<Polinom>& statement_polinoms, 
				bool ndg, std::vector<NDGCondition>& ndg_polinoms);

void simplify_polynomials(std::unordered_map<std::string, Point*>& table_points, std::unordered_map<std::string, Line*>& table_lines, bool solid_exists);


#endif // MAIN_H
