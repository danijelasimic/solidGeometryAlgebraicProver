#ifndef TERM_H
#define TERM_H

#include <vector>
#include <fstream>

#include "exponentiation.h"

class Term
{
public:
	Term(int coefficient = 1)
	:_coefficient {coefficient}
	{}

	Term(const std::vector<Exponentiation>& exponentiations, const int coefficient = 1)
	: _exponentiations {exponentiations}, _coefficient {coefficient}
	{}

	/* makes term coeff*v^power */
	Term(int coefficient, std::shared_ptr<Variable> v, int power=1)
	:_coefficient {coefficient}
	{
		Exponentiation e{v, power};
		_exponentiations.push_back(e);
	}

	/* makes term coeff*v1*v2*v3*... */
	Term(int coefficient, std::initializer_list<std::shared_ptr<Variable>> variables)
	:_coefficient {coefficient}
	{
		for(auto v : variables) {
			_exponentiations.push_back(Exponentiation{v, 1});
		}
	}

	~Term() = default;

	void add_exponentiation(Exponentiation& e);

	bool is_zero() const;
	void print_wolfram(std::ofstream& out) const;
	void print_latex(std::ofstream& out) const;

	void multiply_coefficient_with_scalar(int scalar) { _coefficient *= scalar; }

	// multiply two terms
	Term& operator*=(const Term& rhs);
  	friend Term operator*(const Term& lhs, const Term& rhs)
  	{
  		Term solution{lhs.exponentiations(), lhs.coefficient()};
  		solution._exponentiations.insert(solution._exponentiations.end(), rhs._exponentiations.begin(), rhs._exponentiations.end());
  		solution._coefficient *= rhs._coefficient;
    	return solution;
  	}


  	int coefficient() const { return _coefficient; }
  	std::vector<Exponentiation> exponentiations() const { return _exponentiations; }

private:
	std::vector<Exponentiation> _exponentiations;
	int _coefficient;
};

#endif //TERM_H