#ifndef EXPONENTIATION_H
#define EXPONENTIATION_H

#include "variable.h"
#include <memory>
#include <fstream>

class Exponentiation
{
public:
	Exponentiation()
	:_variable {nullptr}, _power{1}
	{}

	Exponentiation(std::shared_ptr<Variable> v, int power = 1)
	: _variable{v}, _power{power}
	{}

	~Exponentiation()
	{}

	std::shared_ptr<Variable> variable() const { return _variable; }
	int power() const { return _power; }

	void print_wolfram(std::ofstream& out) const;
	void print_latex(std::ofstream& out) const;
	bool is_zero() const;

private:
	std::shared_ptr<Variable> _variable;
	int _power = {1};
};

#endif //EXPONENTIATION