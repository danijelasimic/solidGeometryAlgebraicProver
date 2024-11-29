#include "term.h"

void Term::add_exponentiation(Exponentiation& e)
{
	_exponentiations.push_back(e);
}

bool Term::is_zero() const
{
	if (_coefficient == 0) return true;

	for (auto e : _exponentiations) {
		if (e.is_zero()) return true;
	}

	return false;
}

void Term::print_wolfram(std::ofstream& out) const
{
	if (_exponentiations.size() != 0) {
		if (_coefficient == -1)
			out << "-";

		if (_coefficient != 1 && _coefficient != -1)
			out << _coefficient << "*";

		for (int i=0; i<(int)_exponentiations.size() - 1; i++) {
			_exponentiations[i].print_wolfram(out);
			out << "*";
		}
		_exponentiations[_exponentiations.size()-1].print_wolfram(out);
	} else
		out << _coefficient;

}

void Term::print_latex(std::ofstream& out) const
{
	out << "$";

	if (_coefficient == -1)
		out << "-";

	if (_coefficient != 1 && _coefficient != -1)
		out << _coefficient;

	for (auto e : _exponentiations) 
		e.print_latex(out);

	out << "$";
}


Term& Term::operator*=(const Term& rhs)
{
	_coefficient *= rhs.coefficient();
	_exponentiations.insert(_exponentiations.end(), rhs._exponentiations.begin(), rhs._exponentiations.end());

	return *this;
}
