#include "exponentiation.h"

void Exponentiation::print_wolfram(std::ofstream& out) const
{
	if (_variable->has_value()) {
		out << _variable->value();

		if (_power != 1)
			out << "^" << _power;
	} else {
		out << "x" << _variable->variable_index();

		if (_power != 1)
			out << "^" << _power;
	}
}
	
void Exponentiation::print_latex(std::ofstream& out) const
{
	if (_variable->has_value()) {
		out << _variable->value();

		if (_power != 1)
			out << "^" << _power;
	} else {
		out << "x_{" << _variable->variable_index() << "}";

		if (_power != 1)
			out << "^" << _power;
	}
}

bool Exponentiation::is_zero() const
{
	return _variable->is_zero();
}
