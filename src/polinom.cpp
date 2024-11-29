#include "polinom.h"
#include <type_traits>

void Polinom::add_term(const Term& t)
{
	_terms.push_back(t);
}


bool Polinom::is_zero() const
{
	bool zero = true;
	for (auto t : _terms) 
		if (!t.is_zero()) {
			zero = false;
			break;
		}
	return zero;
}

void Polinom::print_wolfram(std::ofstream& out) const
{
	bool first = true;
	for (auto x : _terms)
		if (!x.is_zero()) {
			if (!first) {
				out << " + ";
			} else 
				first = false;
			x.print_wolfram(out);
		}
}

void Polinom::print_latex(std::ofstream& out) const
{
	bool first = true;
	for (auto x : _terms)
		if (!x.is_zero()) {
			if (first) {
				first = false;
			} else
				out << " + ";
			x.print_latex(out);
		}
}

Polinom& Polinom::operator+=(const Polinom& rhs) 
{
	_terms.insert(_terms.end(), rhs.terms().begin(), rhs.terms().end());
	return *this; 
}

Polinom& Polinom::operator-=(const Polinom& rhs)
{
	Polinom rhs1 = rhs * (-1);
	_terms.insert(_terms.end(), rhs1.terms().begin(), rhs1.terms().end());
	return *this; 	
}

Polinom& Polinom::operator*=(const int scalar)
{
	for(auto& t : _terms) {
		t.multiply_coefficient_with_scalar(scalar);
	}
	return *this;
}


Polinom scalar_product_of_vectors(const Polinom& vx, const Polinom& vy, const Polinom& vz,
									const Polinom& wx, const Polinom& wy, const Polinom& wz)
{
	Polinom rez = vx * wx + vy * wy + vz * wz;
	return rez;
}


std::vector<Polinom> vector_product_of_vectors(const Polinom& vx, const Polinom& vy, const Polinom& vz,
									const Polinom& wx, const Polinom& wy, const Polinom& wz)
{
	std::vector<Polinom> rez;

	Polinom i = vy*wz - vz*wy;
	Polinom j = vz*wx - vx*wz;
	Polinom k = vx*wy - vy*wx;

	rez.insert(rez.end(), {i, j, k});

	return rez;
}									

Polinom mix_product_of_vectors(const Polinom& vx, const Polinom& vy, const Polinom& vz,
									const Polinom& wx, const Polinom& wy, const Polinom& wz,
									const Polinom& tx, const Polinom& ty, const Polinom& tz)
{
	Polinom rez = vx*wy*tz + wx*ty*vz + vy*wz*tx
					- vz*wy*tx - wz*ty*vx - vy*wx*tz;

	return rez;
}


		 