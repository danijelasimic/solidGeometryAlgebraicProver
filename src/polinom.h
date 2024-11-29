#ifndef POLINOM_H
#define POLINOM_H

#include <vector>
#include <fstream>
#include "term.h"

class Polinom 
{
public:
	Polinom() = default;

	/* Enable constructor p{Term{..}, Term{...}, Term{...}, ...} */
	Polinom(std::initializer_list<Term> l) 
	: _terms(l) 
	{}

	/* Enable cast from variable to polinom. */
	Polinom(const std::shared_ptr<Variable>& v)
		: Polinom {Term{1, v}}
	{}

	/* Enable cast from initeger (coefficient) to polinom. */
	Polinom(int coeff)
		: Polinom {Term{coeff}}
	{}

	~Polinom() = default;

	void append(std::initializer_list<Term> l) {
		_terms.insert(_terms.end(), l.begin(), l.end());
	}

	void add_term(const Term& t);

	void print_wolfram(std::ofstream& out) const;
	void print_latex(std::ofstream& out) const;
	bool is_zero() const;

	Polinom& operator+=(const Polinom& rhs);
	friend Polinom operator+(const Polinom& lhs, const Polinom& rhs)
	{
  	Polinom res;
    for(const Term& t1 : lhs.terms())
      res.add_term(t1);
    for(const Term& t1 : rhs.terms()) {
      res.add_term(t1);
    }
    return res;
	}

  Polinom& operator-=(const Polinom& rhs);
	friend Polinom operator-(const Polinom& lhs, const Polinom& rhs) 
  {
      Polinom res;
      for(const Term& t1 : lhs.terms())
        res.add_term(t1);
      for(const Term& t1 : rhs.terms()) {
        Term t2 = t1;
        t2.multiply_coefficient_with_scalar(-1);
        res.add_term(t2);
      }
      return res;
  }

  	// polinom multiply
	friend Polinom operator*(const Polinom& lhs, const Polinom& rhs) 
  {
  		Polinom res;
			for(const Term& t1 : lhs.terms())
				for(const Term& t2 : rhs.terms())
  				res.add_term(t1 * t2);

  		return res;
  }

  	//scalar multiply with polynom
	Polinom& operator*=(const int scalar);
	friend Polinom operator*(const Polinom& lhs, const int rhs) 
	{
    Polinom res;
    for(const Term& t1 : lhs.terms()) {
      Term t2 = t1;
      t2.multiply_coefficient_with_scalar(rhs);
      res.add_term(t2);
    }
    return res;
	}  	
	

  	std::vector<Term> terms() const { return _terms; }

private:
	std::vector<Term> _terms;
};

/*
  V*W = (vx, vy, vz)*(wx, wy, wz) = vx*wx + vy*wy + vz*wz
*/
Polinom scalar_product_of_vectors(const Polinom& vx, const Polinom& vy, const Polinom& vz,
									const Polinom& wx, const Polinom& wy, const Polinom& wz);

/* Ideja sa templatom je: tipovi mogu biti samo std::shared_ptr<Variable> i Polinom
 * i onda hocemo da mozemo da saljemo parametre razlicitog tipa funkciji.
 * Ovo sam jednostavnije resila tako sto sam napravila konstruktor 
 * Polinom(const std::shared_ptr<Variable>& v)
 * I on sada radi automatsko kastovanje.
 */								

/*
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
Polinom scalar_product_of_vectors(T1 vx, T2 vy, T3 vz,
									T4 wx, T5 wy, T6 wz)
{
	Polinom pvx, pvy, pvz, pwx, pwy, pwz;

	if constexpr(std::is_same_v<T1, Polinom> == true)
		pvx = vx;
	else 
		pvx = Polinom{Term{1, vx}};

	if constexpr(std::is_same_v<T2, Polinom>)
		pvy = vy;
	else 
		pvy = Polinom{Term{1, vy}};

	if constexpr(std::is_same_v<T3, Polinom>)
		pvz = vz;
	else 
		pvz = Polinom{Term{1, vz}};

	if constexpr(std::is_same_v<T4, Polinom>)
		pwx = wx;
	else 
		pwx = Polinom{Term{1, wx}};		

	if constexpr(std::is_same_v<T5, Polinom>)
		pwy = wy;
	else 
		pwy = Polinom{Term{1, wy}};		

	if constexpr(std::is_same_v<T6, Polinom>)
		pwz = wz;
	
	if constexpr(std::is_same_v<T6, std::shared_ptr<Variable>>)
		pwz = Polinom{Term{1, wz}};	

	Polinom rez = pvx * pwx + pvy * pwy + pvz * pwz;
	return rez;	
}
*/
/*
  V x W = (vx, vy, vz)x(wx, wy, wz)
*/
std::vector<Polinom> vector_product_of_vectors(const Polinom& vx, const Polinom& vy, const Polinom& vz,
									const Polinom& wx, const Polinom& wy, const Polinom& wz);

/*
  V * (W x P) = (vx, vy, vz) * ((wx, wy, wz) x (tx, ty, tz))
*/
Polinom mix_product_of_vectors(const Polinom& vx, const Polinom& vy, const Polinom& vz,
									const Polinom& wx, const Polinom& wy, const Polinom& wz,
									const Polinom& tx, const Polinom& ty, const Polinom& tz);
						

#endif //POLINOM_H
