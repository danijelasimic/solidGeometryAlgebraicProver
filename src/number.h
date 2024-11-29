#ifndef NUMBER_H
#define NUMBER_H

#include <memory>
#include <vector>

#include "geom3d_object.h"
#include "variable.h"
#include "point.h"
#include "polinom.h"

class Number : public Geom3dObject
{
public:
	Number() = default;

	Number(std::string& name, std::vector<std::shared_ptr<Variable>>& variables, bool dependent)
	:Geom3dObject(name, dependent)
	{
		_variable =  std::make_shared<Variable>();
		variables.push_back(_variable);
	}

	Number(std::string& name, std::vector<std::shared_ptr<Variable>>& variables, int value, bool dependent)
	:Geom3dObject(name, dependent)
	{
		_variable =  std::make_shared<Variable>();
		_variable->set_value(value);
		variables.push_back(_variable);
	}

	~Number() = default;

	void print() const override;

	std::shared_ptr<Variable> variable() const {return _variable; }

	void make_distance(Point& A, Point& B, std::vector<Polinom>& construction_polinoms);
	void make_square_distance(Point& A, Point& B, std::vector<Polinom>& construction_polinoms);
	void make_distance_line(Point& M, Point& A, Point& B, std::vector<Polinom>& construction_polinoms, 
							std::vector<std::shared_ptr<Variable>>& variables);

private:
	std::shared_ptr<Variable> _variable;	
};

#endif //NUMBER_H
