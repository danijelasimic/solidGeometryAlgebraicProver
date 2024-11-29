#ifndef PLANE_H
#define PLANE_H

#include <memory>
#include <vector>

#include "geom3d_object.h"
#include "variable.h"
#include "point.h"
#include "polinom.h"

class Plane : public Geom3dObject
{
public:
	Plane() = default;

	Plane(std::string& name, std::vector<std::shared_ptr<Variable>>& variables, bool dependent = true)
	:Geom3dObject(name, dependent)
	{
		_vector_x = std::make_shared<Variable>();
		_vector_y = std::make_shared<Variable>();
		_vector_z = std::make_shared<Variable>();
        _vector_d = std::make_shared<Variable>();

		variables.insert(variables.end(), {_vector_x, _vector_y, _vector_z, _vector_d});
	}

	~Plane() = default;

	void print() const override;
	void make_plane_trough_points(Point& A, Point& B, Point& C, std::vector<Polinom>& construction_polinoms);

	std::shared_ptr<Variable> x_vec() const {return _vector_x; }
	std::shared_ptr<Variable> y_vec() const {return _vector_y; }
	std::shared_ptr<Variable> z_vec() const {return _vector_z; }
    std::shared_ptr<Variable> d_vec() const {return _vector_d; }

private:
	/* Plane is determined with its vector and free parameter. */
	/* Vector variables. */
	std::shared_ptr<Variable> _vector_x;
	std::shared_ptr<Variable> _vector_y;
	std::shared_ptr<Variable> _vector_z;
    std::shared_ptr<Variable> _vector_d;	
};


#endif