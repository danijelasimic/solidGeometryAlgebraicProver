#ifndef LINE_H
#define LINE_H

#include <memory>
#include <vector>

#include "geom3d_object.h"
#include "variable.h"
#include "point.h"
#include "polinom.h"

class Line : public Geom3dObject
{
public:
	Line() = default;

	Line(const std::string& name, std::vector<std::shared_ptr<Variable>>& variables, bool dependent = true)
	:Geom3dObject(name, dependent)
	{
		_vector_x = std::make_shared<Variable>();
		_vector_y = std::make_shared<Variable>();
		_vector_z = std::make_shared<Variable>();

		variables.insert(variables.end(), {_vector_x, _vector_y, _vector_z});
	}

	Line(const std::string& name, 
		 const std::shared_ptr<Variable>& vx, const std::shared_ptr<Variable>& vy, const std::shared_ptr<Variable>& vz,
	     const std::shared_ptr<Variable>& px, const std::shared_ptr<Variable>& py, const std::shared_ptr<Variable>& pz,
		 bool dependent = true)
	:Geom3dObject(name, dependent), _vector_x {vx}, _vector_y {vy}, _vector_z {vz},
	_point_x {px}, _point_y {py}, _point_z {pz}
	{}

	~Line() = default;

	void print() const override;
	void make_line_trough_points(Point& A, Point& B, std::vector<Polinom>& construction_polinoms);

	std::shared_ptr<Variable> x_vec() const {return _vector_x; }
	std::shared_ptr<Variable> y_vec() const {return _vector_y; }
	std::shared_ptr<Variable> z_vec() const {return _vector_z; }

	std::shared_ptr<Variable> x_p() const {return _point_x; }
	std::shared_ptr<Variable> y_p() const {return _point_y; }
	std::shared_ptr<Variable> z_p() const {return _point_z; }

	// make canonical line (x = 0)
	void set_zeros_for_coordinates();
	void set_point_coordinates(const std::shared_ptr<Variable>& px, const std::shared_ptr<Variable>& py,
								const std::shared_ptr<Variable>& pz);

	Point point()
	{
		Point p;

		p.set_shared_ptr_coordinates(_point_x, _point_y, _point_z);

		return p;
	}

private:
	/* Line is determined with its vector and point on that line. */
	/* Vector variables. */
	std::shared_ptr<Variable> _vector_x;
	std::shared_ptr<Variable> _vector_y;
	std::shared_ptr<Variable> _vector_z;

	/* Point variables. */
	/* TODO; maybe its better to make pointer to Point, but for now its unnecesary complication. */
	std::shared_ptr<Variable> _point_x;
	std::shared_ptr<Variable> _point_y;
	std::shared_ptr<Variable> _point_z;
	
};

#endif //LINE_H
