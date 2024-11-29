#ifndef SPHERE_H
#define SPHERE_H

#include <memory>
#include <vector>

#include "geom3d_object.h"
#include "variable.h"
#include "point.h"
#include "polinom.h"
#include "number.h"

class Sphere : public Geom3dObject
{
public:
	/* Carefful, doesn't do anything! */
	Sphere() = default;

	Sphere(std::string& name, bool dependent = true)
	:Geom3dObject(name, dependent)
	{}

	~Sphere() = default;

	void print() const override;

	std::shared_ptr<Variable> centerX() const {return _centerX; }
    std::shared_ptr<Variable> centerY() const {return _centerY; }
    std::shared_ptr<Variable> centerZ() const {return _centerZ; }

	std::shared_ptr<Variable> radius() const {return _radius; }

	void set_center_and_radius(Point& center, Number& num);

	void make_sphere_through_four_points(Point& A, Point& B, Point& C, Point& D, Point& Center, std::vector<Polinom>& construction_polinoms);

private:
	/* Sphere is determined with its center and radius. */
	/* Center variables. */
    /* TODO: not great, it would be bettter to do shared_ptr<Point>
     * but it's problem with bison, I cannot make shared_ptr in bison
     * file, thus, I cannot make shared_ptr for Point, only Point*
     */
	std::shared_ptr<Variable> _centerX;
    std::shared_ptr<Variable> _centerY;
    std::shared_ptr<Variable> _centerZ;

	// Radius
	std::shared_ptr<Variable> _radius;	
};

#endif //SPHERE_H
