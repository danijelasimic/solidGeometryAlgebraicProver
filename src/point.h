#ifndef POINT_H
#define POINT_H

#include <memory>
#include <vector>

#include "geom3d_object.h"
#include "variable.h"
#include "polinom.h"

class Point : public Geom3dObject
{
public:
    Point() = default;

    Point(const std::string& name, std::vector<std::shared_ptr<Variable>>& variables, bool dependent = true)
    :Geom3dObject(name, dependent)
   	{
      coordinate_x = std::make_shared<Variable>();
      coordinate_y = std::make_shared<Variable>();
      coordinate_z = std::make_shared<Variable>();

      variables.insert(variables.end(), {coordinate_x, coordinate_y, coordinate_z});
    }

    /* This is constructor for point created in implementation in order to 
     * make implementation easier, but it is not a real point in construction.
     */
    Point(const std::shared_ptr<Variable>& x, const std::shared_ptr<Variable>& y, const std::shared_ptr<Variable>& z)
    :Geom3dObject("noname", true)
    {
      coordinate_x = x;
      coordinate_y = y;
      coordinate_z = z;
    }

   	~Point() = default;

   	void print() const override;
    
    std::shared_ptr<Variable> coordinate(int i) const;
    std::shared_ptr<Variable> x() const;
    std::shared_ptr<Variable> y() const;
    std::shared_ptr<Variable> z() const;

    /* coordinate can be 1, 2, 3 */
    void set_point_coordinate_value(int coordinate, int value = 0);
    void set_all_coordinates_values(int val1, int val2, int val3);

    void make_midpoint(Point& A, Point& B, std::vector<Polinom>& construction_polinoms);
    void translate(Point& A, int x_distance, int y_distance, int z_distance, std::vector<Polinom>& construction_polinoms);

    Polinom distance(Point& B);

    void set_shared_ptr_coordinates(std::shared_ptr<Variable>& x, std::shared_ptr<Variable>& y, std::shared_ptr<Variable>& z)
    {
      coordinate_x = x;
      coordinate_y = y;
      coordinate_z = z;
    }

private:
	std::shared_ptr<Variable> coordinate_x;
	std::shared_ptr<Variable> coordinate_y;
	std::shared_ptr<Variable> coordinate_z;

};

Polinom scalar_product_of_vectors(const Point& A, const Point& B, const Point& C, const Point& D);
std::vector<Polinom> vector_product_of_vectors(const Point& A, const Point& B, const Point& C, const Point& D);


#endif //POINT_H