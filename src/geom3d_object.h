#ifndef GEOM_3D_OBJECT_H
#define GEOM_3D_OBJECT_H

#include <iostream>
#include <string>

class Geom3dObject
{
public:
    Geom3dObject()
    {}

    Geom3dObject(const std::string& name, bool dependent = true)
    :_name {name}, _dependent {dependent}
    {
    	// std::cout << "Kreiran objekat: " << _name << std::endl;
    }

    virtual ~Geom3dObject()
    {}

    virtual void print() const = 0;

    std::string name() const { return _name; }

    bool dependent() const { return _dependent; }

private:
	std::string _name;
    bool _dependent; 
    /* If object is depnding on other objects than this value is true. For example, if line is determined by two already
       created points, than for that line this value is true.
       Most likely, this value will be true except for free created point (make_point) and for solid (make_cube).
       This value is important for simlifying polynomialas by adding coordinates.
    */
};

#endif // GEOM_3D_OBJECT_H