#include "HW3_1.h"

Pyramid::Pyramid()
{
	vertices = new Point[4];
}
Pyramid::~Pyramid()
{
	delete[] vertices;
}

Pyramid::Pyramid(const Pyramid &)
{


}

Pyramid & Pyramid::operator=(const Pyramid & P)
{
	Pyramid *P = new Pyramid();

	

	return &P;
}

void Pyramid::SetVertices(Point *)
{


}

Point Pyramid::Center()
{
	return Point();
}

double Pyramid::Area()
{
	return 0.0;
}

double Pyramid::Perimeter()
{
	return 0.0;
}

double Pyramid::Volume()
{
	return 0.0;
}

Cuboid::Cuboid()
{
	vertices = new Point[8];
}

Cuboid::~Cuboid()
{
	delete[] vertices;
}

Cuboid::Cuboid(const Cuboid &)
{
}

Cuboid & Cuboid::operator=(const Cuboid &)
{
	Cuboid C;
	return *C;
}

void Cuboid::SetVertices(Point *)
{
}

double * Cuboid::SideLength()
{
	return nullptr;
}

double * Cuboid::SideArea()
{
	return nullptr;
}

double Cuboid::Area()
{
	return 0.0;
}

double Cuboid::Perimeter()
{
	return 0.0;
}

double Cuboid::Volume()
{
	return 0.0;
}

Cylinder::Cylinder()
{
}

Cylinder::~Cylinder()
{
}

Cylinder::Cylinder(const Cylinder &)
{
}

Cylinder & Cylinder::operator=(const Cylinder &)
{
	Cylinder C;
	return *C;
}

void Cylinder::SetCylinder(Point, Point, double)
{
}

double Cylinder::Height()
{
	return 0.0;
}

double Cylinder::BottomArea()
{
	return 0.0;
}

double Cylinder::SideArea()
{
	return 0.0;
}

double Cylinder::Area()
{
	return 0.0;
}

double Cylinder::Perimeter()
{
	return 0.0;
}

double Cylinder::Volume()
{
	return 0.0;
}
