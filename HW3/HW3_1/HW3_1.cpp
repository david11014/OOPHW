#include "HW3_1.h"

Pyramid::Pyramid()
{
	vertices = new Point[4];
}
Pyramid::Pyramid(const Pyramid & P)
{
	vertices = new Point[4];
	Pyramid(P.vertices);
}
Pyramid::Pyramid(Point* ps)
{
	vertices = new Point[4];
	SetVertices(ps);
}
Pyramid::~Pyramid()
{
	delete[] vertices;
}
Pyramid Pyramid::operator=(const Pyramid & P)
{
	Pyramid p(P);
	return p;
}

void Pyramid::SetVertices(Point * ps)
{	
	vertices[0] = ps[0];
	vertices[1] = ps[1];
	vertices[2] = ps[2];
	vertices[3] = ps[3];
}

Point Pyramid::Center()
{	
	return (vertices[1] + vertices[2] + vertices[3] - vertices[0] * 3) / 4;
}

double Pyramid::Area()
{
	double A = 0;
	A += (vertices[1] - vertices[0]).CrossArea((vertices[2] - vertices[0])) / 2;
	A += (vertices[2] - vertices[0]).CrossArea((vertices[3] - vertices[0])) / 2;
	A += (vertices[3] - vertices[0]).CrossArea((vertices[1] - vertices[0])) / 2;
	A += (vertices[1] - vertices[3]).CrossArea((vertices[2] - vertices[3])) / 2;
	return A;
}

double Pyramid::Perimeter()
{
	double A = 0;
	A += (vertices[1] - vertices[0]).Abs();
	A += (vertices[2] - vertices[0]).Abs();
	A += (vertices[3] - vertices[0]).Abs();
	A += (vertices[1] - vertices[2]).Abs();
	A += (vertices[2] - vertices[3]).Abs();
	A += (vertices[3] - vertices[1]).Abs();
	return	A;
}

double Pyramid::Volume()
{
	return  (vertices[1] - vertices[0]).Dot((vertices[2] - vertices[0]).Cross(vertices[3] - vertices[0])) / 6;
}

Cuboid::Cuboid()
{
	vertices = new Point[8];
}
Cuboid::Cuboid(Point* ps)
{
	vertices = new Point[8];
	SetVertices(ps);
}
Cuboid::Cuboid(const Cuboid & c)
{
	Cuboid(c.vertices);
}
Cuboid::~Cuboid()
{
	delete[] vertices;
}
Cuboid Cuboid::operator=(const Cuboid & c)
{
	Cuboid C(c);
	return C;
}

void Cuboid::SetVertices(Point * ps)
{
	for (int i = 0; i < 8; i++)
		vertices[i] = ps[i];
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

Cylinder Cylinder::operator=(const Cylinder & c)
{
	Cylinder C(c);
	return C;
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
