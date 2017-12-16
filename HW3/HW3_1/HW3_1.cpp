/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 3
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <random>
#include <ctime>
#include "HW3_1.h"
#define PI 3.14159265358979323846
#define R() rand() % 10


//Pyramid建構子 解構子
Pyramid::Pyramid()
{
	vertices = new Point[4];
	
	vertices[0] = Point(R(), R(), R());
	vertices[1] = Point(R(), R(), R());
	vertices[2] = Point(R(), R(), R());
	vertices[3] = Point(R(), R(), R());
}
Pyramid::Pyramid(Point* ps)
{
	vertices = new Point[4];
	SetVertices(ps);
}
Pyramid::Pyramid(const Pyramid &p)
{
	vertices = new Point[4];
	SetVertices(p.vertices);
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
	return  abs((vertices[1] - vertices[0]).Dot((vertices[2] - vertices[0]).Cross(vertices[3] - vertices[0]))) / 6;
}

//Cuboid建構子 解構子
Cuboid::Cuboid()
{
	vertices = new Point[8];
	SetVertices(Point(R(), R(), R()), Point(R(), R(), R()));
}
Cuboid::Cuboid(const Cuboid & c)
{
	vertices = new Point[8];
	SetVertices(c.vertices);
}
Cuboid::Cuboid(Point p1 ,Point p2)
{
	vertices = new Point[8];
	SetVertices(p1, p2);	
}
Cuboid::Cuboid(Point* ps)
{
	vertices = new Point[8];
	SetVertices(ps);
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

void Cuboid::SetVertices(Point * ps) //
{
	for (int i = 0; i < 8; i++)
		vertices[i] = ps[i];
}

void Cuboid::SetVertices(Point p1, Point p2)
{
	vertices[0] = p1;
	vertices[1] = Point(p2.x, p1.y, p1.z);
	vertices[2] = Point(p2.x, p2.y, p1.z);
	vertices[3] = Point(p1.x, p2.y, p1.z);
	vertices[4] = Point(p1.x, p2.y, p2.z);
	vertices[5] = Point(p1.x, p1.y, p2.z);
	vertices[6] = Point(p2.x, p1.y, p2.z);
	vertices[7] = p2;
}

double * Cuboid::SideLength()
{
	Point p = vertices[0] - vertices[7];
	double *ans = new double[3];

	ans[0] = abs(p[0]);
	ans[1] = abs(p[1]);
	ans[2] = abs(p[2]);

	return ans;
}

double * Cuboid::SideArea()
{
	double *l = SideLength();
	double *ans = new double[3];

	ans[0] = l[0] * l[1];
	ans[1] = l[1] * l[2];
	ans[2] = l[2] * l[0];

	return ans;
}

double Cuboid::Area()
{
	double *A = SideArea();
	return (A[0] + A[1] + A[2]) * 2;
}

double Cuboid::Perimeter()
{
	double *l = SideLength();
	return (l[0] + l[1] + l [2]) * 4;
}

double Cuboid::Volume()
{
	double *l = SideLength();
	return l[0] * l[1] * l [2];
}

//Cylinder 建構子 解構子
Cylinder::Cylinder() 
{	
	Top = Point(R(), R(), R());
	Bottom = Point(R(), R(), R());
	r = R();
}

Cylinder::~Cylinder() {};

Cylinder::Cylinder(const Cylinder & c)
{
	Top = c.Top;
	Bottom = c.Bottom;
	r = c.r;
}

Cylinder::Cylinder(Point t, Point b, double R)
{
	SetCylinder(t, b, R);
}

Cylinder Cylinder::operator=(const Cylinder & c)
{
	Cylinder C(c);
	return C;
}

void Cylinder::SetCylinder(Point t, Point b, double R)
{
	Top = t;
	Bottom = b;
	r = R;
}

double Cylinder::Height()
{
	return Top.Distant(Bottom);
}

double Cylinder::BottomArea()
{
	return PI * r * r;
}

double Cylinder::SideArea()
{
	return 2 * PI * r * Height();
}

double Cylinder::Area()
{
	return BottomArea() + SideArea();
}

double Cylinder::Perimeter()
{
	return 4 * PI * r;
}

double Cylinder::Volume()
{
	return BottomArea() * Height();
}
