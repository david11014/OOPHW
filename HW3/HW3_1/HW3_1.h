#include <iostream>
using namespace std;

#ifndef HW3_1_H
#define HW3_1_H

class Point {

public:
	union
	{
		double p[3];
		struct
		{
			double x;
			double y;
			double z;
		};
	};

	Point() 
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	};
	Point(double a, double b, double c) :x(a), y(b), z(c) {};


	double& operator[](unsigned int i)
	{
		return p[i];
	}
	Point operator+(Point P)
	{
		Point p;
		p.x = this->x + P.x;
		p.y = this->y + P.y;
		p.z = this->z + P.z;
		return p;
	}
	Point operator-(Point P)
	{
		Point p;
		p.x = this->x - P.x;
		p.y = this->y - P.y;
		p.z = this->z - P.z;
		return p;
	}
	Point operator*(double a)
	{
		Point p;
		p.x = (this->x) * a;
		p.y = (this->y) * a;
		p.z = (this->z) * a;
		return p;
	}
	Point operator/(double a)
	{
		Point p;
		p.x = (this->x) / a;
		p.y = (this->y) / a;
		p.z = (this->z) / a;
		return p;
	}

	double Distant(Point P)
	{
		return sqrt((p[0] - P[0])*(p[0] - P[0]) + (p[1] - P[1])*(p[1] - P[1]) + (p[2] - P[2])*(p[2] - P[2]));
	}

	void show() {
		std::cout << p[0] << " " << p[1] << " " << p[2];
	}

	Point Cross(Point P)
	{
		return Point(p[1] * P[2] - p[2] * P[1], p[2] * P[0] - p[0] * P[2], p[0] * P[1] - p[1] * P[0]);
	}

	double Dot(Point P)
	{
		return p[0] * P[0] + p[1] * P[1] + p[2] * P[2];
	}

	Point Unit()
	{
		double L = Distant(Point(0, 0, 0));
		return Point(p[0] / L, p[1] / L, p[2] / L);

	}

	friend ostream& operator<<(ostream&, const Point&);

};

ostream& operator<<(ostream& os, const Point& p)
{
	os << p.x << "\t" << p.y << "\t" << p.z;
	return os;
}

class IGeometry
{
public:
	virtual double Area() = 0;
	virtual double Perimeter() = 0;
};

class Pyramid : public IGeometry
{
	
public:
	Pyramid();
	~Pyramid();
	Pyramid(const Pyramid &);
	Pyramid & operator= (const Pyramid & );
	
	void SetVertices(Point*);
	Point Center();

	virtual double Area();
	virtual double Perimeter();
private:
	Point* vertices;
};

class Cuboid : public IGeometry
{
public:
	Cuboid();
	~Cuboid();
	Cuboid(const Cuboid &);
	Cuboid & operator= (const Cuboid & );

	void SetVertices(Point*);
	double* SideLength();
	double* SideArea();	

	virtual double Area();
	virtual double Perimeter();
private:
	Point* vertices;
};

class Cylinder : public IGeometry
{
public:
	Cylinder();
	~Cylinder();
	Cylinder(const Cylinder &);
	Cylinder & operator= (const Cylinder & );

	void SetCylinder(Point,Point,double);
	double Height();
	double BottomArea();
	double SideArea();

	virtual double Area();
	virtual double Perimeter();
private:
	Point Top, Bottom;
	double r;
};


#endif // !HW3_1_H
