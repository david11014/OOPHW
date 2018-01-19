/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 5
Write by david1104
github: https://github.com/david11014
********************************************************/
#ifndef POINT_H
#define POINT_H

#include <math.h>
#include <iostream>
#include <ctime> 

using namespace std;

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

	double Abs()
	{
		return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	}

	void show()
	{
		std::cout << p[0] << " " << p[1] << " " << p[2];
	}

	Point Cross(Point P)
	{
		return Point(p[1] * P[2] - p[2] * P[1], p[2] * P[0] - p[0] * P[2], p[0] * P[1] - p[1] * P[0]);
	}
	double CrossArea(Point P)
	{
		return (this->Cross(P)).Abs();
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

	friend ostream& operator<<(ostream& os, const Point& p);

};

inline
ostream& operator<<(ostream& os, const Point& p)
{
	os << p.x << "\t" << p.y << "\t" << p.z;
	return os;
}

namespace std
{
	template <>
	struct hash< Point>
	{
		size_t operator()(const  Point& k) const
		{
			size_t hx = hash<double>()(k.x);
			size_t hy = hash<double>()(k.y);
			size_t hz = hash<double>()(k.z);

			return ((hx ^ (hy << 1)) >> 1) ^ (hz << 1) ;
		}
	};
}

#endif