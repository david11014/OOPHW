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
		float p[3];
		struct
		{
			float x;
			float y;
			float z;
		};
	};

	Point()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	};
	Point(float a, float b, float c) :x(a), y(b), z(c) {};

	float& operator[](unsigned int i)
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
	Point operator*(float a)
	{
		Point p;
		p.x = (this->x) * a;
		p.y = (this->y) * a;
		p.z = (this->z) * a;
		return p;
	}
	Point operator/(float a)
	{
		Point p;
		p.x = (this->x) / a;
		p.y = (this->y) / a;
		p.z = (this->z) / a;
		return p;
	}

	friend inline bool operator==(const Point &a, const Point &b) {
		return (a.x == b.x) && (a.y == b.y) && (a.z == b.z);
	}
	friend inline bool operator!=(const Point &a, const Point &b) {
		return !(a==b);
	}

	float Distant(Point P)
	{
		return sqrt((p[0] - P[0])*(p[0] - P[0]) + (p[1] - P[1])*(p[1] - P[1]) + (p[2] - P[2])*(p[2] - P[2]));
	}

	float Abs()
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
	float CrossArea(Point P)
	{
		return (this->Cross(P)).Abs();
	}

	float Dot(Point P)
	{
		return p[0] * P[0] + p[1] * P[1] + p[2] * P[2];
	}

	Point Unit()
	{
		float L = Distant(Point(0, 0, 0));
		return Point(p[0] / L, p[1] / L, p[2] / L);

	}

	size_t hash()
	{
		size_t hx = std::hash<float>()(x);
		size_t hy = std::hash<float>()(y);
		size_t hz = std::hash<float>()(z);

		return ((hx ^ (hy << 1)) >> 1) ^ (hz << 1);
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
	struct hash<Point>
	{
		size_t operator()(Point& k) const
		{
			//double x, y, z;
			//x = k.x;
			//y = k.y;
			//z = k.z;

			size_t hx = hash<float>()(k.x);
			size_t hy = hash<float>()(k.y);
			size_t hz = hash<float>()(k.z);

			return ((hx ^ (hy << 1)) >> 1) ^ (hz << 1) ;
		}
	};
}

#endif