/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <math.h>
#include <iostream>
#include <ctime> 
using namespace std;
//#define DEBUG
#ifndef PLA_H
#define PLA_H

class Point2D {

public:
	union
	{
		double p[2];
		struct
		{
			double x;
			double y;
		};
	};
	int l;

	Point2D() {};
	Point2D(double a, double b, int c)
	{
		x = a;
		y = b;
		l = c;
	}

	double& operator[](unsigned int i)
	{
		return p[i];
	}
	Point2D operator+(Point2D P)
	{
		Point2D p;
		p.x = this->x + P.x;
		p.y = this->y + P.y;
		p.l = 0;
		return p;
	}
	Point2D operator-(Point2D P)
	{
		Point2D p;
		p.x = this->x - P.x;
		p.y = this->y - P.y;
		p.l = 0;
		return p;
	}
	Point2D operator*(double x)
	{
		Point2D p;
		p.x = (this->x) * x;
		p.y = (this->y) * x;
		p.l = 0;
		return p;
	}
	Point2D operator/(double x)
	{
		Point2D p;
		p.x = (this->x) / x;
		p.y = (this->y) / x;
		p.l = 0;
		return p;
	}

	double Distant(Point2D P)
	{
		return sqrt((p[0] - P[0])*(p[0] - P[0]) + (p[1] - P[1])*(p[1] - P[1]));
	}

	void show() {
		std::cout << p[0] << " " << p[1] << " " << l;
	}

	double Cross(Point2D P)
	{
		return this->x * P[1] - this->y * P[0];
	}

	double Dot(Point2D P)
	{
		return this->x * P[0] + this->y * P[1];
	}

	Point2D Unit()
	{
		double L = sqrt((this->x * this->x) + (this->y * this->y));
		return Point2D(this->x / L, this->y / L, 0);

	}

	friend ostream& operator<<(ostream&, const Point2D&);

};

ostream& operator<<(ostream& os, const Point2D& p)
{
	os << p.x << ", " << p.y << " " << p.l;
	return os;
}

class PLA
{
public:

	Point2D Mid;
	Point2D V;
	Point2D *trainP;
	int trSize;
	const int Count = 500;

	PLA(Point2D* TP, int n)
	{

		srand(time(NULL));
		trSize = n;
		trainP = TP;
		Mid = Point2D(0.0, 0.0, 0);

		for (int i = 0; i < trSize; i++)
		{
			Mid = Mid + trainP[i];
		}
		Mid = Mid / trSize;	

		V = Point2D(-1,0, 0); //init vector

		for (int i = 0; i < Count;)
		{
			int r = rand() % trSize;
			Point2D v = trainP[r] - Mid;

			double C = V.Cross(v);

			//需要修正方向的狀況
			if ((V.Dot(v) < 0 && trainP[r].l ==-1) || (V.Dot(v) > 0 && trainP[r].l == 1))
			{
				V = V + v * trainP[r].l;
				i++;
#ifdef  DEBUG
				cout << "V: " << V.Unit() << endl;
#endif
			}

		}
		V = V.Unit();
		cout << *this << endl;
	}

	Point2D FindLabel(Point2D P)
	{
		Point2D v = P - Mid;
		double D = V.Dot(v);
		if (D > 0)
			return Point2D(P[0], P[1], 1);
		else 
			return Point2D(P[0], P[1], -1);

	}

	friend ostream& operator<<(ostream&, const  PLA&);
};

ostream& operator<<(ostream&os, const  PLA& pla)
{
	os << "Mid:" << pla.Mid << " V: " << pla.V;
	return os;
}

#endif // !K_DTREE

