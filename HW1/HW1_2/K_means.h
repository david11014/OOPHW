/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <math.h>
#include <iostream>
#include <memory>  
using namespace std;

#ifndef K_MEANS_H
#define K_MEANS_H

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

	friend ostream& operator<<(ostream&, const Point2D&);

};

ostream& operator<<(ostream& os, const Point2D& p)
{
	os << p.x << ", " << p.y << " " << p.l;
	return os;
}


class KMeans
{
public:
	
	Point2D M[2];
	Point2D *trainP;
	int trSize;
	const int Count = 50;

	KMeans(Point2D P1,Point2D P2, Point2D* TP, int n)
	{
		M[0] = P1;
		M[1] = P2;
		trainP = TP;
		trSize = n;

		for (int i = 0; i < Count; i++)
		{
			
			Point2D NewP1(0,0,0), NewP2(0,0,0);
			int cd1 = 0, cd2 = 0;//兩群的記數

			for (int j = 0; j < trSize; j++)
			{
				double d1, d2;
				//紀錄點與兩群中心距離
				d1 = M[0].Distant(TP[j]);
				d2 = M[1].Distant(TP[j]);

				if (d1 < d2)
				{
					NewP1 = NewP1 + TP[j];
					cd1++;
					TP[j].l = -1;					
				}					
				else
				{
					NewP2 = NewP1 + TP[j];
					cd2++;
					TP[j].l = 1;
				}
					
			}
			
			M[0] = NewP1 / cd1;
			M[1] = NewP2 / cd2;		

		}

	}
	

	Point2D FindLabel(Point2D P)
	{
		double d1, d2;
		d1 = M[0].Distant(P);
		d2 = M[1].Distant(P);

		if (d1 < d2)
		{
			P.l = -1;
			return P;
		}
		else
		{
			P.l = 1;
			return P;
		}
	}

	friend ostream& operator<<(ostream&, const  KMeans&);
};

ostream& operator<<(ostream&os , const  KMeans& KM)
{
	os << "Group1 mean point: " << KM.M[0] << endl;
	os << "Group2 mean point: " << KM.M[1] << endl;
	os << "Train Data out" << endl;
	for (int i = 0; i < KM.trSize; i++)
	{
		os << KM.trainP[i] << endl;
	}


	return os;
}

#endif // !K_DTREE
