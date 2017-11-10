/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 2
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include "HW02.h"
#define RANGE 100

//Point
Point::Point() : x(0), y(0) {};//#2
Point::Point(float a, float b) : x(a), y(a) {};//#3

void Point::Set_data(float a, float b)//#5
{
	x = a;
	y = b;
}
float Point::operator[](int i)const//#6
{
	if (i == 0)
		return x;
	else if (i == 1)
		return y;
	else
		return 0;
}
std::istream& operator>>(std::istream& os, Point& p)//#7
{
	float x, y;
	os >> x >> y;
	p.Set_data(x, y);
	return os;
}
std::ostream& operator<<(std::ostream& os, const Point& p)//#8
{
	os << p[0] << "\t" << p[1];
	return os;
}

//QuadtreeNode
QuadtreeNode::QuadtreeNode(const Point& sp, const Point&p, const float s) : size(s), separate_point(sp)//#14
{

	nextNode[0] = nullptr;
	nextNode[1] = nullptr;
	nextNode[2] = nullptr;
	nextNode[3] = nullptr;

	data = new Point(p);

}
QuadtreeNode::QuadtreeNode(const QuadtreeNode& QT) : size(QT.size), separate_point(QT.separate_point), data(QT.data)//#15
{
	//copy child node value
	for (int i = 0; i < 4; i++)
	{
		if (QT.nextNode[i] != nullptr)
			nextNode[i] = new QuadtreeNode(*(QT.nextNode[i]));
	}
	
}
QuadtreeNode::~QuadtreeNode()//#16
{
	//clear data
	if (data != nullptr)
	{
#ifdef DEBUG
		cout << "delete: " << *data << endl;
#endif // DEBUGs
		delete data;
	}

	//clear child node
	if (nextNode[0] != nullptr) nextNode[0]->~QuadtreeNode();
	if (nextNode[1] != nullptr) nextNode[1]->~QuadtreeNode();
	if (nextNode[2] != nullptr) nextNode[2]->~QuadtreeNode();
	if (nextNode[3] != nullptr) nextNode[3]->~QuadtreeNode();


}
bool QuadtreeNode::InsertPoint(const Point& p)//#17
{
	int i = -1;
	bool occupy = false;
	Point *nsp = new Point();

	//find point's place;
	if (p[0] >= separate_point[0] && p[1] >= separate_point[1]) // 0
	{
		i = 0;
		nsp = new Point(separate_point[0] + size / 2, separate_point[1] + size / 2);

		if (data != nullptr)
			if ((*data)[0] >= separate_point[0] && (*data)[1] >= separate_point[1])
				occupy = true;
	}
	else if (p[0] < separate_point[0] && p[1] > separate_point[1]) //1
	{
		i = 1;
		nsp = new Point(separate_point[0] - size / 2, separate_point[1] + size / 2);

		if (data != nullptr)
			if ((*data)[0] < separate_point[0] && (*data)[1] < separate_point[1])
				occupy = true;
	}
	else if (p[0] <= separate_point[0] && p[1] <= separate_point[1]) //2
	{
		i = 2;
		nsp = new Point(separate_point[0] - size / 2, separate_point[1] - size / 2);

		if (data != nullptr)
			if ((*data)[0] <= separate_point[0] && (*data)[1] <= separate_point[1])
				occupy = true;
	}
	else if (p[0] > separate_point[0] && p[1] < separate_point[1]) //3
	{
		i = 3;
		nsp = new Point(separate_point[0] + size / 2, separate_point[1] - size / 2);

		if (data != nullptr)
			if ((*data)[0] > separate_point[0] && (*data)[1] < separate_point[1])
				occupy = true;
	}


	
	if (nextNode[i] == nullptr) //add node
	{
		nextNode[i] = new QuadtreeNode(*nsp, p, size / 2);

		if (occupy)
		{
			nextNode[i]->InsertPoint(*data);
			data = nullptr;
		}

		return true;
	}
	else //find in next level
		return nextNode[i]->InsertPoint(p);

	return false;
}
Point QuadtreeNode::FindClosestPoint(const Point & p) const //#18
{
	int i = -1;
	
	//find point's place;
	if (p[0] >= separate_point[0] && p[1] >= separate_point[1]) // 0
		i = 0;
	else if (p[0] < separate_point[0] && p[1] > separate_point[1]) //1
		i = 1;
	else if (p[0] <= separate_point[0] && p[1] <= separate_point[1]) //2
		i = 2;
	else if (p[0] > separate_point[0] && p[1] < separate_point[1]) //3
		i = 3;

	if (nextNode[i] == nullptr && data != nullptr)
	{
		return *data;
	}
	else
		return nextNode[i]->FindClosestPoint(p);

}

