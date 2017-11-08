/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <math.h>
#include <iostream>
#include <memory>  
using namespace std;

#define DEBUG
#ifndef QUAD_TREE_H
#define QUAD_TREE_H

class Point
{
private: //#1
	float x;
	float y;
public:
	Point();//#2
	Point(float, float);//#3
	~Point() {}//#4
	void Set_data(float, float);//#5
	float operator[](int i)const;//#6
	friend std::istream& operator >> (std::istream&, Point&);//#7
	friend std::ostream& operator<<(std::ostream&, const Point&);//#8
};
class QuadtreeNode
{
private:
	const Point* data; //#9
	QuadtreeNode* nextNode[4]; //#10
	const Point separate_point; //#11
	const float size;//#12

public:
	QuadtreeNode() = delete;//#13
	QuadtreeNode(const Point&, const Point&, const float);//#14
	QuadtreeNode(const QuadtreeNode&);//#15
	~QuadtreeNode();//#16
	bool InsertPoint(const Point&);//#17
	Point FindClosestPoint(¡K¡K¡K.)const;//#18
	
};

#endif // !QUAD_TREE_H
