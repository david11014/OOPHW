/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 2
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "HW02.h"
#define RANGE 100
//#define DEBUG

using namespace std;

int main()
{

	Point* P;
	int PNum;

	//Read point file
	ifstream PointsFile;

	PointsFile.open("Point_HW2.txt");

	if (PointsFile.is_open())
	{
		cout << "Open points data" << endl;
		char s[13];
		PointsFile.get(s, 13);
		cout << s;
		PointsFile >> s;

		PNum = stoi(s);
		cout << PNum << endl;
		float x, y;
		P = new Point[PNum];

		for (int i = 0; i < PNum; i++)
		{
			PointsFile >> P[i];
			cout << P[i] << endl;
		}
		cout << "Done load points data" << endl;
	}
	Point sp(0, 0);

	QuadtreeNode* QTRoot = new QuadtreeNode(sp, P[0], RANGE);

	for (int i = 1; i < PNum; i++)
	{
#ifdef DEBUG
		cout << P[i] << QTRoot->InsertPoint(P[i]) << endl;
#else
		QTRoot->InsertPoint(P[i]);
#endif // DEBUG		
	}

	cout << "Done Build Quadtree" << endl;

	cout << "============================================" << endl;
	cout << "�i�����ϥιw�]�������禡" << endl;
	cout << "~Point(): �����������Ŷ�" << endl;
	cout << "Point(const Point&): ���O�����L���Хi�����ƻs" << endl;
	cout << "============================================" << endl;

	Point testP;
	QuadtreeNode *QTCopy = new QuadtreeNode(*QTRoot);
		
	while (1)
	{
		cout << "Please enter the point:" << endl;
		cin >> testP;
		cout << "Quadtree closest point: " << QTRoot->FindClosestPoint(testP) << endl;
		cout << "Copy Quadtree closest point: " << QTCopy->FindClosestPoint(testP) << endl;
	}

	return 0;
}
