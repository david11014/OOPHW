/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 2
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <fstream>
#include "Tree.h"
//#define DEBUG
#define TR_SIZE 100
#define TES_SIZE 20

using namespace std;

Point2D Find(Point2D P, Point2D *trainP, int n)
{
	double D = DBL_MAX;
	int I;
	for (int i = 0; i < n; i++)
	{
		if (P.Distant(trainP[i]) < D)
		{
			D = P.Distant(trainP[i]);
			I = i;
		}
			
	}
	return trainP[I];
}

int main()
{
	
	Point2D trainP[TR_SIZE], testP[TES_SIZE];
	

	/*Read train file and test file*/
	ifstream PointsFile;
	ofstream outFile;
	outFile.open("test-result-1.txt");
	PointsFile.open("Point_HW2.txt");
	
	
	if (PointsFile.is_open())
	{
		cout << "Open train data" << endl;
		for (int i = 0; i < TR_SIZE; i++)
		{
			PointsFile >> trainP[i].x;
			PointsFile >> trainP[i].y;
			PointsFile >> trainP[i].l;
			//cout << trainP[i].x << " " << trainP[i].y << " " << trainP[i].l << endl;
		}
		cout << "Done load train data" << endl;
	}
	
	
	//init K-D tree
	KDTree trainTree(trainP,TR_SIZE);	

	cout << "train done" << endl;
	
	cout << "find label" << endl;
	for (int i = 0; i < TES_SIZE; i++)
	{	
		Point2D TP = Find(testP[i], trainP, TR_SIZE);
		Point2D P = trainTree.FindNear(testP[i]);

#ifdef DEBUG

		cout << P << " "<< P.l + testP[i].l << endl;
		outFile << P << "\t" << P.l + testP[i].l << endl;
#else
		cout << P << endl;
		outFile << P << endl;
#endif 			
		
	}

	
	getchar();
	return 0;
}