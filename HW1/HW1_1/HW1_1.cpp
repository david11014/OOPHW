/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <fstream>
#include "K_DTree.h"

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
	ifstream trainFile,testFile;
	ofstream outFile;
	outFile.open("test-result-1.txt");
	trainFile.open("train-data.txt");
	testFile.open("test-data.txt");
	
	if (trainFile.is_open())
	{
		cout << "Open train data" << endl;
		for (int i = 0; i < TR_SIZE; i++)
		{
			trainFile >> trainP[i].x;
			trainFile >> trainP[i].y;
			trainFile >> trainP[i].l;
			//cout << trainP[i].x << " " << trainP[i].y << " " << trainP[i].l << endl;
		}
		cout << "Done load train data" << endl;
	}
	if (testFile.is_open())
	{
		cout << "Open test data" << endl;
		for (int i = 0; i < TES_SIZE; i++)
		{
			testFile >> testP[i].x;
			testFile >> testP[i].y;
			testFile >> testP[i].l;
			//cout << testP[i].x << " " << testP[i].y << " " << testP[i].l << endl;
		}
		cout << "Done load test data" << endl;
	}	
	
	//init K-D tree
	KDTree trainTree(trainP,TR_SIZE);	

	cout << "train done" << endl;
	
	cout << "find label" << endl;
	for (int i = 0; i < TES_SIZE; i++)
	{
		cout << "find:" << testP[i] << endl;
		Point2D P = trainTree.FindNear(testP[i]);
		Point2D TP = Find(testP[i], trainP, TR_SIZE);	

		//P.x = testP[i].x;
		//P.y = testP[i].y;
		double D = P.Distant(testP[i]);
		if (P.l != testP[i].l)
		{
			cout << testP[i].l << endl;
			P.l = 0;
		}
			

		cout << P << " " << D << endl;		
		cout << TP << " " << TP.Distant(testP[i]) << endl;

		outFile << P << " " << D << endl;
		outFile << TP << " " << TP.Distant(testP[i]) << endl;
	}
	//
	//Point2D t(35.38895558336125, -24.516727931346917, 1);

	//Point2D P = trainTree.FindNear(t);
	//Point2D TP = Find(t, trainP, TR_SIZE);
	//double D = P.Distant(t);
	//cout << P << " " << D << endl;
	//cout << TP << " " << TP.Distant(t) << endl;

	cout << "all done." << endl;
	outFile.close();
	trainFile.close();
	testFile.close();
	
	getchar();
	return 0;
}