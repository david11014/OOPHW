/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <fstream>
#include "K_means.h"

#define TR_SIZE 100
#define TES_SIZE 20

using namespace std;

int main()
{

	Point2D trainP[TR_SIZE], testP[TES_SIZE];


	/*Read train file and test file*/
	ifstream trainFile, testFile;
	ofstream outFile;
	outFile.open("test-result-2.txt");
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
	Point2D P1(140.57414041910056,18.0524769614646,0), P2(16.20521473399179,42.37661357896893,0);
	KMeans KM(P1, P2, trainP, TR_SIZE);
	//cout << KM << endl;
	cout << "train done" << endl;

	cout << "find label" << endl;
	for (int i = 0; i < TES_SIZE; i++)
	{
		Point2D P = KM.FindLabel(testP[i]);
		cout << P << endl;
		outFile << P << endl;
	}

	outFile.close();
	trainFile.close();
	testFile.close();

	getchar();
	return 0;
}