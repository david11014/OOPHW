/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <fstream>
#include "PLA.h"

//#define DEBUG

#define TR_SIZE 100
#define TES_SIZE 20

using namespace std;

int main()
{

	Point2D trainP[TR_SIZE], testP[TES_SIZE];


	/*Read train file and test file*/
	ifstream trainFile, testFile;
	ofstream outFile;
	outFile.open("test-result-3.txt");
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
		}
		cout << "Done load test data" << endl;
	}
	
	PLA pla(trainP, TR_SIZE);
	cout << "train done" << endl;

	cout << "find label" << endl;
	for (int i = 0; i < TES_SIZE; i++)
	{
		Point2D P = pla.FindLabel(testP[i]);

#ifdef DEBUG
		cout << P << " " << P.l + testP[i].l << endl;
		outFile << P << "\t" << P.l + testP[i].l << endl;
#else
		cout << P << endl;
		outFile << P << endl;
#endif 
	}

	outFile.close();
	trainFile.close();
	testFile.close();

	getchar();
	return 0;
}