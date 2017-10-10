/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1

********************************************************/
#include <iostream>
#include <fstream>
#include "K_DTree.h"

#define TR_SIZE 100
#define TES_SIZE 20

using namespace std;

int main()
{
	
	Point2D trainP[TR_SIZE], testP[TES_SIZE];
	

	/*Read train file and test file*/
	ifstream outFile,trainFile,testFile;
	outFile.open("test-result-1.txt", std::ifstream::out);
	trainFile.open("train-data.txt", std::ifstream::in);
	testFile.open("test-data.txt", std::ifstream::in);
	
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

	trainTree.Root->show();

	cout << "\n" <<trainTree.Root->up << " " << trainTree.Root->layer << endl;

	//trainTree.show();
	//Point2D P = trainTree.FindNear(testP[0]);
	
	//P.show();

	outFile.close();
	trainFile.close();
	testFile.close();

	getchar();
	return 0;
}