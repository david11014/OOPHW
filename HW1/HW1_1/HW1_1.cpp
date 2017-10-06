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
	
	KDTree trainTree;
	Node *R = trainTree.Root;
	int mid = (int)(TR_SIZE / 2);
	R->P = trainP[mid];
	R->layer = 0;
	R->up = nullptr;
	
	trainTree.MakeTree(R, trainP, trainP + mid, 0);
	trainTree.MakeTree(R, trainP + mid + 1, trainP + TR_SIZE + 1, 1);



	outFile.close();
	trainFile.close();
	testFile.close();

	getchar();
	return 0;
}