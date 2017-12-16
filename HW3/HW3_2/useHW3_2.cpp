/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 3
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include "HW3_2.h"
#include "HW3_2.cpp"

int main()
{
	Queue<int> Q(10);
	

	cout << "add 0~9\n";
	for (int i = 0; i < 10; i++)
		Q.Add(i);

	Queue<int> Qcopy(Q); //使用拷貝建構子
	Queue<int> Qcopy2 = Q; //使用指定運算子
	cout << "queue Q:\t" << Q << endl;
	cout << "queue Qcopy:\t" << Qcopy << endl;
	cout << "queue Qcopy2:\t" << Qcopy2 << endl;
	cout << endl;
	
	cout << "add 22\n";
	Q.Add(22);
	cout << Q << endl;
	cout << endl;

	cout << "Del: " << *(Q.Delete()) << endl;
	cout << "queue Q:\t" << Q << endl;
	cout << "queue Qcopy:\t" << Qcopy << endl;
	cout << "queue Qcopy2:\t" << Qcopy2 << endl;
	cout << endl;

	cout << "add 23\n";
	Q.Add(23);
	cout << Q << endl;
	cout << endl;

	cout << "add 24\n";
	Q.Add(24);
	cout << Q << endl;
	cout << endl;

	cout << "Delete 12 time\n";
	for (int i = 0; i < 12; i++)
	{
		int*  value = Q.Delete();
		if(value != 0)
			cout << "Del: " << *value << endl;
	}
		

	getchar();


	return 0;
}