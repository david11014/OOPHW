/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 5
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "Point.h"

Point* LoadPoint(const char*, unsigned int&);

int main()
{
	std::chrono::steady_clock::time_point start, end;//�����ɶ����ܼ�

	start = std::chrono::steady_clock::now();//����Ū�ɶ}�l�ɶ�
	unsigned int nPoint;
	Point* point_array = LoadPoint("HW5-3.bin", nPoint);//Ū�I
	end = std::chrono::steady_clock::now();//����Ū�ɵ����ɶ�
	std::cout << "Time " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
	std::cout << "Number of point " << nPoint << std::endl;

	start = std::chrono::steady_clock::now();//�������������I�}�l�ɶ�

											 //�бN���������I�t��k���g��

	end = std::chrono::steady_clock::now();//�������������I�����ɶ�

	std::cout << "Time " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
	std::cout << /*���������I�᪺�I�Ӽ�*/ << std::endl;

	//�b���[�J�g�ɵ{���X

	delete[]point_array;
	system("pause");
	return 0;
}


Point* LoadPoint(const char* filename, unsigned int& nPoint)		// �q��r��Ū�J�h���I���
{
	using namespace std;
	Point* pPoint;

	try
	{
		if (!filename)
			throw new exception("�޼ƿ��~");

		ifstream fin;
		fin.open(filename, ios::in | ios::binary);
		if (!fin.good())
			throw new exception("�ɦW�θ��|���~! �L�k�}���ɮ�!");

		const unsigned int tmpLen = 80;
		char tmpBuff[tmpLen];

		//read Header
		fin.read(tmpBuff, tmpLen);
		if (!fin.good())
			throw new exception("�榡���X(header���~)");

		//how many point ?
		fin.read((char *)&nPoint, 4);//unsigned long, must be 4 bytes
		if (!fin.good())
			throw new exception("�榡���X(�T������ƿ��~)");

		//allocate array memory
		if (nPoint == 0)
		{
			throw new exception("NO Point!");
		}
		pPoint = new Point[nPoint];

		//read triangles
		fin.seekg(84, ios::beg);
		for (int i = 0; i < nPoint; i++)
		{
			fin.read(tmpBuff, 14);
			if (!fin.good())
				throw new exception("�榡���X");
			int gc = fin.gcount();
			if (gc < 14)
				throw;

			memcpy(&pPoint[i], tmpBuff, 12);
		}
		fin.clear();
		fin.close();
	}
	catch (exception *)
	{
		pPoint = 0;
		return (Point *)0;
	}

	return pPoint;
}
