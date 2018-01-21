/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 5
Write by david1104
github: https://github.com/david11014
********************************************************/
/*請使用release 組態編譯*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <unordered_set>
#include <ppl.h>
#include <concurrent_unordered_set.h>
#include "Point.h"
#define PARALLEL

using namespace std;
using namespace concurrency;

Point* LoadPoint(const char*, unsigned int&);

int main()
{
	std::chrono::steady_clock::time_point start, end;//紀錄時間的變數

	start = std::chrono::steady_clock::now();//紀錄讀檔開始時間
	unsigned int nPoint;
	Point* point_array = LoadPoint("HW5-3.bin", nPoint);//讀點
	end = std::chrono::steady_clock::now();//紀錄讀檔結束時間
	std::cout << "Time " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
	std::cout << "Number of point " << nPoint << std::endl;
		
	start = std::chrono::steady_clock::now();//紀錄移除重複點開始時間

#ifdef PARALLEL
	concurrent_unordered_set<Point> point_set;
	parallel_for_each(point_array, point_array + nPoint, [&point_set](Point &p){
		
#else
	unordered_set<Point> point_set;
	for (auto it = point_array; it < point_array + nPoint; it++){
		Point p = *it;		
#endif
	if (point_set.find(p) == point_set.end())
		point_set.insert(p);
}
#ifdef PARALLEL
);
#endif

	end = std::chrono::steady_clock::now();//紀錄移除重複點結束時間

	std::cout << "Time " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
	std::cout << point_set.size() << std::endl;

	//在此加入寫檔程式碼
	cout << "Start write file...\n";
	ofstream outFile;
	outFile.open("OutPoint.txt");
	int count = 0;
	for (auto p : point_set)
	{
		outFile << p << endl;
		count++;
	}
	cout << "End write file\n";
	delete[]point_array;
	system("pause");
	return 0;
}


Point* LoadPoint(const char* filename, unsigned int& nPoint)		// 從文字檔讀入多個點資料
{
	using namespace std;
	Point* pPoint;

	try
	{
		if (!filename)
			throw new exception("引數錯誤");

		ifstream fin;
		fin.open(filename, ios::in | ios::binary);
		if (!fin.good())
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");

		const unsigned int tmpLen = 80;
		char tmpBuff[tmpLen];

		//read Header
		fin.read(tmpBuff, tmpLen);
		if (!fin.good())
			throw new exception("格式不合(header錯誤)");

		//how many point ?
		fin.read((char *)&nPoint, 4);//unsigned long, must be 4 bytes
		if (!fin.good())
			throw new exception("格式不合(三角網格數錯誤)");

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
				throw new exception("格式不合");
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

