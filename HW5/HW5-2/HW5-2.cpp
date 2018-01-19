/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 5
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <initializer_list>>

using namespace std;

double average_list(const initializer_list<double> il)
{
	double sum = 0;
	int n = il.size();
	for (auto it : il)
		sum += it;

	return sum / n;
}

template <typename T>
T average_list(const initializer_list<T> il)
{
	T sum = 0;
	int n = il.size();
	for (auto it : il)
		sum += it;

	return sum / n;
}

int main()
{
	using namespace std;
	// list of double deduced from list contents
	auto q = average_list({ 15.4, 10.7, 9.0 });
	cout << q << endl;
	// list of int deduced from list contents
	cout << average_list({ 20, 30, 19, 17, 45, 38 }) << endl;
	// forced list of double
	auto ad = average_list<double>({ 'A', 70, 65.33 });
	cout << ad << endl;
	cin.get();
	return 0;
}