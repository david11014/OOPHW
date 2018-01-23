/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 5
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>

template<class T>  // functor class defines operator()()
class TooBig
{
private:
	T cutoff;
public:
	TooBig(const T & t) : cutoff(t) {}
	bool operator()(const T & v) { return v > cutoff; }
};

int main()
{
	TooBig<int> f100(100); // limit = 100
	int vals[10] = { 50, 100, 90, 180, 60, 210, 415, 88, 188, 201 };
	std::list<int> yadayada(vals, vals + 10); // range constructor
	std::list<int> etcetera(vals, vals + 10);

	auto outint_lb = [](int n) { std::cout << n << " "; };
	auto f100_lb = [&f100](int n) {return f100(n);};
	

	// C++0x can use the following instead
	//  list<int> yadayada = {50, 100, 90, 180, 60, 210, 415, 88, 188, 201};
	//  list<int> etcetera {50, 100, 90, 180, 60, 210, 415, 88, 188, 201};

	std::cout << "Original lists:\n";
	for_each(yadayada.begin(), yadayada.end(), outint_lb);
	std::cout << std::endl;
	for_each(etcetera.begin(), etcetera.end(), outint_lb);
	std::cout << std::endl;
	yadayada.remove_if(f100_lb);               // use a named function object
	etcetera.remove_if([](int n) {return TooBig<int>(200)(n); });   // construct a function object
	std::cout << "Trimmed lists:\n";
	for_each(yadayada.begin(), yadayada.end(), outint_lb);
	std::cout << std::endl;
	for_each(etcetera.begin(), etcetera.end(), outint_lb);
	std::cout << std::endl;
	std::cin.get();
	return 0;
}
