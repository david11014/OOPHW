#pragma once

#ifndef __time_measure_H__
#define __time_measure_H__

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
namespace mytime
{
	using namespace std::chrono;
	class mytime
	{
	public:

		mytime() :dt(high_resolution_clock::duration::zero())
		{
		}
		mytime(high_resolution_clock::duration dt_) :dt(dt_)
		{
		}
		void start()
		{
			t1 = high_resolution_clock::now();
		}
		void end()
		{
			dt += high_resolution_clock::now() - t1;
		}
		void clear()
		{
			dt = high_resolution_clock::duration::zero();
		}
		void print(const std::string& s = std::string())
		{
			std::cout << s << duration_cast<milliseconds>(dt).count() << "ms.\n";
		}
		void print(const char* const s)
		{
			std::cout << s << duration_cast<milliseconds>(dt).count() << "ms.\n";
		}
		void printfps()
		{
			std::cout << fps() << "\n";
		}
		double fps()
		{
			if (dt != high_resolution_clock::duration::zero())
			{
				return 1000.0 / duration_cast<milliseconds>(dt).count();
			}
			return 0;
		}
		std::string toString()
		{
			return std::to_string(duration_cast<milliseconds>(dt).count());
		}
	private:
		high_resolution_clock::duration dt;
		high_resolution_clock::time_point t1;
	};
}
#endif