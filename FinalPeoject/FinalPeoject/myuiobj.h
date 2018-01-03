#pragma once

#ifndef __myuiobj_H__
#define __myuiobj_H__

#include <algorithm>
#include <functional>
#include <map>
#include <vector>

struct point2d
{
	int X, Y;
	bool operator==(const point2d& p) const
	{
		return X == p.X && Y == p.Y;
	}
	bool operator!=(const point2d& p) const
	{
		return !(*this == p);
	}
};

struct bbox2d
{
	point2d m, M;
	bool iscross(const point2d& p) const
	{
		return m.X <= p.X && m.Y <= p.Y && p.X <= M.X && p.Y <= M.Y;
	}
};

namespace myapp
{
	/*
	0: "Default_thickness"
	1: "Default_allowance"
	2: "Default_smooth_level"
	3: "Default_High_Limit"
	4: "Default_Low_Limit"
	*/
	class myDefaultValue
	{
	public:
		
		std::string keyname;
		std::map<std::string, double> key_container;
		std::vector<std::string> IndexToKey = { "Default_thickness","Default_allowance","Default_smooth_level", "Default_High_Limit","Default_Low_Limit" };
		std::vector<double> OriginalSetting = { 3.0,2.0,10,6.0,0.0 };

		void Add(const std::string& str, double v)
		{
			key_container[str] = v;
		}
		void Add(int i, double v)
		{
			Add(IndexToKey[i], v);
		}

		double Get(const std::string& str) const
		{
			auto it = key_container.find(str);
			if (it != key_container.end())
			{
				return it->second;
			}
			return std::nan("");
		}

		double Get(int i)
		{
			return Get(IndexToKey.at(i));
		}

		int KeyNum()
		{
			return IndexToKey.size();
		}

		void SetToOriginal()
		{
			for (int i = 0; i < KeyNum(); i++)
				Add(IndexToKey[i], OriginalSetting[i]);
		}

		myDefaultValue(const std::string& name) : keyname(name)
		{			
			SetToOriginal();
		}
	};

	class myButton :public bbox2d
	{
	public:
		myButton(std::string t, std::function<void()> dth) :text(t), dosomething(dth), isdown(false)
		{

		}

		std::string text;
		std::function<void()> dosomething;
		bool isdown;
	};

	class myButtons :public bbox2d
	{
	public:
		int btn_w = 200;
		int btn_h = 20;
		int btn_spacing_w = 10;
		//int btn_spacing_h;
		std::vector<myButton> btnarr;

		bool checkcross(const point2d& p)
		{
			bool flag = false;
			for (auto&& btn : btnarr)
			{
				btn.isdown = false;
			}
			if (iscross(p))
			{
				for (auto&& btn : btnarr)
				{
					if (btn.iscross(p))
					{
						btn.isdown = true;
						flag = true;
					}
				}
			}

			return flag;
		}
		void resize(int w, int h)
		{
			int wM = 0;
			for (auto&& btn : btnarr)
			{
				wM = std::max(wM, (int)btn.text.size() * 9);
			}
			btn_w = int(wM * 1.25);

			int cx = w / 2;
			int cy = 30;

			int allw = (btn_w + btn_spacing_w)*((int)btnarr.size()) - btn_spacing_w;
			int allh = btn_h;//(btn_h + btn_spacing_h)*btnarr.size() - btn_spacing_h;

			m = {cx - allw / 2, cy - allh / 2};
			M = {cx + allw / 2, cy + allh / 2};

			auto tempm = m;
			for (auto&& btn : btnarr)
			{
				btn.m = tempm;
				btn.M = {tempm.X + btn_w, tempm.Y + btn_h};

				tempm.X += btn_w + btn_spacing_w;
				//tempm.Y += btn_h + btn_spacing_h;
			}
		}
	};

	class myTrackBar
	{
	public:
		int value;
		int minimum;
		int maximum;
		int largechange;
		int smallchange;

		std::string text;
		std::function<void(int)> _valuechange;
		void valuechange()
		{
			_valuechange(value);
		}
	};
}
#endif