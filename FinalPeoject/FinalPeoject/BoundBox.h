#pragma once

#ifndef __BoundBox_H__
#define __BoundBox_H__

#include "myMath.h"

namespace myModel
{
	using namespace myMath;
	using namespace std;

	class BoundBox
	{
	public:
		Point3D m;
		Point3D M;
		
		bool iscross(const Point3D& i) const;
		bool iscross(const BoundBox& i) const;
		bool iscross(const Line3D& i) const;
		pair<Point3D, Point3D> Intersection(const Line3D& i) const;
		double distance(const Point3D& p) const;
		double distancesq(const Point3D& p) const;
		Point3D bboxcenter() const;
		Point3D bboxcorner(const int ind) const;
		vector<Point3D> bboxcorner() const;
		Point3D bboxcorneradj(const int ind, double fac = 1.0 / 4) const;

		static BoundBox genbbox(const vector<Point3D>& ps);
		static BoundBox genbbox(const baseTriangle& bt);
	};

	class Prism
	{
	public:
		myMath::Polygon poly;
		Point3D cen;
		double height;
		//bool twoside;

		bool isinside(const Point3D& p) const
		{
			if (std::abs(p.z) > height)
			{
				return false;
			}

			return poly.isinside(p);
		}

		BoundBox getBoundBox() const
		{
			return {{poly.xm*poly.scale, poly.ym*poly.scale, -height}, {poly.xM*poly.scale, poly.yM*poly.scale, height}};
		}
	};
	
	
}

#endif