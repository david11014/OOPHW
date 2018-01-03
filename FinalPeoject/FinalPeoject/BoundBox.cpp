#include "stdafx.h"

#include "BoundBox.h"

using namespace std;
using namespace myMath;
using namespace myModel;

bool BoundBox::iscross(const Point3D& i) const
{
	return m.x <= i.x && m.y <= i.y && m.z <= i.z && i.x <= M.x && i.y <= M.y && i.z <= M.z;
}

bool BoundBox::iscross(const BoundBox& i) const
{
	return m.x <= i.M.x && m.y <= i.M.y && m.z <= i.M.z && i.m.x <= M.x && i.m.y <= M.y && i.m.z <= M.z;
}

bool BoundBox::iscross(const Line3D& i) const
{
	auto mM = (m + M) / 2;
	auto p = i.p + i.n*Vector3D::dot(i.p - mM, i.n);
	return iscross(p);
}

pair<Point3D, Point3D> BoundBox::Intersection(const Line3D& i) const
{
	double sm = DBL_MAX;
	double sM = DBL_MAX;

	if (i.n.x != 0.0)
	{
		double nm, nM;
		nm = (m.x - i.p.x) / i.n.x;
		nM = (M.x - i.p.x) / i.n.x;
		if (i.n.x < 0)
		{
			swap(nm, nM);
		}
		if (iscross(i.p + i.n*nm)/*abs(sm) > abs(nm)*/)
		{
			sm = nm;
		}
		if (iscross(i.p + i.n*nM)/*abs(sM) > abs(nM)*/)
		{
			sM = nM;
		}
	}
	if (i.n.y != 0.0)
	{
		double nm, nM;
		nm = (m.y - i.p.y) / i.n.y;
		nM = (M.y - i.p.y) / i.n.y;
		if (i.n.y < 0)
		{
			swap(nm, nM);
		}
		if (iscross(i.p + i.n*nm)/*abs(sm) > abs(nm)*/)
		{
			sm = nm;
		}
		if (iscross(i.p + i.n*nM)/*abs(sM) > abs(nM)*/)
		{
			sM = nM;
		}
	}
	if (i.n.z != 0.0)
	{
		double nm, nM;
		nm = (m.z - i.p.z) / i.n.z;
		nM = (M.z - i.p.z) / i.n.z;
		if (i.n.z < 0)
		{
			swap(nm, nM);
		}
		if (iscross(i.p + i.n*nm)/*abs(sm) > abs(nm)*/)
		{
			sm = nm;
		}
		if (iscross(i.p + i.n*nM)/*abs(sM) > abs(nM)*/)
		{
			sM = nM;
		}
	}

	return make_pair(i.p + i.n*sm, i.p + i.n*sM);
}

double BoundBox::distance(const Point3D& p) const
{
	return sqrt(distancesq(p));
}

double BoundBox::distancesq(const Point3D& p) const
{
	Point3D pp;
	pp.x = p.x < m.x ? m.x : p.x < M.x ? p.x : M.x;
	pp.y = p.y < m.y ? m.y : p.y < M.y ? p.y : M.y;
	pp.z = p.z < m.z ? m.z : p.z < M.z ? p.z : M.z;

	return p.distancesq(pp);
}

Point3D BoundBox::bboxcenter() const
{
	return (m + M) / 2;
}

Point3D BoundBox::bboxcorner(const int ind) const
{
	return {(ind & 1 ? m[0] : M[0]), (ind & 2 ? m[1] : M[1]), (ind & 4 ? m[2] : M[2])};
}

vector<Point3D> BoundBox::bboxcorner() const
{
	vector<Point3D> result(8);
	for (unsigned int i = 0; i < 8; i++)
	{
		result[i] = bboxcorner(i);
	}
	return result;
}

Point3D BoundBox::bboxcorneradj(const int ind, double fac) const
{
	int ind1 = ind % 3;
	int ind2 = ind / 3;
	double len[3] = {(M[0] - m[0]) * fac, (M[1] - m[1]) * fac, (M[2] - m[2]) * fac};

	auto cor = bboxcorner(ind2);
	cor.p[ind1] += (~ind2 & (1 << ind1)) ? -len[ind1] : len[ind1];
	return cor;
}

BoundBox BoundBox::genbbox(const vector<Point3D>& ps)
{
	Point3D m = {DBL_MAX, DBL_MAX, DBL_MAX};
	Point3D M = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
	for (auto&& p : ps)
	{
		m = Point3D::min(m, p);
		M = Point3D::max(M, p);
	}

	return {m, M};
}

BoundBox BoundBox::genbbox(const baseTriangle& bt)
{
	auto& a = bt[0].p;
	auto& b = bt[1].p;
	auto& c = bt[2].p;
	auto xx = minmax({a[0], b[0], c[0]});
	auto yy = minmax({a[1], b[1], c[1]});
	auto zz = minmax({a[2], b[2], c[2]});

	return {{xx.first, yy.first, zz.first}, {xx.second, yy.second, zz.second}};
}