#include "stdafx.h"

#include <cmath>
#include "myFloatingPoint.h"
#include "myMath.h"

// unsigned int myMath::Point3D::operator < (const Point3D& P) const
// {
// 	return (x < P.x) | ((y < P.y) << 1) | ((z < P.z) << 2);
// }

myMath::Vector3D myMath::Point3D::operator-(const Point3D& p) const
{
	return {x - p.x, y - p.y, z - p.z};
}

myMath::Point3D myMath::Point3D::operator+(const Point3D& p) const
{
	return {x + p.x, y + p.y, z + p.z};
}

myMath::Point3D& myMath::Point3D::operator+=(const Point3D& p)
{
	x += p.x;
	y += p.y;
	z += p.z;
	return *this;
}

myMath::Point3D myMath::Point3D::operator+(const Vector3D& v) const
{
	return {x + v.x, y + v.y, z + v.z};
}

myMath::Point3D& myMath::Point3D::operator+=(const Vector3D& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

myMath::Point3D myMath::Point3D::operator-(const Vector3D& v) const
{
	return {x - v.x, y - v.y, z - v.z};
}

myMath::Point3D& myMath::Point3D::operator-=(const Vector3D& p)
{
	x -= p.x;
	y -= p.y;
	z -= p.z;
	return *this;
}

myMath::Point3D myMath::Point3D::operator*(const double d) const
{
	return {x*d, y*d, z*d};
}

myMath::Point3D& myMath::Point3D::operator*=(const double d)
{
	x *= d;
	y *= d;
	z *= d;
	return *this;
}

myMath::Point3D myMath::Point3D::operator/(const double d) const
{
	return {x * (1 / d), y * (1 / d), z * (1 / d)};
}

myMath::Point3D& myMath::Point3D::operator/=(const double d)
{
	x *= 1 / d;
	y *= 1 / d;
	z *= 1 / d;
	return *this;
}

myMath::Point3D::operator myMath::Vector3D() const
{
	return {x, y, z};
} 

bool myMath::Point3D::operator==(const Point3D& p) const
{
	return equal_to(p, 0);
}

bool myMath::Point3D::equal_to(const Point3D& p, double th) const
{
	if (th == 0)
	{
		return x == p.x && y == p.y && z == p.z;
	}
	return abs(x - p.x) <= th && abs(y - p.y) <= th && abs(z - p.z) <= th;
}

bool myMath::Point3D::equal_tof(const Point3D& p) const
{
	return (float)x == (float)p.x && (float)y == (float)p.y && (float)z == (float)p.z;
}

bool myMath::Point3D::near_to(const Point3D& p) const
{
	myFloat md[3] = {(float)x, (float)y, (float)z};
	return md[0] == (float)p[0] && md[1] == (float)p[1] && md[2] == (float)p[2];
}

double myMath::Point3D::distance(const Point3D& p) const
{
	return sqrt((x - p.x)*(x - p.x) + (y - p.y)*(y - p.y) + (z - p.z)*(z - p.z));
}

double myMath::Point3D::distancesq(const Point3D& p) const
{
	return (x - p.x)*(x - p.x) + (y - p.y)*(y - p.y) + (z - p.z)*(z - p.z);
}
//////////////////////////////////////////////////////////////////////////
myMath::Vector3D myMath::Vector3D::operator+(const myMath::Vector3D& V) const
{
	return {x + V.x, y + V.y, z + V.z};
}

myMath::Vector3D& myMath::Vector3D::operator+=(const myMath::Vector3D& V)
{
	x += V.x;
	y += V.y;
	z += V.z;
	return *this;
}

myMath::Vector3D myMath::Vector3D::operator-() const
{
	return {-x, -y, -z};
}

myMath::Vector3D myMath::Vector3D::operator-(const myMath::Vector3D& V) const
{
	return {x - V.x, y - V.y, z - V.z};
}

myMath::Vector3D& myMath::Vector3D::operator-=(const myMath::Vector3D& V)
{
	x -= V.x;
	y -= V.y;
	z -= V.z;
	return *this;
}

myMath::Vector3D myMath::Vector3D::operator*(const myMath::Vector3D& V) const
{
	return cross(*this, V);
}

myMath::Vector3D& myMath::Vector3D::operator*=(const myMath::Vector3D& V)
{
	return *this = *this * V;
}

myMath::Vector3D myMath::Vector3D::operator*(const double d) const
{
	return {x*d, y*d, z*d};
}

myMath::Vector3D& myMath::Vector3D::operator*=(const double d)
{
	x *= d;
	y *= d;
	z *= d;
	return *this;
}

myMath::Vector3D myMath::Vector3D::operator/(const double d) const
{
	return {x*(1 / d), y*(1 / d), z*(1 / d)};
}

myMath::Vector3D& myMath::Vector3D::operator/=(const double d)
{
	x *= (1 / d);
	y *= (1 / d);
	z *= (1 / d);
	return *this;
}

void myMath::Vector3D::normalize()
{
	double l = length();
	if (l > 0.0)
	{
		x *= 1 / l;
		y *= 1 / l;
		z *= 1 / l;
	}	
}

double myMath::Vector3D::length() const
{
	return sqrt(x*x+y*y+z*z);
}

double myMath::Vector3D::lengthsq() const
{
	return x*x+y*y+z*z;
}

double myMath::Vector3D::dot(const Vector3D& va, const Vector3D& vb)
{
	return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
}

myMath::Vector3D myMath::Vector3D::cross(const Vector3D& va, const Vector3D& vb)
{
	return {va[1] * vb[2] - vb[1] * va[2], va[2] * vb[0] - vb[2] * va[0], va[0] * vb[1] - vb[0] * va[1]};
}

bool myMath::Vector3D::is_parallel(const Vector3D& va, const Vector3D& vb)
{
	//uVector3D v1(va);
	//uVector3D v2(vb);
	//return Vector3D::cross(v1, v2).length() < 1e-15;

	myDouble md[6] = {va[1] * vb[2], vb[1] * va[2], va[2] * vb[0], vb[2] * va[0], va[0] * vb[1], vb[0] * va[1]};
	return md[0] == md[1] && md[2] == md[3] && md[4] == md[5];
}

myMath::Vector3D myMath::Vector3D::genVVec(const Vector3D& V)
{
	if (abs(V.x) < abs(V.y) && abs(V.x) < abs(V.z))
	{
		//x min
		return Vector3D::cross(V, {1,0,0});
	}
	else if(abs(V.y) < abs(V.x) && abs(V.y) < abs(V.z))
	{
		//y min
		return Vector3D::cross(V, {0,1,0});
	}
	else if (abs(V.z) < abs(V.x) && abs(V.z) < abs(V.y))
	{
		//z min
		return Vector3D::cross(V, {0,0,1});
	}
	else
	{
		//x==y, z max
		return Vector3D::cross(V, {1,0,0});
	}
}
//////////////////////////////////////////////////////////////////////////
std::pair<double, double> myMath::Line3D::base_height(const Point3D& P) const
{
	Vector3D hy = P - p;
	double base = Vector3D::dot(hy, n);
	double height = Vector3D::cross(hy, n).length();

	return std::make_pair(base, height);
}

double myMath::Line3D::distance(const Point3D& P) const
{
	return Vector3D::cross(P - p, n).length();
}

bool myMath::Line3D::Intersection(const Point3D& P) const
{
	return P.near_to(p + n*Vector3D::dot(P - p, n));
}

std::pair<bool, myMath::Point3D> myMath::Line3D::Intersection(const Line3D& L) const
{
	if (Vector3D::is_parallel(L.n, n))
	{
		std::cout << "wtf line intersection parallel\n";
		return std::make_pair(false, Point3D());
	}

	Vector3D pp = p - L.p;
	uVector3D nn = Vector3D::cross(n, L.n);
	Plane pl1(Vector3D::cross(nn, L.n), L.p);
	Plane pl2(Vector3D::cross(nn, n), p);
	auto ppp = pl1.Intersection(*this);

	if (!ppp.near_to(pl2.Intersection(L)))
	{
		std::cout << "wtf line intersection askew " << Vector3D::dot(pp, Vector3D::cross(L.n, n)) << "\n";
		return std::make_pair(false, Point3D());
	}

	return std::make_pair(true, ppp);
}
//////////////////////////////////////////////////////////////////////////
double myMath::baseTriangle::distancesq(const Point3D& P) const
{
	if (isinside(P))
	{
		auto d = Plane(*this).distance(P);
		return d*d;
	}
	else
	{
		return std::min({Segment3D(at(0), at(1)).distancesq(P), Segment3D(at(1), at(2)).distancesq(P), Segment3D(at(2), at(0)).distancesq(P)});
	}
}

double myMath::baseTriangle::distance(const Point3D& P) const
{
	return sqrt(distancesq(P));
}

bool myMath::baseTriangle::isinside(const Point3D& P) const
{
	Vector3D pa = P - at(0);
	Vector3D pb = P - at(1);
	Vector3D pc = P - at(2);

	double pab = Vector3D::dot(Vector3D::cross(pa, pb), n);
	double pbc = Vector3D::dot(Vector3D::cross(pb, pc), n);
	double pca = Vector3D::dot(Vector3D::cross(pc, pa), n);

	if (pab >= 0 && pbc >= 0 && pca >= 0)
	{
		return true;
	}
	//std::cout << pab << ", " << pbc << ", " << pca << std::endl;
	return false;
}

bool myMath::baseTriangle::onplane(const Point3D& P) const
{
	return Plane(*this).onplane(P);
}

bool myMath::baseTriangle::isslim(double limit) const
{
	uVector3D u1 = at(1) - at(0);
	uVector3D u2 = at(2) - at(0);
	uVector3D u3 = at(2) - at(1);

	double dd1 = Vector3D::cross(u1, u2).lengthsq();
	double dd2 = Vector3D::cross(u1, u3).lengthsq();
	double dd3 = Vector3D::cross(u3, u2).lengthsq();

	return std::sqrt(std::max({dd1, dd2, dd3})) < limit;
}

myMath::Point3D myMath::baseTriangle::getcen() const
{
	return (at(0) + at(1) + at(2)) / 3;
}
//////////////////////////////////////////////////////////////////////////
myMath::Plane::Plane(const baseTriangle& t) :n(t.n), cen({ (t[0].x + t[1].x + t[2].x) / 3, (t[0].y + t[1].y + t[2].y) / 3, (t[0].z + t[1].z + t[2].z) / 3 })
{
	d = Vector3D::dot(n, cen);
}

double myMath::Plane::distance(const Point3D& P) const
{
	return Vector3D::dot(n, P) - d;
}

myMath::Point3D myMath::Plane::Intersection(const Line3D& L) const
{
	double t = (d - Vector3D::dot(n, Vector3D(L.p))) / (Vector3D::dot(n, L.n));
	
	//return {L.n.x*t + L.p.x, L.n.y*t + L.p.y, L.n.z*t + L.p.z};
	return {fma(L.n.x, t, L.p.x), fma(L.n.y, t, L.p.y), fma(L.n.z, t, L.p.z)};
}

myMath::Point3D myMath::Plane::Projection(const Point3D& P) const
{
	return P + n*(-distance(P));
}

myMath::Point3D myMath::Plane::Mirror(const Point3D& P) const
{
	return P + n*(-distance(P) * 2);
}

bool myMath::Plane::onplane(const Point3D& P) const
{
	return P./*near_to*/equal_tof(Projection(P));
}

bool myMath::Plane::onplane(const baseTriangle& T) const
{
	return onplane(T.at(0)) && onplane(T.at(1)) && onplane(T.at(2));
}

bool myMath::Plane::is_parallel(const Plane& pa, const Plane& pb)
{
	return Vector3D::is_parallel(pa.n, pb.n);
}

bool myMath::Plane::is_coincide(const Plane& pa, const Plane& pb)
{
	return is_parallel(pa, pb) && pb.cen.near_to(pb.cen + pa.n * Vector3D::dot(pa.cen - pb.cen, pa.n));
}
//////////////////////////////////////////////////////////////////////////
bool myMath::Segment3D::isnull() const
{
	return pstart.near_to(pend);
}

std::pair<int, myMath::Point3D> myMath::Segment3D::Intersection(const Segment3D& L) const
{
	auto r = ((Line3D*)this)->Intersection((Line3D)L);
	if (r.first)
	{
		auto& rp = r.second;
		auto d1 = rp.near_to(pstart) || rp.near_to(pend) ? 0 : Vector3D::dot(pstart - rp, rp - pend);
		auto d2 = rp.near_to(L.pstart) || rp.near_to(L.pend) ? 0 : Vector3D::dot(L.pstart - rp, rp - L.pend);

		if (d1 >= 0 && d2 >= 0)
		{
			if (d1 == 0 && d2 == 0)
			{
				//std::cout << "3\n";
				return std::make_pair(3, r.second);
			}
			else if (d1 == 0)
			{
				//std::cout << "1\n";
				return std::make_pair(1, r.second);
			}
			else if (d2 == 0)
			{
				//std::cout << "2\n";
				return std::make_pair(2, r.second);
			}
			//std::cout << "0\n";
			return std::make_pair(0, r.second);
		}
	}

	return std::make_pair(-1, Point3D());
}

bool myMath::Segment3D::Intersection(const Point3D& P) const
{
	return ((Line3D*)this)->Intersection(P) && (Vector3D::dot(pstart - P, P - pend) >= 0 || P.near_to(pstart) || P.near_to(pend));
}

double myMath::Segment3D::distancesq(const Point3D& P) const
{
	double rr = 1;
	Point3D pp = pstart + n * Vector3D::dot(n, P - pstart);
	auto d1 = Vector3D::dot(pstart - pp, pp - pend);
	if (d1 > 0)
	{
		rr = P.distancesq(pp);
	}
	else if (d1 <= 0)
	{
		rr = std::min(P.distancesq(pstart), P.distancesq(pend));
	}

	return rr;
}

double myMath::Segment3D::distance(const Point3D& P) const
{
	return sqrt(distancesq(P));
}

double myMath::Segment3D::length() const
{
	return (pstart - pend).length();
}
//////////////////////////////////////////////////////////////////////////
std::pair<bool, myMath::Segment3D> myMath::Triangle::Intersection(const Plane& P) const
{
	if (P.onplane(*this))
	{
		std::cout << __FUNCTION__ << __LINE__ << " coincide\n";
		return std::make_pair(false, Segment3D());
	}
	for (int i = 0; i < 3; i++)
	{
		if (P.onplane(p[i]))
		{
			std::cout << __FUNCTION__ << __LINE__ << " p onplane\n";
			return std::make_pair(false, Segment3D());
		}
	}

	int flag = 0;

	flag += P.distance(a) > 0 ? 1 : 0;
	flag += P.distance(b) > 0 ? 2 : 0;
	flag += P.distance(c) > 0 ? 4 : 0;

	Line3D L1(a, b);//!3,4
	Line3D L2(a, c);//!2,5
	Line3D L3(b, c);//!1,6

	Point3D P1;
	Point3D P2;
	switch (flag)
	{
	case 1:
	case 6:
		P1 = P.Intersection(L1);
		P2 = P.Intersection(L2);
		break;
	case 2:
	case 5:
		P1 = P.Intersection(L3);
		P2 = P.Intersection(L1);
		break;
	case 3:
	case 4:
		P1 = P.Intersection(L2);
		P2 = P.Intersection(L3);
		break;
	case 0:
	case 7:
	default:
		return std::make_pair(false, Segment3D());
	}

	return std::make_pair(true, Segment3D(P1, P2));
}

double myMath::Triangle::area() const
{
	return 0.5*Vector3D::cross(b - a, c - a).length();
}
//////////////////////////////////////////////////////////////////////////
myMath::Cylinder::Cylinder(const Point3D& cen_, const uVector3D& n_, double r, double h) :cen(cen_), n(n_), radius(r), height(h)
{}

std::vector<std::shared_ptr<myMath::baseTriangle>> myMath::Cylinder::Triangulation(bool twoside) const
{
	std::vector<std::shared_ptr<myMath::baseTriangle>> resultarr;
	
	double delta = 2 * M_PI / nedge;

	Vector3D v = {0, 0, 1};
	uVector3D rr = Vector3D::cross(v, n);
	
	double ctheta = Vector3D::dot(v, n);
	double stheta = sqrt(1-ctheta*ctheta);
	double vtheta = 1 - ctheta;
	Vector3D rmat[3] = 
	{
		{rr[0] * rr[0] * vtheta + ctheta, rr[0] * rr[1] * vtheta - rr[2] * stheta, rr[0] * rr[2] * vtheta + rr[1] * stheta},
		{rr[0] * rr[1] * vtheta + rr[2]*stheta, rr[1] * rr[1] * vtheta + ctheta, rr[1] * rr[2] * vtheta - rr[0] * stheta},
		{rr[0] * rr[2] * vtheta - rr[1]*stheta, rr[1] * rr[2] * vtheta + rr[0] * stheta, rr[2] * rr[2] * vtheta + ctheta}
	};

	auto rxp = [&](const Point3D& p)->Point3D
	{
		double x = Vector3D::dot(rmat[0], p) + cen.x;
		double y = Vector3D::dot(rmat[1], p) + cen.y;
		double z = Vector3D::dot(rmat[2], p) + cen.z;

		return {x, y, z};
	};

	std::vector<Point3D> basepoint(nedge);
	std::vector<Point3D> uppoint(nedge);
	Point3D basecen = twoside ? cen - n*height : cen;
	Point3D upcen = cen + n*height;
	double z = twoside ? -height : 0;
	for (unsigned int i = 0; i < nedge; i++)
	{
		double deg = delta*i;
		double x = radius*cos(deg);
		double y = radius*sin(deg);
		basepoint[i] = rxp({x, y, z});
		uppoint[i] = rxp({x, y, height});
	}

	for (unsigned int i1 = 0; i1 < nedge; i1++)
	{
		unsigned int i2 = (i1 + 1) % nedge;
		resultarr.push_back(std::make_shared<Triangle>(basecen, basepoint[i2], basepoint[i1]));
		resultarr.push_back(std::make_shared<Triangle>(upcen, uppoint[i1], uppoint[i2]));
		resultarr.push_back(std::make_shared<Triangle>(basepoint[i1], basepoint[i2], uppoint[i2]));
		resultarr.push_back(std::make_shared<Triangle>(basepoint[i1], uppoint[i2], uppoint[i1]));
	}

	return resultarr;
}
//////////////////////////////////////////////////////////////////////////
void myMath::TransformationMatrix::LoadMat(std::string str)
{
 	using namespace std;
 	Init();
 	ifstream fin(str.c_str(), ios::in | ios::binary);
 	if (!fin.good())
 	{
 		throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
 	}
 	fin.seekg(0, ios::end);
 	string contents((size_t)fin.tellg(), '0');
 	fin.seekg(0, ios::beg);
 	fin.read(&contents[0], contents.size());
 	fin.close();
 	for (auto&& c : contents)
 	{
 		if (c == ',')
 		{
 			c = ' ';
 		}
 	}
 
 	stringstream ss(contents);
 	for (unsigned int i = 0; i < 16; i++)
 	{
 		ss >> r[i];
 	}
 	t.p[0] = r[12]; r[12] = 0;
 	t.p[1] = r[13]; r[13] = 0;
 	t.p[2] = r[14]; r[14] = 0;
}

void myMath::TransformationMatrix::SaveMat(std::string str) const
{
 	using namespace std;
 	std::ofstream fout(str.c_str(), ios::out | ios::binary);
 
 	if (!fout.good())
 	{
 		throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
 	}
 	auto arr = getMat();
 	for (unsigned int i = 0; i < 16; i++)
 	{
 		fout << arr[i] << (i % 4 == 3 ? "\r\n" : ", ");
 	}
 	fout.close();
}
//////////////////////////////////////////////////////////////////////////
double myMath::Determinant(const Point3D& A, const Point3D& B, const Point3D& C, const myMath::Point3D& D)
{
	//return ((B.x - D.x)*(C.y - D.y) - (B.y - D.y)*(C.x - D.x))*(A.z - D.z) - ((B.x - D.x)*(C.z - D.z) - (B.z - D.z)*(C.x - D.x))*(A.y - D.y) + ((B.y - D.y)*(C.z - D.z) - (B.z - D.z)*(C.y - D.y))*(A.x - D.x);

// 	double bdx = B.x - D.x;
// 	double bdy = B.y - D.y;
// 	double bdz = B.z - D.z;
// 	double cdx = C.x - D.x;
// 	double cdy = C.y - D.y;
// 	double cdz = C.z - D.z;
// 	return (bdx*cdy - bdy*cdx)*(A.z - D.z) - (bdx*cdz - bdz*cdx)*(A.y - D.y) + (bdy*cdz - bdz*cdy)*(A.x - D.x);

	return Vector3D::dot(A - D, Vector3D::cross(B - D, C - D));
}

double myMath::Determinant(const Vector3D& A, const Vector3D& B, const Vector3D& C)
{
	return Vector3D::dot(A, Vector3D::cross(B, C));
}

