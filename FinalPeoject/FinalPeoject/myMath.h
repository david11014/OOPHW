#pragma once

#ifndef __myMath_H__
#define __myMath_H__

#include <algorithm>
#include <array>
#include <assert.h>
#include <memory>
#include <vector>

#define union_T_xyz_p(T, x, y, z, p) \
union\
{\
	T p[3]; \
	struct\
	{\
		T x;\
		T y;\
		T z;\
	};\
}\

namespace myMath
{
	class Vector3D
	{
	public:
		union_T_xyz_p(double, x, y, z, p);

		Vector3D() = default;
		
		const double& operator[](size_t ind) const
		{
			assert(ind<3);
			return p[ind];
		}
		double& operator[](size_t ind)
		{
			assert(ind<3);
			return p[ind];
		}
		const double& at(size_t ind) const
		{
			assert(ind < 3);
			return p[ind];;
		}
		double& at(size_t ind)
		{
			assert(ind < 3);
			return p[ind];;
		}
		Vector3D operator+(const Vector3D& V) const;
		Vector3D& operator+=(const Vector3D& V);
		Vector3D operator-() const;
		Vector3D operator-(const Vector3D& V) const;
		Vector3D& operator-=(const Vector3D& V);
		Vector3D operator*(const Vector3D& V) const;
		Vector3D& operator*=(const Vector3D& V);
		Vector3D operator*(const double d) const;
		Vector3D& operator*=(const double d);
		Vector3D operator/(const double d) const;
		Vector3D& operator/=(const double d);

		void normalize();
		double length() const;
		double lengthsq() const;
		static double dot(const Vector3D& va, const Vector3D& vb);
		static Vector3D cross(const Vector3D& va, const Vector3D& vb);
		static bool is_parallel(const Vector3D& va, const Vector3D& vb);
		static Vector3D genVVec(const Vector3D& V);

		std::array<double, 3> todarr() const
		{
			return {p[0], p[1], p[2]};
		}
		std::array<float, 3> tofarr() const
		{
			return {(float)p[0], (float)p[1], (float)p[2]};
		}
	};

	class uVector3D : public Vector3D
	{
	public:
		uVector3D() :Vector3D({ 1, 0, 0 }){}

		uVector3D(const Vector3D& V)
		{
			Vector3D v(V);
			v.normalize();
			x = v.x;
			y = v.y;
			z = v.z;
		}
		uVector3D& operator=(const Vector3D& V)
		{
			Vector3D v(V);
			v.normalize();
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}
		uVector3D& operator=(Vector3D&& v)
		{
			v.normalize();
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}
	};

	class Point3D
	{
	public:
		union_T_xyz_p(double, x, y, z, p);

		Point3D() = default;
// 		Point3D(double X, double Y, double Z) :x(X), y(Y), z(Z){}

		const double& operator[](size_t ind) const
		{
			assert(ind < 3);
			return p[ind];
		}
		double& operator[](size_t ind)
		{
			assert(ind < 3);
			return p[ind];
		}
		const double& at(size_t ind) const
		{
			assert(ind < 3);
			return p[ind];
		}
		double& at(size_t ind)
		{
			assert(ind < 3);
			return p[ind];
		}
		//inline unsigned int operator < (const Point3D& P) const
		unsigned int operator < (const Point3D& P) const
		{
			return (x < P.x) + ((y < P.y) << 1) + ((z < P.z) << 2);
		}
		Vector3D operator-(const Point3D& p) const;
		Point3D operator+(const Point3D& p) const;
		Point3D& operator+=(const Point3D& v);
		Point3D operator+(const Vector3D& v) const;
		Point3D& operator+=(const Vector3D& v);
		Point3D operator-(const Vector3D& v) const;
		Point3D& operator-=(const Vector3D& v);
		Point3D operator*(const double d) const;
		Point3D& operator*=(const double d);
		Point3D operator/(const double d) const;
		Point3D& operator/=(const double d);

		operator Vector3D() const;
		
		bool operator==(const Point3D& p) const;
		bool equal_to(const Point3D& p, double th = 1e-6) const;
		bool equal_tof(const Point3D& p) const;
		bool near_to(const Point3D& p) const;
		double distance(const Point3D& p) const;
		double distancesq(const Point3D& p) const;
		static bool erasezero(Point3D& p)
		{
			bool r = false;
			for (int i = 0; i < 3;i++)
			{
				if (p.p[i] != 0 && abs(p.p[i]) < FLT_EPSILON)
				{
					//std::cout << i << " " << p.p[i] << " " << std::hexfloat << p.p[i] << std::defaultfloat << std::endl;
					p.p[i] = 0;
					r = true;
				}				
			}
			return r;
		}
		static Point3D min(const Point3D& p1, const Point3D& p2)
		{
			return {std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)};
		}
		static Point3D max(const Point3D& p1, const Point3D& p2)
		{
			return {std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)};
		}

		std::array<double, 3> todarr() const
		{
			return {p[0], p[1], p[2]};
		}
		std::array<float, 3> tofarr() const
		{
			return {(float)p[0], (float)p[1], (float)p[2]};
		}

		char* ToString()
		{
			std::stringstream strStream;

			strStream << "x:" << p[0] << " y:" << p[1] << " z:" << p[2];			

			std::string str = strStream.str();
			std::vector<char> cstr(str.c_str(), str.c_str() + str.size() + 1);
			char* ca = &cstr[0];
			return ca;
		}
	};

	class Line3D
	{
	public:
		uVector3D n;
		Point3D p;

		Line3D(){};
		Line3D(const Vector3D& N, const Point3D& P) :n(N), p(P){}
		Line3D(const Point3D& P1, const Point3D& P2) :n(P2 - P1), p(P1){}

		std::pair<double, double> base_height(const Point3D& P) const;
		double distance(const Point3D& P) const;
		bool Intersection(const Point3D& P) const;
		std::pair<bool, Point3D> Intersection(const Line3D& L) const;
	};
	
	class baseTriangle
	{
	public:
		uVector3D n;

		baseTriangle(const Vector3D& V) :n(V){}

		virtual const Point3D& operator[](size_t ind) const = 0;
		virtual Point3D& operator[](size_t ind) = 0;

		virtual const Point3D& at(size_t ind) const = 0;
		virtual Point3D& at(size_t ind) = 0;
		virtual double area() const = 0;
		double distancesq(const Point3D& P) const;
		double distance(const Point3D& P) const;
		bool isinside(const Point3D& P) const;
		bool onplane(const Point3D& P) const;
		bool isslim(double limit = 0.000291) const;
		Point3D getcen() const;
	};

	class Plane
	{
	public:
		double d;//ax+by+cz=d
		uVector3D n;
		Point3D cen;		

		explicit Plane(const baseTriangle&);
		Plane(const Vector3D& N, const Point3D& C) :n(N), cen(C){ d = Vector3D::dot(n, cen); }
		Plane(double A, double B, double C, double D) :d(D), n({ A, B, C }), cen({ d/A/3, d/B/3, d/C/3 })
		{
		}

		double distance(const Point3D& P) const;
		Point3D Intersection(const Line3D& L) const;
		Point3D Projection(const Point3D& P) const;
		Point3D Mirror(const Point3D& P) const;
		bool onplane(const Point3D& P) const;
		bool onplane(const baseTriangle& T) const;

		static bool is_parallel(const Plane& pa, const Plane& pb);
		static bool is_coincide(const Plane& pa, const Plane& pb);
	};

	class Segment3D :public Line3D
	{
	public:
		Point3D pstart;
		Point3D pend;
		
		Segment3D(){};
		Segment3D(const Point3D& PStart, const Point3D& PEnd) :Line3D(PStart, PEnd), pstart(PStart), pend(PEnd)
		{			
		}

		bool isnull() const;
		std::pair<int, Point3D> Intersection(const Segment3D& L) const;
		bool Intersection(const Point3D& P) const;
		double distancesq(const Point3D& P) const;
		double distance(const Point3D& P) const;
		double length() const;
		
		void reversal()
		{
			std::swap(pstart, pend);
			n = -n;
		}
		void updateSegment3D()
		{
			n = pend - pstart;
			p = pstart;
		}
	};
	
	class Triangle : public baseTriangle
	{
	public:
		union_T_xyz_p(Point3D, a, b, c, p);

		Triangle(const Point3D& A, const Point3D& B, const Point3D& C) :baseTriangle(Vector3D::cross(B - A, C - A)), a(A), b(B), c(C){}

		const Point3D& operator[](size_t ind) const
		{
			assert(ind < 3);
			return p[ind];
		}
		Point3D& operator[](size_t ind)
		{
			assert(ind < 3);
			return p[ind];
		}

		const Point3D& at(size_t ind) const
		{
			assert(ind < 3);
			return p[ind];
		}
		Point3D& at(size_t ind)
		{
			assert(ind < 3);
			return p[ind];
		}

		std::pair<bool, Segment3D> Intersection(const Plane& P) const;
		double area() const;
	};

	class QSegment3D :public Segment3D
	{
	public:
		size_t tind1;
		size_t tind2;

		QSegment3D(const Point3D& PStart, const Point3D& PEnd, const size_t i1 = 0, const size_t i2 = 0) : Segment3D(PStart, PEnd), tind1(i1), tind2(i2)
		{
		}
	};

	class Cylinder
	{
	public:
		static const unsigned int nedge = 36;
		Point3D cen;
		uVector3D n;
		double radius;
		double height;

		Cylinder()
		{
		}
		Cylinder(const Point3D& cen_, const uVector3D& n_, double r = 2.5, double h = 10);

		bool isinside(const Point3D& p) const
		{
			auto v = p - cen;
			double h = Vector3D::dot(v, n);

			return h <= height && v.lengthsq() - h*h <= radius*radius;
		}

		std::vector<std::shared_ptr<baseTriangle>> Triangulation(bool twoside) const;
	};



	class TransformationMatrix
	{
	private:
		void Init()
		{
			InitR();
			InitT();
		}
		void InitR()
		{
			memset(r, 0, sizeof(r));

			r[0] = r[5] = r[10] = r[15] = 1;
		}
		void InitT()
		{
			t = {0, 0, 0};
		}
	public:
		double r[16];
		Point3D t;

		TransformationMatrix()
		{
			Init();
		}
		TransformationMatrix(uVector3D x, uVector3D y, uVector3D z)
		{
			Init();

			r[0] = x.x, r[1] = y.x, r[2] = z.x;
			r[4] = x.y, r[5] = y.y, r[6] = z.y;
			r[8] = x.z, r[9] = y.z, r[10] = z.z;
		}

		Vector3D getXaxis() const
		{
			return {r[0], r[4], r[8]};
		}
		Vector3D getYaxis() const
		{
			return {r[1], r[5], r[9]};
		}
		Vector3D getZaxis() const
		{
			return {r[2], r[6], r[10]};
		}

		static Vector3D getXaxis(const double r[16])
		{
			return {r[0], r[4], r[8]};
		}
		static Vector3D getYaxis(const double r[16])
		{
			return {r[1], r[5], r[9]};
		}
		static Vector3D getZaxis(const double r[16])
		{
			return {r[2], r[6], r[10]};
		}
		static Point3D getTranslation(const double r[16])
		{
			return {r[12], r[13], r[14]};
		}

		Point3D applyMat(const Point3D& p) const
		{
			return t + Point3D {Vector3D::dot(getXaxis(), p), Vector3D::dot(getYaxis(), p), Vector3D::dot(getZaxis(), p)};
		}
		Point3D removeMat(const Point3D& p) const
		{
			Point3D pmc = p + t*-1;

			return {Vector3D::dot({r[0], r[1], r[2]}, pmc), Vector3D::dot({r[4], r[5], r[6]}, pmc), Vector3D::dot({r[8], r[9], r[10]}, pmc)};
		}
		Vector3D applyMat(const Vector3D& v) const
		{
			return {Vector3D::dot(getXaxis(), v), Vector3D::dot(getYaxis(), v), Vector3D::dot(getZaxis(), v)};
		}
		Vector3D removeMat(const Vector3D& v) const
		{
			return {Vector3D::dot({r[0], r[1], r[2]}, v), Vector3D::dot({r[4], r[5], r[6]}, v), Vector3D::dot({r[8], r[9], r[10]}, v)};
		}
		Triangle applyMat(const baseTriangle& t) const
		{
			return Triangle(applyMat(t.at(0)), applyMat(t.at(1)), applyMat(t.at(2)));
		}
		Triangle removeMat(const baseTriangle& t) const
		{
			return Triangle(removeMat(t.at(0)), removeMat(t.at(1)), removeMat(t.at(2)));
		}
		TransformationMatrix applyMat(const TransformationMatrix& tm) const
		{
			auto p = applyMat(tm.t);
			auto vx = applyMat(tm.getXaxis());
			auto vy = applyMat(tm.getYaxis());
			auto vz = applyMat(tm.getZaxis());

			TransformationMatrix tm2;
			tm2.r[0] = vx.x;
			tm2.r[1] = vy.x;
			tm2.r[2] = vz.x;
			tm2.r[4] = vx.y;
			tm2.r[5] = vy.y;
			tm2.r[6] = vz.y;
			tm2.r[8] = vx.z;
			tm2.r[9] = vy.z;
			tm2.r[10] = vz.z;
			tm2.t = p;
			return tm2;
		}
		TransformationMatrix removeMat(const TransformationMatrix& tm) const
		{
			auto p = removeMat(tm.t);
			auto vx = removeMat(tm.getXaxis());
			auto vy = removeMat(tm.getYaxis());
			auto vz = removeMat(tm.getZaxis());

			TransformationMatrix tm2;
			tm2.r[0] = vx.x;
			tm2.r[1] = vy.x;
			tm2.r[2] = vz.x;
			tm2.r[4] = vx.y;
			tm2.r[5] = vy.y;
			tm2.r[6] = vz.y;
			tm2.r[8] = vx.z;
			tm2.r[9] = vy.z;
			tm2.r[10] = vz.z;
			tm2.t = p;
			return tm2;
		}

		static TransformationMatrix genmat(double angle, double x, double y, double z)
		{
			TransformationMatrix tm;
			auto result = tm.r;

			double l = 1.0 / std::sqrt(x*x + y*y + z*z);
			double xx = x*l;
			double yy = y*l;
			double zz = z*l;

			double dtor = M_PI / 180;
			double ct = std::cos(angle * dtor);
			double st = std::sin(angle * dtor);
			double ctt = 1 - ct;

			result[0] = ct + xx*xx*ctt;
			result[1] = xx*yy*ctt + zz*st;
			result[2] = xx*zz*ctt - yy*st;
			result[4] = xx*yy*ctt - zz*st;
			result[5] = ct + yy*yy*ctt;
			result[6] = yy*zz*ctt + xx*st;
			result[8] = xx*zz*ctt + yy*st;
			result[9] = yy*zz*ctt - xx*st;
			result[10] = ct + zz*zz*ctt;

			return tm;
		}
		static std::array<double, 16> getMat(const double r[16], const double t[3], const double c[3])
		{
			std::array<double, 16> arr;
			memcpy(arr.data(), r, sizeof(arr));
			Point3D cen {c[0], c[1], c[2]};
			arr[12] = t[0] + cen.x - Vector3D::dot({r[0], r[4], r[8]}, cen);
			arr[13] = t[1] + cen.y - Vector3D::dot({r[1], r[5], r[9]}, cen);
			arr[14] = t[2] + cen.z - Vector3D::dot({r[2], r[6], r[10]}, cen);

			return arr;
		}
		std::array<double, 16> getMat(const Point3D& cen = Point3D{0, 0, 0}) const
		{
			//-c*r + c + t;
			std::array<double, 16> arr;
			memcpy(arr.data(), r, sizeof(arr));
			auto tt = applyMat(cen*-1) + cen;
			arr[12] = tt.x;
			arr[13] = tt.y;
			arr[14] = tt.z;

			return arr;
		}
		void setMat(double(&mat)[16])
		{
			memcpy(r, mat, sizeof(mat));
			t[0] = r[12]; r[12] = 0;
			t[1] = r[13]; r[13] = 0;
			t[2] = r[14]; r[14] = 0;
		}

		void LoadMat(std::string str);
		void SaveMat(std::string str) const;
		void ResetMat()
		{
			Init();
		}
		void ResetMatR()
		{
			InitR();
		}
		void ResetMatT()
		{
			InitT();
		}

		void Transposes()
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = i + 1; j < 4; j++)
				{
					std::swap(r[i * 4 + j], r[j * 4 + i]);
				}
			}
		}

		static TransformationMatrix genmatbyZaxis(uVector3D z)
		{
			uVector3D x = Vector3D::genVVec(z);
			uVector3D y = Vector3D::cross(z, x);
			return TransformationMatrix(x, y, z);
		}
	};

	class Circle
	{

	};

	class Polygon
	{
	public:
		std::vector<Point3D> parr;
		double scale = 1;
		bool isregular = false;

		double xm = DBL_MAX, xM = -DBL_MAX, ym = DBL_MAX, yM = -DBL_MAX;

		void setparr(const std::vector<Point3D>& p)
		{
			isregular = false;
			scale = 1;
			parr = p;
			if (!checkccw())
			{
				std::reverse(parr.begin(), parr.end());
			}
			updatemM();
		}

		void setparr(std::vector<Point3D>&& p)
		{
			isregular = false;
			scale = 1;
			parr.swap(p);
			if (!checkccw())
			{
				std::reverse(parr.begin(), parr.end());
			}
			updatemM();
		}

		bool isinside(const Point3D& p) const
		{
			return isinside(p, scale);
		}

		bool isinside(const Point3D& p, double s) const
		{
			if (isregular && parr.size() > 36)
			{
				return isinside_circle(p, s);
			}

			auto p2 = p / s;
			p2.z = 0;
			double sumangle = 0;
			auto n = parr.size();
			for (size_t i = 0; i < n; i++)
			{
				uVector3D aaa = parr[i] - p2;
				uVector3D bbb = parr[(i + 1) % n] - p2;
				
				double d1 = Vector3D::dot(aaa, bbb);
				d1 = d1 > 0 ? std::min(d1, 1.0) : std::max(d1, -1.0);
				sumangle += std::copysign(acos(d1), Vector3D::cross(aaa, bbb).z);
			}

			return abs(sumangle / M_PI / 2) > 0.99;
		}

		bool isinside_circle(const Point3D& p, double s) const
		{
			auto p2 = p / s;
			
			return hypot(p2.x, p2.y) <= 1;
		}

		void make_regular(size_t n, double radius = 1.0, double offset = 0.0)
		{
			isregular = true;
			parr.clear();
			if (n < 3)
			{
				return;
			}
			if (n > 36)
			{
				make_circle(radius);
				return;
			}

			parr.resize(n);
			double TWOPI = 2.0f * M_PI;
			double delta = TWOPI / n;

			for (size_t count = 0; count < n; count++)
			{
				double d = count * delta + offset;
				parr[count] = {radius * cos(d) , radius * sin(d), 0};
			}
			updatemM();
		}

		void make_circle(double radius)
		{
			parr.resize(40);
			double TWOPI = 2.0f * M_PI;
			double delta = TWOPI / 40;
			for (unsigned int count = 0; count < 40; count++)
			{
				double d = count * delta;
				parr[count] = {radius * cos(d) , radius * sin(d), 0};
			}
			updatemM();
		}

		bool checkccw() const
		{
			double deg = 0;
			auto n = parr.size();
			for (size_t i = 0; i < n; i++)
			{
				uVector3D v1 = parr[(i + 1) % n] - parr[i];
				uVector3D v2 = parr[(i + 2) % n] - parr[(i + 1) % n];
				double d1 = Vector3D::dot(v1, v2);
				d1 = d1 > 0 ? std::min(d1, 1.0) : std::max(d1, -1.0);
				deg += std::copysign(acos(d1), Vector3D::cross(v1, v2).z);
			}

			return deg > 0;
		}

		void updatemM()
		{
			for (auto&& pp : parr)
			{
				xm = std::min(xm, pp.x);
				xM = std::max(xM, pp.x);
				ym = std::min(ym, pp.y);
				yM = std::max(yM, pp.y);
			}
		}
	};	

	double Determinant(const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D);
	double Determinant(const Vector3D& A, const Vector3D& B, const Vector3D& C);
}

namespace std
{
#define ROT32(x, y) ((x << y) | (x >> (32 - y))) // avoid effort
	inline uint32_t murmur_hash3_32(const char *key, uint32_t len, uint32_t seed = 0)
	{
		static const uint32_t c1 = 0xcc9e2d51;
		static const uint32_t c2 = 0x1b873593;
		static const uint32_t r1 = 15;
		static const uint32_t r2 = 13;
		static const uint32_t m = 5;
		static const uint32_t n = 0xe6546b64;

		uint32_t hash = seed;

		const int nblocks = len / 4;
		const uint32_t *blocks = (const uint32_t *)key;
		for (int i = 0; i < nblocks; i++) {
			uint32_t k = blocks[i];
			k *= c1;
			k = ROT32(k, r1);
			k *= c2;

			hash ^= k;
			hash = ROT32(hash, r2) * m + n;
		}

		const uint8_t *tail = (const uint8_t *)(key + nblocks * 4);
		uint32_t k1 = 0;

		switch (len & 3) {
		case 3:
			k1 ^= tail[2] << 16;
		case 2:
			k1 ^= tail[1] << 8;
		case 1:
			k1 ^= tail[0];

			k1 *= c1;
			k1 = ROT32(k1, r1);
			k1 *= c2;
			hash ^= k1;
		}

		hash ^= len;
		hash ^= (hash >> 16);
		hash *= 0x85ebca6b;
		hash ^= (hash >> 13);
		hash *= 0xc2b2ae35;
		hash ^= (hash >> 16);

		return hash;
	}
#undef  ROT32


	template <>
	struct hash<myMath::Point3D>
	{
		size_t operator()(const myMath::Point3D& p) const
		{
			float ff[3];
			for (int i = 0; i < 3; i++)
			{
				ff[i] = p.p[i] == 0 ? 0.0f : abs(p.p[i]) < FLT_EPSILON ? 0.0f : (float)p.p[i];
			}
			
			return murmur_hash3_32((const char*)&ff, sizeof(float) * 3);
			//return _Hash_seq((const unsigned char*)&ff, sizeof(float) * 3);
		}
	};

	template <>
	struct hash<myMath::baseTriangle>
	{
		size_t operator()(const myMath::baseTriangle& t) const
		{
			hash<myMath::Point3D> hp;
			return hash<size_t>()((hp(t.at(0)) ^ hp(t.at(1)) ^ hp(t.at(2))));
		}
	};
	template <>
	struct hash<myMath::baseTriangle*>
	{
		size_t operator()(const myMath::baseTriangle* const t) const
		{
			hash<myMath::baseTriangle> ht;
			return ht(*t);
		}
	};
	template <>
	struct hash<shared_ptr<myMath::baseTriangle>>
	{
		size_t operator()(const std::shared_ptr<myMath::baseTriangle>& t) const
		{
			hash<myMath::baseTriangle> ht;
			return ht(*t);
		}
	};
}
#endif