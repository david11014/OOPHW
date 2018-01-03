#pragma once

#ifndef __STLmodel_H__
#define __STLmodel_H__

#include <bitset>
#include <functional>
#include <set>
#include <map>
#include <unordered_map>

#include "AABBtree.h"
#include "OctCloud.h"
#include "OBJrender.h"
#include "myuiobj.h"

void glTranslatedv(const double* v);
void glTranslatedvm(const double* v);
std::string my_to_string(double val, size_t decimal_places = 1);

namespace myModel
{
	using namespace myMath;
	using namespace std;

	static void LoadStlBinary(const string& fileName, tri_ind_t& trianglecount, float*& dest);
	static float* LoadStlBinary(const string& fileName, tri_ind_t& trianglecount);
	static void SaveStlBinary(const string& fileName, const spvector<baseTriangle>& V);

	class ModelTriangle : public baseTriangle
	{
	public:
		union_T_xyz_p(Point3D*, a, b, c, p);
		union_T_xyz_p(Point3DNode*, na, nb, nc, np);
		bool hasEdge;

		ModelTriangle() :baseTriangle({1,0,0}), a(nullptr), b(nullptr), c(nullptr), na(nullptr), nb(nullptr), nc(nullptr), hasEdge(false)
		{
		}

		const Point3D& operator[](size_t ind) const
		{
			assert(ind<3);
			return *p[ind];//*((&a)[ind]);
		}
		Point3D& operator[](size_t ind)
		{
			assert(ind<3);
			return *p[ind];//*((&a)[ind]);
		}
		const Point3D& at(size_t ind) const
		{
			assert(ind<3);
			return *p[ind];
		}
		Point3D& at(size_t ind)
		{
			assert(ind<3);
			return *p[ind];
		}
		void caln()
		{
			n = Vector3D::cross(*b - *a, *c - *a);
		}
		double area() const
		{
			return 0.5*Vector3D::cross(*b - *a, *c - *a).length();
		}

		inline double lab() const
		{
			return (*a - *b).length();
		}
		inline double lbc() const
		{
			return (*b - *c).length();
		}
		inline double lca() const
		{
			return (*c - *a).length();
		}
		inline double dab() const
		{
			return acos(Vector3D::dot(*b - *c, *a - *c) / lbc() / lca());
		}
		inline double dbc() const
		{
			return acos(Vector3D::dot(*b - *a, *c - *a) / lab() / lca());
		}
		inline double dca() const
		{
			return acos(Vector3D::dot(*a - *b, *c - *b) / lab() / lbc());
		}

		bool isintri(const Point3D& P)
		{
			static auto& VC = Vector3D::cross;
			static auto& VD = Vector3D::dot;

			Vector3D pa = *a - P;
			Vector3D pb = *b - P;
			Vector3D pc = *c - P;
		
			return VD(VC(pa, pb), n) > 0 && VD(VC(pb, pc), n) > 0 && VD(VC(pc, pa), n) > 0;
		}
	};

	class SubTriangles
	{
		vector<double> dtri_ind;
	public:
		SubTriangles(const ModelTriangle* T, const spvector<QSegment3D>& S);
		SubTriangles(const ModelTriangle* T, const spvector<ModelTriangle>& S);
		SubTriangles(const ModelTriangle* T, const double edgelength);

		spvector<baseTriangle> Triangulation() const;
	};

	class ModelPCA
	{
	public:
		Point3D O, basis0, basis1, basis2;
		double eigValues[3];
		ModelPCA() {};
		ModelPCA(spvector<ModelTriangle> *tarr, BoundBox *bbox);
		~ModelPCA() {};

		char* ToString()
		{
			stringstream strStream;

			strStream << "O:" << O[0] << "," << O[1] << "," << O[2] << "\n";
			strStream << "eigveter:\n";
			strStream << basis0.x << " " << basis0.y << " " << basis0.z << " eigvalue: " << eigValues[0] << "\n";
			strStream << basis1.x << " " << basis1.y << " " << basis1.z << " eigvalue: " << eigValues[1] << "\n";
			strStream << basis2.x << " " << basis2.y << " " << basis2.z << " eigvalue: " << eigValues[2] << "\n";

			string str = strStream.str();
			std::vector<char> cstr(str.c_str(), str.c_str() + str.size() + 1);
			char* ca = &cstr[0];
			return ca;
		}

		void CalModelPCA(spvector<ModelTriangle> *tarr, BoundBox *bbox);
	};

	class Cells
	{
	public:
		struct data
		{
			double deep;
			sp<ModelTriangle> Tri;
		};

		Point3D O;
		uVector3D direction;		
		vector<data> IntsDatas;
		Cells() {};
		Cells(Point3D o,uVector3D d) : O(o), direction(d){};

		void SortCell() { sort(IntsDatas.begin(), IntsDatas.end(), DataCompare); };

		char* ToString()
		{
			stringstream strStream;

			strStream << "O: " << O.ToString() << " D:" << direction.x << " " << direction.y << " " << direction.z;

			if (IntsDatas.size() > 0)
			{
				strStream << "\n";
				for (auto a : IntsDatas)
				{
					strStream << a.deep << " " << a.Tri << "|";
				}
			}					

			string str = strStream.str();
			std::vector<char> cstr(str.c_str(), str.c_str() + str.size() + 1);
			char* ca = &cstr[0];
			return ca;
		}

	private:
		struct DataCompareObject {
			bool operator()(const data &a, const data &b)
			{
				return a.deep < b.deep;
			}
		}DataCompare;
		
	};

	class LDI
	{
		/** 由XYZ平面發出射線採樣 **/
	public:
		const double D = 1.0; //sample distance

		map<pair<int,int>, sp<Cells>> LDI_Data[3];
		
		int BoxMax[3];
		int BoxMin[3];		
			
		spvector<ModelTriangle> *ModelTr;
		BoundBox *bbbox;

		LDI(BoundBox& BBox, spvector<ModelTriangle>& ModelTr) :bbbox(&BBox), ModelTr(&ModelTr)
		{	
			UpdateLDI();
		};
		~LDI() {};
		void UpdateLDI();		

	private:

		//Moller–Trumbore intersection algorithm
		//code From https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
		bool RayIntersectsTriangle(Vector3D rayOrigin, Vector3D rayVector, ModelTriangle* inTriangle, Point3D& outIntersectionPoint)
		{
			const double EPSILON = 0.0000001;
			Vector3D vertex0 = (Vector3D)(*(inTriangle->a));
			Vector3D vertex1 = (Vector3D)(*(inTriangle->b));
			Vector3D vertex2 = (Vector3D)(*(inTriangle->c));
			Vector3D edge1, edge2, h, s, q;
			float a, f, u, v;
			edge1 = vertex1 - vertex0;
			edge2 = vertex2 - vertex0;
			h = Vector3D::cross(rayVector, edge2);
			a = Vector3D::dot(edge1, h);
			if (a > -EPSILON && a < EPSILON)
				return false;
			f = 1 / a;
			s = rayOrigin - vertex0;
			u = f * (Vector3D::dot(s, h));
			if (u < 0.0 || u > 1.0)
				return false;
			q = Vector3D::cross(s, edge1);
			v = f * Vector3D::dot(rayVector, q);
			if (v < 0.0 || u + v > 1.0)
				return false;
			// At this stage we can compute t to find out where the intersection point is on the line.
			float t = f * Vector3D::dot(edge2, q);
			if (t > EPSILON) // ray intersection
			{
				Vector3D temp = rayOrigin + rayVector * t;
				outIntersectionPoint.x = temp.x;
				outIntersectionPoint.y = temp.y;
				outIntersectionPoint.z = temp.z;
				return true;
			}
			else // This means that there is a line intersection but not a ray intersection.
				return false;
		}

		
	};

	class Model :public VBOobj
	{
	public:
		typedef vector<triid_pid> adjtri;
		typedef map<tri_ind_t, spvector<QSegment3D>> TriCrossResult;
		typedef map<tri_ind_t, spvector<baseTriangle>> TriSeparateResult;
		
		struct triangleinfo
		{
			adjtri a;
			adjtri b;
			unsigned char c;
		};//28
		struct Edgeind
		{
			tri_ind_t a;
			vertex_ind_t b;
			vertex_ind_t c;
		};
		
		double colorR, colorG, colorB;
		
		TransformationMatrix tm;

		double MinTriSize;
		double AvgTriSize;
		double MaxTriSize;
		BoundBox bbox;
		sp<OctCloud> ocd;
		sp<AABBtree> ABt;
		LDI *ldi;
		ModelPCA PAxis;


		bool showmesh;
		bool dynamicshowmesh;
	private:
		vector<Point3D> sourceP;
	public:
		spvector<ModelTriangle> tarr;
		spvector<triangleinfo> tarrinfo;

		Model();
		Model(const Model& M);
		Model(Model&& M);
		Model& operator=(const Model& M) = delete;
		virtual ~Model();

		explicit Model(const spvector<baseTriangle>& V);
		explicit Model(const string& fileName);
		
		static Model LoadObjASCI(const string& fileName);
	private:
		void SetModel(float* xyz, tri_ind_t numofTri);
	public:
		void SaveStlBinary(const string& fileName) const;
		void SavePointCloudASCII(const string& fileName) const;
		void SavePointCloudBinary(const string& fileName) const;
		const ModelTriangle& operator[](size_t ind) const
		{
			if (ind < tarr.size())
			{
				return *(tarr[ind]); 
			}
			return *(tarr[0]);
		}
		ModelTriangle& operator[](size_t ind)
		{
			if (ind < tarr.size())
			{
				return *(tarr[ind]); 
			}
			return *(tarr[0]);
		}
		
		void Add_tri(const spvector<baseTriangle>& bt);

		void ApplyTMat();
		void RebuildOctCloud();
		void Triangle_update();
		void SetMinTriSize();
		void BuildAABB();

		tuple<bool, char, char, double> checkDegenerationTriangle(const ModelTriangle& T, double llimit) const;
		spvector<baseTriangle> searchDegenerationTriangle(double llimit = 1e-6) const;
		void removeDegenerationTriangle(double llimit = 1e-6);

		static pair<bool, QSegment3D*> checkTriangleCross(const ModelTriangle& T1, const ModelTriangle& T2);
		static pair<bool, QSegment3D*> checkTriangleCross_onesame(const ModelTriangle& T1, const ModelTriangle& T2);

		TriCrossResult searchTriangleCross();
		TriCrossResult searchTriangleCross(const Model* rhs, const inddatas& testind) const;
		TriSeparateResult separateTriangle(const TriCrossResult& x);

		static spvector<Segment3D> searchTriangleCross_fordisplay(const TriCrossResult& x);
		static spvectors<Segment3D> searchTriangleCross_fordisplay2(const TriCrossResult& x);
		static spvector<baseTriangle> subtri_fordisplay(const TriSeparateResult& x);
		spvectors<Segment3D> searchEdgecontour_fordisplay() const;
		static spvectors<Segment3D> searchEdgecontour_fordisplay(const Model& m, const vectors<Edgeind>& vei);
		static void repairtsr(vector<reference_wrapper<spvector<ModelTriangle>>>& modelarr, vector<Model::TriSeparateResult>& tsr);

		spvector<baseTriangle> searchEdgeTriangle() const;
		spvector<Segment3D> searchEdge() const;
		vector<Edgeind> searchEdgeind() const;
		vectors<Edgeind> searchEdgecontour() const;
		vectors<Edgeind> searchEdgecontour(vector<Edgeind>& sourceEdge) const;

		sp<baseTriangle> searchTriangle(const Point3D& P) const;
		tri_ind_t searchTriangleind(const Point3D& P) const;
		spvector<baseTriangle> getTriangles(const inddata& ind) const;
		spvectors<baseTriangle> getTriangles(const inddatas& ind) const;

		double calcMeanCurvature(const Point3DNode* pn) const;
		double calcGaussCurvature(const Point3DNode* pn) const;

		std::pair<double, double> calcCurvature(const Point3DNode* pn) const;
		bool signCurvature(const Point3DNode* pn) const;
		double calcCurvature(const Point3D& P) const;

		Point3D PointinModelCoor(const Point3D& P) const;
		Vector3D VectorinModelCoor(const Vector3D& V) const;
		Point3D PointinWorldCoor(const Point3D& P) const;
		Vector3D VectorinWorldCoor(const Vector3D& V) const;

		spvector<baseTriangle> shiftTriangle(double D, unsigned int dir) const;
		spvector<baseTriangle> shiftTriangle2(double D, unsigned int dir, Vector3D n) const;
		spvector<baseTriangle> shiftTriangle3(const map<Point3DNode*, double>& mpnd, double shiftd, unsigned int dir) const;
		spvector<baseTriangle> shiftTriangle4(const map<Point3DNode*, double>& mpnd, double D, unsigned int dir, Vector3D n) const;

		pair<Point3D, Point3D> Intersection(const Line3D& i) const;
		sp<Model> filpmodel() const;

		void modifyvbo(const inddata& ind) const;

		sp<Model> refine() const;
		sp<Model> refine(double edgelength) const;

		Point3D find_biggest_tri_cen_without_self_cross();
	};

	class ModelBrush
	{
	public:
		TransformationMatrix tm;
		Prism pm;

		void settm(double(&modelview)[16])
		{
			tm.setMat(modelview);
		}
		void setpm_cen(const Point3D& p)
		{
			pm.cen = p;
		}
		void setpm_poly(const vector<Point3D>& parr)
		{
			pm.poly.setparr(parr);
		}
		void setpm_poly(vector<Point3D>&& parr)
		{
			pm.poly.setparr(std::move(parr));
		}
		void setpm_height(double h)
		{
			pm.height = h;
		}
// 		void setpm_twoside(bool b)
// 		{
// 			pm.twoside = b;
// 		}

		bool isinside(const Point3D& p) const
		{
			return pm.isinside(tm.applyMat(p-(Vector3D)pm.cen));
		}
		bool isinside(const baseTriangle& bt) const
		{
			return isinside(bt[0]) && isinside(bt[1]) && isinside(bt[2]);
		}
		bool iscross(const baseTriangle& bt) const
		{
			return isinside(bt[0]) || isinside(bt[1]) || isinside(bt[2]);
		}
		bool isinside_s(const Point3D& p) const
		{
			return (p - pm.cen).length() / pm.poly.scale < 1;
		}
		bool iscross_s(const baseTriangle& bt) const
		{
			return isinside_s(bt[0]) || isinside_s(bt[1]) || isinside_s(bt[2]);
		}
	};

	class ModelSmooth
	{
		vector<Point3D> parr;
		sp<Model> mmodel;
		sp<Model> result;
		size_t iter = 5;

		template<typename algo, typename... arg>
		void applyalgo(algo a, arg... p)
		{
			applyalgo(result->ocd->pnarr, parr, a, p...);
		}

		void update()
		{
			update(result->ocd->pnarr, parr);
		}
	public:
		void setmodel(sp<Model> mmodel_)
		{
			mmodel = mmodel_;
			result = make_shared<Model>(*mmodel);
			parr.resize(result->ocd->pnarr.size());
		}

		template<typename algo, typename... arg>
		static void applyalgo(const vector<Point3DNode*> pnarr, vector<Point3D>& parr, algo a, arg... p)
		{
			for (size_t i = 0; i < pnarr.size(); i++)
			{
				parr[i] = a(pnarr[i], p...);
			}
		}
		static void update(const vector<Point3DNode*> pnarr, const vector<Point3D>& parr)
		{
			for (size_t i = 0; i < pnarr.size(); i++)
			{
				*pnarr[i]->pdata = parr[i];
			}
		}

		static Point3D Laplacian_average(const Point3DNode* pn)
		{
			Point3D p = {0,0,0};
			for (auto&& adjpn : pn->adj_Point)
			{
				p += *adjpn->pdata;
			}
			return p / (double)pn->adj_Point.size();
		}
		static Point3D ScaleDependentLaplacian_average(const Point3DNode* pn, double fac)
		{
			Point3D p = {0,0,0};
			double w = 0;
			for (auto&& adjpn : pn->adj_Point)
			{
				double tempw = 1.0 / (*pn->pdata - *adjpn->pdata).length();
				p += (*adjpn->pdata - *pn->pdata) * tempw;
				w += tempw;
			}
			if (isinf(w))
			{
				p = {0,0,0};
			}
			else
			{
				p /= w;
			}
			return *pn->pdata + p*fac;
		};
		static double cott(sp<Model> m, const Point3DNode* pn, const Point3DNode* pn2)
		{
			auto a = find_if(pn->adj_tri.begin(), pn->adj_tri.end(), [&](const triid_pid& ind)
			{
				return m->tarr[ind.tid()]->np[(ind.pid() + 2) % 3] == pn2;
			});
			auto b = find_if(pn->adj_tri.begin(), pn->adj_tri.end(), [&](const triid_pid& ind)
			{
				return m->tarr[ind.tid()]->np[(ind.pid() + 1) % 3] == pn2;
			});

			auto& t1 = m->tarr[a->tid()];
			auto& t2 = m->tarr[b->tid()];

			auto v11 = t1->at((a->pid() + 1) % 3) - *pn->pdata;
			auto v12 = t1->at((a->pid() + 1) % 3) - *pn2->pdata;

			auto v21 = t2->at((b->pid() + 2) % 3) - *pn2->pdata;
			auto v22 = t2->at((b->pid() + 2) % 3) - *pn->pdata;

			return Vector3D::dot(v11, v12) / Vector3D::cross(v12, v11).length() + Vector3D::dot(v21, v22) / Vector3D::cross(v22, v21).length();
		};
		static Point3D CurvatureFlow_average(const Point3DNode* pn, sp<Model> m)
		{
			if (pn->adj_Point.size() != pn->adj_tri.size())
			{
				return *pn->pdata;
			}

			Point3D p = {0,0,0};
			double w = 0;
			for (auto&& adjpn : pn->adj_Point)
			{
				double tempw = cott(m, pn, adjpn);
				p += (*adjpn->pdata - *pn->pdata) * tempw;
				w += tempw;
			}
			
			p /= w;
			
			return *pn->pdata + p;
		};
		static Point3D Taubin_average(const Point3DNode* pn, double fac)
		{
			Point3D p = {0,0,0};
			for (auto&& adjpn : pn->adj_Point)
			{
				p += (*adjpn->pdata - *pn->pdata);
			}
			p /= (double)pn->adj_Point.size();
			return *pn->pdata + p*fac;
		};
		static Point3D HCLaplacian_average(const Point3DNode* pn, Point3D& srcp, double fac1)
		{
			Point3D p = {0,0,0};
			for (auto&& adjpn : pn->adj_Point)
			{
				p += *adjpn->pdata;
			}
			p /= (double)pn->adj_Point.size();
			auto tempp = p - Vector3D(srcp * fac1 + *pn->pdata * (1 - fac1));

			srcp = p;
			return tempp;
		};
		static Point3D HCLaplacian_average2(const Point3DNode* pn, Point3D& srcp, double fac2)
		{
			Point3D p = {0,0,0};
			for (auto&& adjpn : pn->adj_Point)
			{
				p += *adjpn->pdata;
			}
			p *= (1 - fac2) / pn->adj_Point.size();

			return srcp - Vector3D(*pn->pdata * fac2 + p);
		};

		sp<Model> Laplacian()
		{
			applyalgo(Laplacian_average);
			update();

			return result;
		}
		sp<Model> ScaleDependentLaplacian(double lda = 0.6307)
		{
			applyalgo(ScaleDependentLaplacian_average, lda);
			update();

			return result;
		}
		sp<Model> CurvatureFlow()
		{
			for (size_t t = 0; t < iter; t++)
			{
				applyalgo(CurvatureFlow_average, result);
				update();
			}

			return result;
		}
		sp<Model> Taubin(double lda = 0.6307, double kpb = 0.1)
		{

			if (kpb <= 0)
			{
				cout << "need \"kpb > 0\"\n";
				return result;
			}
			double mu = 1 / (kpb - 1 / lda);

			for (size_t t = 0; t < iter; t++)
			{
				applyalgo(Taubin_average, lda);
				update();
				applyalgo(Taubin_average, mu);
				update();
			}

			return result;
		}
		sp<Model> HCLaplacian(double fac1 = 0.1, double fac2 = 0.6)
		{
			vector<Point3D> srcparr(result->ocd->pnarr.size());
			transform(result->ocd->pnarr.begin(), result->ocd->pnarr.end(), srcparr.begin(), [](auto&& pn) {return *pn->pdata; });
			
			for (size_t t = 0; t < iter; t++)
			{
				for (size_t i = 0; i < result->ocd->pnarr.size(); i++)
				{
					parr[i] = HCLaplacian_average(result->ocd->pnarr[i], srcparr[i], fac1);
				}
				update();
				for (size_t i = 0; i < result->ocd->pnarr.size(); i++)
				{
					parr[i] = HCLaplacian_average2(result->ocd->pnarr[i], srcparr[i], fac2);
				}
				update();
			}

			return result;
		}
	};

	class ModelSplitPlane
	{
	public:
		sp<Model> mmodel;
		spvectors<Segment3D> sarrs;

		ModelSplitPlane(sp<Model> m) :mmodel(m)
		{
			sarrs = mmodel->searchEdgecontour_fordisplay();
		}
	};

	class ModelCutter
	{
	public:
		Point3D lastPoint;
		vector<Point3D> parr;
		vector<Point3D> sparr;
		int setstep;
		
		ModelCutter() :parr(), sparr(), setstep(1)/*, splitplane(nullptr)*/ {}

		void add(double p_dis_in_pixel = 0)
		{
			if (parr.size() == 0 || parr.back().distance(lastPoint) > p_dis_in_pixel)
			{
				parr.push_back(lastPoint);
			}
		}
		void reset()
		{
			setstep = 1;
			parr.clear();
			sparr.clear();
		}

		sp<Model> cutModel(const Model& M, double(&mdvm)[16]);
		sp<ModelSplitPlane> generateplane(double(&mdvm)[16], const Point3D& mz, const Point3D& Mz, double minl) const
		{
			spvector<baseTriangle> tarr;
			if (sparr.size() < 3)
			{
				return nullptr;
			}

			auto ps = sparr.size();
			
			uVector3D zaxis = Vector3D({mdvm[2], mdvm[6], mdvm[10]});
			const size_t zdirc = min(50, int(Vector3D::dot(Mz - mz, zaxis) / minl / 10) + 1);
			
			Plane pl1(zaxis, mz);
			pl1.d -= 2;
			double zl = (pl1.distance(Mz)+2) / zdirc;

			vector<Point3D> transp(ps);
			Vector3D vx {mdvm[0], mdvm[1], mdvm[2]};
			Vector3D vy {mdvm[4], mdvm[5], mdvm[6]};
			Vector3D vz {mdvm[8], mdvm[9], mdvm[10]};
			Point3D tp {mdvm[12], mdvm[13], mdvm[14]};
			for (size_t i = 0; i < ps; i++)
			{
				auto tempp = sparr[i] - tp;
				transp[i].x = Vector3D::dot(vx, tempp);
				transp[i].y = Vector3D::dot(vy, tempp);
				transp[i].z = Vector3D::dot(vz, tempp);
			}

			transp.push_back(transp[0]);

			for (size_t i = 0; i < ps; i++)
			{
				Point3D startp = pl1.Projection(transp[i]);
				Point3D endp = pl1.Projection(transp[(i + 1)]);
				
				Vector3D l = endp - startp;
				const size_t ldirc = min(50, int(l.length() / minl / 10) + 1);
				
				l *= 1.0 / ldirc;
				vectors<Point3D> pp;

				pp.resize(ldirc+1);
				
				size_t a = 0;
				double b = 1;

				for (size_t j = 0; j < ldirc; j++)
				{
					pp[j].push_back(startp + l*(double)j);
					for (size_t k = 0; k < zdirc; k++)
					{
						pp[j].push_back(pp[j][0] + zaxis*zl*(k + 1.0));
					}
				}

				pp[ldirc].push_back(endp);
				for (size_t k = 0; k < zdirc; k++)
				{					
					pp[ldirc].push_back(pp[ldirc][0] + zaxis*zl*(k + 1.0));
				}

				for (size_t j = 0; j < ldirc; j++)
				{
					for (size_t k = 0; k < zdirc; k++)
					{
						tarr.push_back(make_shared<Triangle>(pp[j][k], pp[j + 1][k], pp[j][k + 1]));
						tarr.push_back(make_shared<Triangle>(pp[j + 1][k], pp[j + 1][k + 1], pp[j][k + 1]));
					}
				}
			}
			return make_shared<ModelSplitPlane>(make_shared<Model>(tarr));
		}
	};

	class ModelManger
	{
	public:
		static const unsigned int n = 9;
		shared_ptr<Model> mmodel;
		string modelfilefullname;
		string modelfilepath;
		string modelname;
		string modelnamewtoext;
		array<OBJrender::test2::layerbase<world_layer>*, n> renderarr;
		bitset<n> renderflags;
		bitset<n> visibleflags;

		bool isshown;
		bool isMatchange;
		int transparent_value;
		~ModelManger();
		ModelManger(const string& filename, bool compucross = true);
		ModelManger(const string& filename, shared_ptr<Model> Mmodel, bool compucross = true);
		
		bool istransparent() const
		{
			return transparent_value != 10;
		}
		void initialthis(bool compucross = true);
		
		void updateMat();
		void modifyname(const string& name);
		void setModelcolor(double r, double g, double b)
		{
			mmodel->colorR = r;
			mmodel->colorG = g;
			mmodel->colorB = b;
		}

		void reinitial(bool compucross = true);
		void updatevbo(const inddata& ind);
		void releaseVBO();
	};
	
	class ModelSelect
	{
	public:
		shared_ptr<Model> mmodel;
		pair<vector<unsigned int>, Point3D> tid;
		vector<unsigned int> checktid;

		void setmodel(const shared_ptr<Model>& sm)
		{
			mmodel = sm;

			tid.first.clear();
			checktid.resize(mmodel->tarr.size());

			unsigned int c = 0;
			generate(checktid.begin(), checktid.end(), [&c]{ return c++; });
		}

		void setcolor(double r, double g, double b)
		{
			tid.second = {r, g, b};
		}

		void addselect(const Point3D& p, const Vector3D& dir)
		{
			for (auto&& i : checktid)
			{
				auto& t = mmodel->tarr.at(i);
				if (Vector3D::dot(t->n, dir) > 0)
				{
					for (unsigned int j = 0; j < 3; j++)
					{
						auto ap = t->at(j) - p;
						double zz = abs(Vector3D::dot(ap, dir));
						if (zz < 2 && ap.lengthsq() - zz*zz < 25)
						{
							tid.first.push_back(i * 3);
							tid.first.push_back(i * 3 + 1);
							tid.first.push_back(i * 3 + 2);
							i = -1;
							break;
						}
					}
				}
			}
			checktid.erase(remove_if(checktid.begin(), checktid.end(), [](unsigned int ind){ return ind == -1; }), checktid.end());
		}

	};

	class ColorHistogram :public VBOobj
	{
	public:
		float center;
		float halfrange;
		float Saturation;
		float Value;

		ColorHistogram(float C, float H) :center(C), halfrange(H), Saturation(1.0f), Value(1.0){}
		virtual ~ColorHistogram(){};

		void fontall(float x, float y) const;
		array<float, 3> getRGB(float p) const
		{
			float temp = p - center;

			float temp2 = min(1.0f, abs(temp) / halfrange);//-1~1

			temp2 = copysign(temp2, temp)*150.0f + 150.0f;//0~300
			int ccc = (int)(temp2) / 60;//0~4
			float X = (1 - abs(fmod(temp2 / 60.0f, 2.0f) - 1));

			float rr, gg, bb;
			switch (ccc)
			{
			case 0:
				rr = 1.0f, gg = X, bb = 0.0f;
				break;
			case 1:
				rr = X, gg = 1.0f, bb = 0.0f;
				break;
			case 2:
				rr = 0.0f, gg = 1.0f, bb = X;
				break;
			case 3:
				rr = 0.0f, gg = X, bb = 1.0f;
				break;
			case 4:
				rr = X, gg = 0.0f, bb = 1.0f;
				break;
			case 5:
				rr = 1.0f, gg = 0.0f, bb = X;
				break;
			default:
				rr = 0.0f, gg = 0.0f, bb = 0.0f;
				break;
			}

			return {rr, gg, bb};
		}

		double getHue(const Point3D& colorp) const
		{
			auto m = std::minmax({colorp.x, colorp.y, colorp.z});
			double deltap = m.second - m.first;


			if (m.second == 0 || deltap == 0)
			{
				return 0;
			}
			double h;
			if (colorp.x == m.second)
			{
				h = (colorp.y - colorp.z) / deltap;
			}
			else if (colorp.y == m.second)
			{
				h = 2 + (colorp.z - colorp.x) / deltap;
			}
			else
			{
				h = 4 + (colorp.x - colorp.y) / deltap;
			}

			h *= 60;

			if (h < 0)
			{
				h += 360;
			}

			return h;
		}
		double getD(const Point3D& colorp) const
		{
			auto h = getHue(colorp);

			return (min(h, 300.0) - 150) / 150 * halfrange + center;
		}
	};

	class ColorHistogram2 :public VBOobj, public OBJrender::test2::moveableui
	{
	public:
		float lowlim;
		float hilim;
		int currtc;
		int changelim;

		ColorHistogram2(float l, float h) :moveableui({55.0, 55.0}), lowlim(l), hilim(h), currtc(0), changelim(0)
		{
			if (hilim < 0)
			{
				hilim = 0;
			}
			if (lowlim > hilim)
			{
				lowlim = hilim;
			}
		}
		virtual ~ColorHistogram2(){};

		static std::array<float, 3> hsv2rgb(float h, float s, float v)
		{
			h = 300 - h;//for fjj
			int hi = (int)(h) / 60;//0~4
			float f = h / 60 - hi;

			float p = v*(1 - s);
			float q = v*(1 - f*s);
			float t = v*(1 - (1 - f)*s);

			float rr, gg, bb;
			switch (hi)
			{
			case 0:
				rr = v, gg = t, bb = p;
				break;
			case 1:
				rr = q, gg = v, bb = p;
				break;
			case 2:
				rr = p, gg = v, bb = t;
				break;
			case 3:
				rr = p, gg = q, bb = v;
				break;
			case 4:
				rr = t, gg = p, bb = v;
				break;
			case 5:
				rr = v, gg = p, bb = q;
				break;
			default:
				rr = 0.0f, gg = 0.0f, bb = 0.0f;
				break;
			}

			return {rr, gg, bb};
		}
		static std::array<float, 3> rgb2hsv(float r, float g, float b)
		{
			auto m = std::minmax({r, g, b});
			float deltap = m.second - m.first;

			if (m.second == 0 || deltap == 0)
			{
				return {0, 0, 0};
			}

			float h;
			if (r == m.second)
			{
				h = (g - b) / deltap;
			}
			else if (g == m.second)
			{
				h = 2 + (b - r) / deltap;
			}
			else
			{
				h = 4 + (r - g) / deltap;
			}

			h *= 60;

			if (h < 0)
			{
				h += 360;
			}

			h = 300 - h;//for fjj

			return {h, 1 - m.first / m.second, m.second};
		}

		void fontall(float x, float y) const;
		void fontcurrt(float x, float y) const;
		array<float, 3> getRGB(float p) const
		{
			float temp = p - (hilim + lowlim) / 2;

			float temp2 = min(1.0f, abs(temp) / ((hilim - lowlim) / 2));//-1~1

			temp2 = copysign(temp2, temp)*150.0f + 150.0f;//0~300

			return hsv2rgb(temp2, 1, 1);
		}
		array<float, 4> getRGBA(float p, float A = 1.0f)
		{
			auto RGB = getRGB(p);
			return {RGB[0], RGB[1], RGB[2], A};
		}
		array<float, 3> getcurrtc() const
		{
			return hsv2rgb((float)currtc, 1.0f, 1.0f);
		}
		array<float, 3> getcurrtcc() const
		{
			return hsv2rgb((float)((currtc+180)%360), 1.0f, 1.0f);
		}
		double getcurrtd() const
		{
			return currtc / 300.0*(hilim - lowlim) + lowlim;
		}

		double getHue(const array<float, 3>& colorp) const
		{
			return rgb2hsv(colorp[0], colorp[1], colorp[2])[0];
		}
		double getHue(const Point3D& colorp) const
		{
			return getHue(colorp.tofarr());
		}
		double getD(const array<float, 3>& colorp) const
		{
			auto h = getHue(colorp);

			return (min(h, 300.0) - 150) / 150 * (hilim - lowlim) / 2 + (hilim + lowlim) / 2;
		}
		double getD(const Point3D& colorp) const
		{
			return getD(colorp.tofarr());
		}

		void setcurrtc(int c)
		{
			currtc = c;
			if (currtc < 0)
			{
				currtc = 0;
			}
			else if (currtc > 300)
			{
				currtc = 300;
			}
		}
		void addcurrtc(int c)
		{
			setcurrtc(currtc + c);
		}

		bool checkchangelim(int x, int y)
		{
			int xx = x - 55;
			int yy = y - 55;

			int temp = changelim;
			changelim = 0;
			if (xx < 0)
			{
				if (yy < 315 && yy > 0)
				{
					if (yy > 300 && (xx > -9 * (int(log10(hilim)) + 3)))
					{
						changelim = 2;
					}
					else if (yy < 15 && (xx > -9 * (max(0, int(log10(lowlim))) + 3)))
					{
						changelim = 1;
					}
				}
			}

			return changelim != temp;
		}
		bool checkoncolor(int x, int y)
		{
			int xx = x - 55;
			int yy = y - 55;

			if (xx > 0 && xx < 20)
			{
				if (yy < 315 && yy > 0)
				{
					return true;
				}
			}
			
			return false;
		}

		void setlim(float c)
		{
			if (changelim == 1)
			{
				setlowlim(c);
			}
			else if (changelim == 2)
			{
				sethilim(c);
			}
		}
		void addlim(float c)
		{
			if (changelim == 1)
			{
				addlowlim(c);
			}
			else if (changelim == 2)
			{
				addhilim(c);
			}
		}

		void sethilim(float c)
		{
			hilim = max(c, lowlim + 0.1f);
			if (hilim < 0)
			{
				hilim = 0;
			}
			else if (hilim > 100)
			{
				hilim = 100;
			}
		}
		void addhilim(float c)
		{
			sethilim(hilim + c);
		}

		void setlowlim(float c)
		{
			lowlim = min(c, hilim - 0.1f);
			if (lowlim < 0)
			{
				lowlim = 0;
			}
			else if (lowlim > 300)
			{
				lowlim = 300;
			}
		}
		void addlowlim(float c)
		{
			setlowlim(lowlim + c);
		}
	};

	class RegionProperty
	{
	public:
		double color;
		double thinkness;
		double allowance;		
				
		RegionProperty(double t = 0, double a = 0, double c = 0) : thinkness(t), allowance(a), color(c) {};
		static RegionProperty MakeProperty(double t, double a, double c)
		{
			RegionProperty r(t, a, c);
			return r;
		}
		void SetProperty(double t, double a, double c)
		{
			thinkness = t;
			allowance = a;
			color = c;
		}
	};

	class ModelSelect2
	{
	public:
		typedef double vt;
		shared_ptr<Model> mmodel;
		shared_ptr<ColorHistogram2> colorht;
		map<tri_ind_t, vt> tcid;
		map<vt, set<tri_ind_t>> ctid;
		double radius, x, y, s;

		ModelSelect2() :mmodel(nullptr), radius(5), x(0), y(0), colorht(nullptr){}

		void setmodel(const shared_ptr<Model>& sm)
		{
			if (mmodel == sm)
			{
				return;
			}

			mmodel = sm;
			if (mmodel->ABt == nullptr)
			{
				mmodel->BuildAABB();
			}
			tcid.clear();
			ctid.clear();
		}

		void clear()
		{
			mmodel.reset();

			tcid.clear();
			ctid.clear();
		}

		void setxys(double X, double Y, double S)
		{
			x = X;
			y = Y;
			s = S;
		}

		template<typename Fn>
		void selecttri(const Point3D& p, const Vector3D& dir, Fn f)
		{
			auto pp = mmodel->PointinModelCoor(p);
			auto ddir = mmodel->VectorinModelCoor(dir);
			auto tind = mmodel->ABt->returnind({{pp.x - radius, pp.y - radius, pp.z - radius}, {pp.x + radius, pp.y + radius, pp.z + radius}});

			for (auto&& i : tind)
			{
				auto& t = mmodel->tarr.at(i);
				if (Vector3D::dot(t->n, ddir) > 0)
				{
					for (vertex_ind_t j = 0; j < 3; j++)
					{
						auto ap = t->at(j) - pp;
						double zz = abs(Vector3D::dot(ap, ddir));
						if (ap.lengthsq() - zz*zz < radius*radius)
						{
							if (tcid.find(i) != tcid.end())
							{
								ctid[tcid[i]].erase(i);
							}
							f(i);
							break;
						}
					}
				}
			}
		}

		void addselect(const Point3D& p, const Vector3D& dir)
		{
			vt d = colorht->getcurrtd();
			if (ctid.find(d) == ctid.end())
			{
				ctid.insert(make_pair(d, set<tri_ind_t>()));
			}
			selecttri(p, dir, [&](tri_ind_t i)
			{
				tcid[i] = d;
				ctid[d].insert(i);
			});
		}

		void removeselect(const Point3D& p, const Vector3D& dir)
		{
			selecttri(p, dir, [&](tri_ind_t i)
			{
				tcid.erase(i);
			});
		}

		void swapselect(const Point3D& p)
		{
			auto pp = mmodel->PointinModelCoor(p);
			auto ind = mmodel->searchTriangleind(pp);
			if (ind < mmodel->tarr.size() && tcid.find(ind) != tcid.end())
			{
				vt d = colorht->getcurrtd();
				auto v = tcid[ind];
				if (v != d)
				{
					if (ctid.find(d) == ctid.end())
					{
						ctid.insert(make_pair(d, move(ctid[v])));
					}
					else
					{
						ctid[d].insert(ctid[v].begin(), ctid[v].end());
					}
					ctid.erase(v);
					for (auto ii : ctid[d])
					{
						tcid[ii] = d;
					}
				}				
			}
		}

		map<Point3DNode*, double> get_mpnd(double shiftd, unsigned int level = 10) const;

		bool test(const shared_ptr<Model>& sm) const
		{
			return mmodel != sm || tcid.size() == 0;
		}
	};

	class ModelSelect3
	{
	public:
		typedef double vt;
		shared_ptr<Model> mmodel;
		shared_ptr<ColorHistogram2> colorht;
		map<tri_ind_t, vt> Tr2ID; //三角形對區塊編號對應
		map<vt, set<tri_ind_t>> ID2Tr; //區塊編號對三角形對應		
		map<vt, RegionProperty> ID2Pr; //區塊編號與區塊參數對應
		map<vt, pair<double, double>> atid; //區塊編號裕度厚度對應

		double radius, x, y, s;
		sp<ModelBrush> brush;
		myapp::myDefaultValue *DefaultValue;	
				
		double DefTh, DefAw;
		

		ModelSelect3() :mmodel(nullptr), radius(5), x(0), y(0), colorht(nullptr) {};

		~ModelSelect3()
		{
			mmodel.reset();
			colorht.reset();
			brush.reset();
		}

		void setmodel(const shared_ptr<Model>& sm)
		{
			if (mmodel == sm)
			{
				return;
			}

			mmodel = sm;
			if (mmodel->ABt == nullptr)
			{
				mmodel->BuildAABB();
			}
			Tr2ID.clear();
			ID2Tr.clear();
			ID2Pr.clear();
			atid.clear();
			
			DefTh = DefaultValue->Get(0);
			DefAw = DefaultValue->Get(1);			
			
		}

		void clear()
		{
			mmodel.reset();

			Tr2ID.clear();
			ID2Tr.clear();
			ID2Pr.clear();
			atid.clear();

			DefTh = DefaultValue->Get(0);
			DefAw = DefaultValue->Get(1);
		}

		void setxys(double X, double Y, double S)
		{
			x = X;
			y = Y;
			s = S;
		}

		BoundBox getbrushbbox(const Point3D& p, double(&modelview)[16])
		{
			brush->setpm_cen(mmodel->PointinModelCoor(p));
			TransformationMatrix tm;
			tm.setMat(modelview);
			brush->tm = mmodel->tm.removeMat(tm);
			brush->tm.ResetMatT();

			auto bbarr = brush->pm.getBoundBox().bboxcorner();
			for (auto&& pp : bbarr)
			{
				pp = brush->tm.removeMat(pp) + brush->pm.cen;
			}

			return BoundBox::genbbox(bbarr);
		}

		template<typename Fn, typename Fn2>
		void selecttri(const Point3D& p, double(&modelview)[16], Fn f, Fn2 f2, const bool checkfacefront)
		{
			auto tind = mmodel->ABt->returnind(getbrushbbox(p, modelview));
			Vector3D ddir = brush->tm.getZaxis();
			for (auto&& i : tind)
			{
				auto& t = mmodel->tarr.at(i);
				if (!checkfacefront || Vector3D::dot(t->n, ddir) > 0)
				{
					if (f2(*brush, *t))
					{
						f(i);
					}
				}
			}
		}

		template<typename Fn>
		void selecttri(const Point3D& p, double(&modelview)[16], Fn f, const bool checkfacefront = true)
		{
			selecttri(p, modelview, f, std::mem_fn(&ModelBrush::iscross), checkfacefront);
		}

		inddata smoothselect(const Point3D& p, double(&modelview)[16])
		{
			inddata tind;
			vector<Point3DNode*> pnarr;
			auto sm = [&](tri_ind_t i)
			{
				tind.push_back(i);
				pnarr.push_back(mmodel->tarr[i]->na);
				pnarr.push_back(mmodel->tarr[i]->nb);
				pnarr.push_back(mmodel->tarr[i]->nc);
			};
			selecttri(p, modelview, sm, std::mem_fn(&ModelBrush::iscross_s), false);

			sort(pnarr.begin(), pnarr.end());
			pnarr.erase(remove_if(pnarr.begin(), pnarr.end(), [&](Point3DNode* pn) {return !brush->isinside_s(*pn->pdata); }), pnarr.end());

			vector<Point3D> parr(pnarr.size());
			//ModelSmooth::applyalgo(pnarr, parr, ModelSmooth::ScaleDependentLaplacian_average, 0.6307);
			//ModelSmooth::applyalgo(pnarr, parr, ModelSmooth::CurvatureFlow_average, mmodel);
			ModelSmooth::applyalgo(pnarr, parr, ModelSmooth::Taubin_average, 0.6307);
			ModelSmooth::update(pnarr, parr);
			
			return tind;
		}

		template<typename algo, typename... arg>
		inddata smoothselect(algo a, arg... p)
		{
			if (Tr2ID.empty())
			{
				return inddata();
			}

			inddata tind;
			vector<Point3DNode*> pnarr;
			for (auto&& t : Tr2ID)
			{
				auto i = t.first;
				pnarr.push_back(mmodel->tarr[i]->na);
				pnarr.push_back(mmodel->tarr[i]->nb);
				pnarr.push_back(mmodel->tarr[i]->nc);
			}
			sort(pnarr.begin(), pnarr.end());
			pnarr.erase(unique(pnarr.begin(), pnarr.end()), pnarr.end());

			for (auto&& pn : pnarr)
			{
				for (auto&& adjt : pn->adj_tri)
				{
					tind.push_back(adjt.tid());
				}
			}

			sort(tind.begin(), tind.end());
			tind.erase(unique(tind.begin(), tind.end()), tind.end());

			vector<Point3D> parr(pnarr.size());
			ModelSmooth::applyalgo(pnarr, parr, a, p...);
			ModelSmooth::update(pnarr, parr);

			return tind;
		}

		void addselect(const Point3D& p, double(&modelview)[16])
		{
			vt d = colorht->getcurrtd();
			if (ID2Tr.find(d) == ID2Tr.end())
			{
				
				ID2Tr.insert(make_pair(d, set<tri_ind_t>()));
				RegionProperty R(DefaultValue->Get(0), DefaultValue->Get(1), d);
				ID2Pr.insert(make_pair(d, RegionProperty::MakeProperty(DefaultValue->Get(0), DefaultValue->Get(1), d)));
				atid.insert(make_pair(d, make_pair(DefaultValue->Get(1), DefaultValue->Get(0))));
			}
			selecttri(p, modelview, [&](tri_ind_t i)
			{
				if (Tr2ID.find(i) != Tr2ID.end())
				{
					ID2Tr[Tr2ID[i]].erase(i);
				}
				Tr2ID[i] = d;
				ID2Tr[d].insert(i);
			});
		}

		void removeselect(const Point3D& p, double(&modelview)[16])
		{
			selecttri(p, modelview, [&](tri_ind_t i)
			{
				if (Tr2ID.find(i) != Tr2ID.end())
				{
					ID2Tr[Tr2ID[i]].erase(i);
					if (ID2Tr[Tr2ID[i]].empty())
					{
						ID2Tr.erase(Tr2ID[i]);
						ID2Pr.erase(Tr2ID[i]);
						atid.erase(Tr2ID[i]);
					}
					
				}
				Tr2ID.erase(i);
			});
		}

		void swapselect(const Point3D& p)
		{
			auto pp = mmodel->PointinModelCoor(p);
			auto ind = mmodel->searchTriangleind(pp);
			if (ind < mmodel->tarr.size() && Tr2ID.find(ind) != Tr2ID.end())
			{
				vt d = colorht->getcurrtd();
				auto v = Tr2ID[ind];
				if (v != d)
				{
					if (ID2Tr.find(d) == ID2Tr.end())
					{
						ID2Tr.insert(make_pair(d, move(ID2Tr[v])));
						RegionProperty R(DefaultValue->Get(0), DefaultValue->Get(1), d);
						ID2Pr.insert(make_pair(d, R));
						atid.insert(make_pair(d, make_pair(DefaultValue->Get(1), DefaultValue->Get(0))));
					}
					else
					{
						ID2Tr[d].insert(ID2Tr[v].begin(), ID2Tr[v].end());						
					}
					ID2Tr.erase(v);
					atid.erase(v);
					for (auto ii : ID2Tr[d])
					{
						Tr2ID[ii] = d;
					}

				}
			}
		}

		void reverseselect()
		{
			if (mmodel == nullptr)
			{
				return;
			}

			vt d = colorht->getcurrtd();
			std::vector<tri_ind_t> id;
			std::vector<tri_ind_t> id2;
			std::vector<tri_ind_t> id3;

			id.reserve(Tr2ID.size());
			for (auto&& tid : Tr2ID)
			{
				id.push_back(tid.first);
			}
			id2.reserve(mmodel->tarr.size());
			for (size_t i = 0; i < mmodel->tarr.size(); i++)
			{
				id2.push_back(i);
			}

			id3.reserve(id2.size() - id.size());
			std::set_difference(id2.begin(), id2.end(), id.begin(), id.end(), std::back_inserter(id3));

			Tr2ID.clear();
			ID2Tr.clear();
			ID2Pr.clear();
			atid.clear();

			ID2Tr.insert(make_pair(d, set<tri_ind_t>()));
			std::vector<std::pair<tri_ind_t, vt>> totcid;
			totcid.reserve(id3.size());

			for (auto&& tid : id3)
			{
				totcid.push_back({tid, d});
			}

			Tr2ID.insert(totcid.begin(), totcid.end());
			ID2Tr[d].insert(id3.begin(), id3.end());
			RegionProperty R(DefaultValue->Get(0), DefaultValue->Get(1), d);
			ID2Pr.insert(make_pair(d, R));
			atid.insert(make_pair(d, make_pair(DefaultValue->Get(1), DefaultValue->Get(0))));
		}

		map<Point3DNode*, double> get_mpnd(double shiftd, unsigned int level = 10) const;
		map<Point3DNode*, double> get_PointLoc(double shiftd, int mode, unsigned int level = 10) const;

		bool test(const shared_ptr<Model>& sm) const
		{
			return mmodel != sm || Tr2ID.size() == 0;
		}

		spvector<baseTriangle> extract_select() const
		{
			spvector<baseTriangle> result;
			for (auto&& t:Tr2ID)
			{
				result.push_back(mmodel->tarr[t.first]);
			}
			return result;
		}

		spvector<baseTriangle> filp_select() const
		{
			spvector<baseTriangle> result(mmodel->tarr.size());

			copy(mmodel->tarr.begin(), mmodel->tarr.end(), result.begin());
			for (auto&& ind : Tr2ID)
			{
				auto& t = result[ind.first];
				result[ind.first] = make_shared<myMath::Triangle>(t->at(0), t->at(2), t->at(1));
			}
			return result;
		}
	};

	class ColoredModel :public VBOobj
	{
	public:
		shared_ptr<Model> mmodel;
		shared_ptr<Model> mmodel2;
		shared_ptr<ColorHistogram2> colorht;
		
		ColoredModel(const shared_ptr<Model>& ma, const shared_ptr<Model>& mb, const shared_ptr<ColorHistogram2>& ct);
		virtual ~ColoredModel(){};
	};

	class ColoredModel2 :public VBOobj
	{
	public:
		shared_ptr<Model> mmodel;
		shared_ptr<ColorHistogram2> colorht;

		ColoredModel2(const shared_ptr<Model>& m, const shared_ptr<ColorHistogram2>& ct);
		virtual ~ColoredModel2() {};
	};

	class RegionGrowing
	{
	private:
		typedef size_t groupind;
		struct seedtemp1
		{
			tri_ind_t tind;
			groupind group;
			Point3D cen;			
			double dis;
		};
		struct seedtemp2
		{
			double cur;
			groupind group;
			tri_ind_t tind;

			bool operator< (const seedtemp2& s2) const
			{
				return cur < s2.cur;
			}
			bool operator> (const seedtemp2& s2) const
			{
				return cur > s2.cur;
			}
		};
	public:
		RegionGrowing(sp<Model>& m) :mmodel(m), groupcount(0), growstep(0), growcheck(nullptr)
		{

		}
		~RegionGrowing()
		{

		}

		sp<Model> mmodel;
		groupind groupcount;
		int growstep;
		bool* growcheck;

		vector<pair<groupind, tri_ind_t>> seed0;
		vectors<tri_ind_t> submodelind;
		vectors<gl_tri_ind_t> renderind;

		Point3D tricen(tri_ind_t tind) const;
		double tricur(tri_ind_t tind1, tri_ind_t tind2, vertex_ind_t pind) const;
		void addseed(const Point3D& p, bool newseed = false);
		void Growing1();
		void Growing();

		spvector<Model> resultmodel();
	};

	class RegionGrowingControl
	{
	public:
		sp<Model> mmodel;
		typedef vectors<Point3D> seed_cont;
		sp<RegionGrowing> rg;
		sp<OBJrender::test2::layerbase<world_layer>> rgr;
		seed_cont seed;
		bool showseed;

		RegionGrowingControl(sp<Model>& m) :mmodel(m), showseed(true)
		{

		}

		bool iscurrentmodel(sp<Model>& m)
		{
			return m == mmodel;
		}
		void setmodel(sp<Model>& m);
		
		void addseed(const Point3D& p, bool newflag);
		
		void removeallseed()
		{
			seed.clear();
			update();
		}
		void removelastseed()
		{
			if (!seed.empty())
			{
				if (!seed.back().empty())
				{
					seed.back().pop_back();
				}

				if (seed.back().empty())
				{
					seed.pop_back();
				}
				update();
			}
		}
		void update();
		spvector<Model> resultmodel()
		{
			if (rg!=nullptr)
			{
				return rg->resultmodel();
			}
			else
			{
				return spvector<Model>();
			}
		}
	};

	class VaildRegionGrowing
	{
	private:
		struct eihv
		{
			size_t hv;
			float lensq;
			shared_ptr<baseTriangle> t;
			tri_ind_t tid;
			edge_ind_t eid;
			bool b;
			bool operator<(const eihv& e2) const
			{
				return hv < e2.hv || (!(e2.hv < hv) && lensq < e2.lensq);
			}
		};

		int* growcheck;
		hash<Point3D> hp;
		vector<eihv> earr;
		vector<eihv> estack;
		vector<tri_ind_t> growstack;
		void addsubtri(sp<baseTriangle>& t);
		void subgrow(triid_pid t, tri_ind_t currseed);
	public:
		VaildRegionGrowing(sp<Model>& m) :mmodel(m), growstep(0), growcheck(nullptr), seedtri(0)
		{
		}
		sp<Model> mmodel;
		Model::TriCrossResult tcr;
		Model::TriSeparateResult tsr;

		int growstep;
		vector<tri_ind_t> seedtri;
		vector<gl_tri_ind_t> renderind;
		vector<tri_ind_t> resultind;
		spvector<baseTriangle> subtri;

		void setseedtri(const Point3D& p);
		void Growing();
		sp<Model> resultmodel();
	};

	class ModelDigHole
	{
	public:
		struct cydata
		{
			Cylinder cy;
			Point3D pu;
			Point3D pr;
			//double theta;
			double k;

			TransformationMatrix tm;

			cydata()
			{
			}
			cydata(const Cylinder& cy_) :cy(cy_)//, theta(0)
			{
				update();
			}

			void update_pu()
			{
				pu = cy.cen + cy.n * cy.height;
			}

			void update_pr()
			{
				uVector3D rn = Vector3D::genVVec(cy.n);
				pr = cy.cen + rn * cy.radius;
			}

			void update()
			{
				update_pu();
				update_pr();
			}

			bool isinside(const Point3D& p) const
			{
				return cy.isinside(p);
			}
		};
		vector<cydata> cyarr;
		cydata cytemp;
		bool hastemp;

		sp<ModelBrush> brush;

		ModelDigHole() :hastemp(false)
		{
			
		}

		void addtempCylinder()
		{
			cyarr.push_back(cytemp);
			hastemp = false;
		}
		void addCylinder(const Cylinder& cy_);
		
		sp<Model> Triangulation(bool twoside) const;
		sp<Model> Triangulation2(bool twoside) const;
		sp<Model> Triangulation3() const;
		sp<Model> Triangulation4() const;
		sp<Model> Triangulation5() const;
	};

	class ModelPlane :public Plane
	{
	public:
		Point3D upcen;
		Point3D cor[4];

		ModelPlane(const Vector3D& N, const Point3D& C) :Plane(N, C)
		{
			upcen = cen + n*10;


			Vector3D v = {0, 0, 1};
			uVector3D rr = Vector3D::cross(v, n);

			double ctheta = Vector3D::dot(v, n);
			double stheta = sqrt(1 - ctheta*ctheta);
			double vtheta = 1 - ctheta;
			Vector3D rmat[3] =
			{
				{rr[0] * rr[0] * vtheta + ctheta, rr[0] * rr[1] * vtheta - rr[2] * stheta, rr[0] * rr[2] * vtheta + rr[1] * stheta},
				{rr[0] * rr[1] * vtheta + rr[2] * stheta, rr[1] * rr[1] * vtheta + ctheta, rr[1] * rr[2] * vtheta - rr[0] * stheta},
				{rr[0] * rr[2] * vtheta - rr[1] * stheta, rr[1] * rr[2] * vtheta + rr[0] * stheta, rr[2] * rr[2] * vtheta + ctheta}
			};

			auto rxp = [&](const Point3D& p)->Point3D
			{
				double x = Vector3D::dot(rmat[0], p) + cen.x;
				double y = Vector3D::dot(rmat[1], p) + cen.y;
				double z = Vector3D::dot(rmat[2], p) + cen.z;

				return {x, y, z};
			};

			cor[0] = rxp({-50, -50, 0});
			cor[1] = rxp({50, -50, 0});
			cor[2] = rxp({50, 50, 0});
			cor[3] = rxp({-50, 50, 0});
		}
	};

	class ModelPlane2
	{
	public:
		myMath::TransformationMatrix tm;
		myMath::Point3D cen;
		myMath::Vector3D n;
		myMath::Point3D p[4];
	};

	class ModelEdge
	{	
		struct ctt
		{
			double cur;
			tri_ind_t tind;
			triid_pid adjt;
			bool operator< (const ctt& c) const
			{
				return cur < c.cur;
			}
		};
	public:
		shared_ptr<Model> mmodel;
		vector<ctt> edgearr;
		vector<ctt>::const_iterator eit;
		double lim;

		ModelEdge(shared_ptr<Model> M) :mmodel(M)
		{
			M->showmesh = false;
			for (tri_ind_t i = 0; i < M->tarr.size(); i++)
			{
				for (auto&& adj2t : M->tarrinfo[i]->b)
				{
					if (i < adj2t.tid())
					{
						auto cur = triang(i, adj2t.tid());
						
						edgearr.push_back({cur, i, adj2t});
					}
				}
			}
			sort(edgearr.begin(), edgearr.end());

			setlim(5);
		}

		double triang(tri_ind_t tind1, tri_ind_t tind2) const
		{
			auto& n1 = mmodel->tarr.at(tind1)->n;
			auto& n2 = mmodel->tarr.at(tind2)->n;
			auto nn = Vector3D::dot(n1, n2);
			nn = nn > 0 ? std::min(nn, 1.0) : std::max(nn, -1.0);			
			return std::acos(nn);
		}

		spvectors<Segment3D> EdgeSeparate() const
		{
			vector<Model::Edgeind> cttarrit;
			cttarrit.reserve(std::distance(eit, edgearr.end()));
			for (auto it = eit; it != edgearr.end(); it++)
			{
				auto pind = it->adjt.pid();
				cttarrit.push_back({it->tind, (pind + 1) % 3u, (pind + 2) % 3u});
			}

			auto cttarr = mmodel->searchEdgecontour(cttarrit);
			
			auto se2 = spvectors<Segment3D>();
			se2.reserve(cttarr.size());
			for (auto&& i : cttarr)
			{
				spvector<Segment3D> tempss;
				tempss.reserve(i.size());
				for (auto&& j : i)
				{
					auto& t = mmodel->tarr[j.a];
					tempss.push_back(make_shared<Segment3D>(t->at(j.b), t->at(j.c)));
				}
				se2.push_back(tempss);
			}
			return se2;
		}

		spvectors<QSegment3D> EdgeSeparateQ() const
		{
			vector<Model::Edgeind> cttarrit;
			cttarrit.reserve(std::distance(eit, edgearr.end()));
			for (auto it = eit; it != edgearr.end(); it++)
			{
				auto pind = it->adjt.pid();
				cttarrit.push_back({it->tind, (pind + 1) % 3u, (pind + 2) % 3u});
			}

			auto cttarr = mmodel->searchEdgecontour(cttarrit);

			auto se2 = spvectors<QSegment3D>();
			se2.reserve(cttarr.size());
			for (auto&& i : cttarr)
			{
				spvector<QSegment3D> tempss;
				tempss.reserve(i.size());
				for (auto&& j : i)
				{
					auto& t = mmodel->tarr[j.a];
					auto& adjt = mmodel->tarrinfo[j.a]->b;
					auto it = find_if(adjt.begin(), adjt.end(), [&j](const triid_pid& aa){return aa.pid() == (3 - j.b - j.c); });
					tempss.push_back(make_shared<QSegment3D>(t->at(j.b), t->at(j.c), j.a, it->tid()));
				}
				se2.push_back(tempss);
			}
			return se2;
		}

		spvectors<baseTriangle> EdgetoTube(double r) const
		{
			auto edgess = EdgeSeparate();

			auto checkstart = [](Segment3D& a, Segment3D& b)
			{
				if (a.pstart == b.pstart || a.pstart == b.pend)
				{
					a.reversal();
				}

				if (!(a.pend == b.pstart))
				{
					b.reversal();
				}
			};
			const size_t nedge = 41;
			auto gencircle = [r, nedge](uVector3D n, Point3D cen)->vector<Point3D>
			{
				const double delta = 2 * M_PI / nedge;
				Vector3D v = {0, 0, 1};
				uVector3D rr = Vector3D::cross(v, n);

				double ctheta = Vector3D::dot(v, n);
				double stheta = sqrt(1 - ctheta*ctheta);
				double vtheta = 1 - ctheta;
				Vector3D rmat[3] =
				{
					{rr[0] * rr[0] * vtheta + ctheta, rr[0] * rr[1] * vtheta - rr[2] * stheta, rr[0] * rr[2] * vtheta + rr[1] * stheta},
					{rr[0] * rr[1] * vtheta + rr[2] * stheta, rr[1] * rr[1] * vtheta + ctheta, rr[1] * rr[2] * vtheta - rr[0] * stheta},
					{rr[0] * rr[2] * vtheta - rr[1] * stheta, rr[1] * rr[2] * vtheta + rr[0] * stheta, rr[2] * rr[2] * vtheta + ctheta}
				};

				auto rxp = [&](const Point3D& p)->Point3D
				{
					double x = Vector3D::dot(rmat[0], p) + cen.x;
					double y = Vector3D::dot(rmat[1], p) + cen.y;
					double z = Vector3D::dot(rmat[2], p) + cen.z;

					return {x, y, z};
				};

				std::vector<Point3D> basepoint(nedge);
				for (size_t i = 0; i < nedge; i++)
				{
					double deg = delta*i;
					double x = r*cos(deg);
					double y = r*sin(deg);
					basepoint[i] = rxp({x, y, 0});					
				}
				return basepoint;
			};

			spvectors<baseTriangle> resulttriarr;
			for (auto&& edges : edgess)
			{
				if (edges.empty())
				{
					continue;
				}
				spvector<baseTriangle> resulttri;
				//////////////////////////////////////////////////////////////////////////
				for (size_t i = 0; i < edges.size()-1; i++)
				{
					checkstart(*edges[i], *edges[i + 1]);
				}
				bool isloop = edges.front()->pstart == edges.back()->pend;

				vector<Vector3D> narr;
				vectors<Point3D> parr;
				
				narr.push_back(isloop ? edges.front()->n + edges.back()->n : edges.front()->n);
				for (size_t i = 1; i < edges.size(); i++)
				{
					narr.push_back(edges[i - 1]->n + edges[i]->n);
				}

				for (size_t i = 0; i < edges.size(); i++)
				{
					parr.push_back(gencircle(narr[i], edges[i]->pstart));
				}
				if (isloop)
				{
					parr.push_back(parr.front());
				}
				else
				{
					parr.push_back(gencircle(edges.back()->n, edges.back()->pend));
				}

				for (size_t i = 0; i < edges.size(); i++)
				{
					for (size_t i1 = 0; i1 < nedge; i1++)
					{
						size_t i2 = (i1 + 1) % nedge;
						resulttri.push_back(make_shared<Triangle>(parr[i][i1], parr[i][i2], parr[i + 1][i2]));
						resulttri.push_back(make_shared<Triangle>(parr[i][i1], parr[i + 1][i2], parr[i + 1][i1]));
					}
				}
				if (!isloop)
				{
					for (size_t i1 = 0; i1 < nedge; i1++)
					{
						size_t i2 = (i1 + 1) % nedge;
						resulttri.push_back(make_shared<Triangle>(edges.front()->pstart, parr.front()[i2], parr.front()[i1]));
						resulttri.push_back(make_shared<Triangle>(edges.back()->pend, parr.back()[i1], parr.back()[i2]));
					}
				}

				//////////////////////////////////////////////////////////////////////////
				resulttriarr.push_back(move(resulttri));
// 				if (resulttriarr.size() == 2)
// 				{
// 					break;
// 				}
			}

			return resulttriarr;
		}

		void setlim(double ddeg)
		{
			lim = ddeg;
			if (lim < 0)
			{
				lim = 0;
			}
			else if (lim > 180)
			{
				lim = 180;
			}

			//eit = std::upper_bound(edgearr.begin(), edgearr.end(), lim*M_PI / 180.0, [](double c1, const ctt& c2){return c1 < c2.cur; });
			eit = std::lower_bound(edgearr.begin(), edgearr.end(), lim*M_PI / 180.0, [](const ctt& c1, double c2){return c1.cur < c2; });
		}
		void addlim(double ddeg)
		{
			setlim(lim + ddeg);
		}
		void sublim(double ddeg)
		{
			setlim(lim - ddeg);
		}
	};

	class ModelHoleFill
	{
	public:
		static spvector<baseTriangle> getresult(const spvector<Segment3D>& ss);
		static spvector<baseTriangle> getresult(const sp<Model> mmodel, const vector<Model::Edgeind>& ss);
		static spvector<baseTriangle> getresult(const sp<Model> mmodel, const vectors<Model::Edgeind>& sss);

		static spvector<baseTriangle> getresult(const vector<Point3D>& pp);
	};
	
	class ModelLaserCutter
	{
	public:
		double laserlen;
		double safed;
		uVector3D laserdir;
		vectors<Segment3D> sarrs;
		bool laserenable;
		Point3D lasercen;
		Point3D s;
		Point3D e;

		TransformationMatrix tm;

		ModelLaserCutter()
		{
			laserlen = 10.0;
			safed = 5.0;
			laserdir = Vector3D {0,0,-1};
			//laserdir = Vector3D {-1, 1, -1};
			//laserdir = Vector3D {-1, 1, -sqrt(2.0)};
			laserenable = false;
		}

		sp<ModelSplitPlane> generateplane()
		{
			if (sarrs.empty())
			{
				return nullptr;
			}

			//auto s = make_shared<Splitplane>();
			spvector<baseTriangle> tarr;
			for (auto&& sarr : sarrs)
			{
				if (sarr.size() < 2)
				{
					continue;
				}

				for (size_t i = 0; i < sarr.size() - 1; i++)
				{
					tarr.push_back(make_shared<Triangle>(sarr[i].pstart, sarr[i + 1].pstart, sarr[i].pend));
					tarr.push_back(make_shared<Triangle>(sarr[i + 1].pstart, sarr[i + 1].pend, sarr[i].pend));
// 					s->splane.push_back(std::make_shared<myMath::Triangle>(sarr[i].pstart, sarr[i + 1].pstart, sarr[i].pend));
// 					s->splane.push_back(std::make_shared<myMath::Triangle>(sarr[i + 1].pstart, sarr[i + 1].pend, sarr[i].pend));
// 					s->psarr.push_back(sarr[i].pstart);
// 					s->pearr.push_back(sarr[i].pend);
				}
			}
			return make_shared<ModelSplitPlane>(make_shared<Model>(tarr));
			//return s;
		}
	};
}

namespace std
{
	template <>
	struct hash < myModel::ModelTriangle >
	{
		std::size_t operator()(const myModel::ModelTriangle& t) const
		{
			using std::hash;
			hash<myMath::Point3D> hp;
			return hash<size_t>()(hp(*t.a) ^ hp(*t.b) ^ hp(*t.c));
		}
	};
}

namespace OBJrender
{
// 	typedef template_Render<std::string> string_Render;
// 	typedef template_Render<myMath::Point3D> Point3D_Render;
// 	typedef template_Render<myMath::baseTriangle> Triangle_Render;
// 	typedef template_Render<myMath::Segment3D> Segment3D_Render;
// 	typedef template_Render<spvector<myMath::Segment3D>> Segment3Ds_Render;
// 	typedef template_Render<spvectors<myMath::Segment3D>> Segment3Dss_Render;
// 	typedef template_Render<std::pair<sp<myModel::Model>, spvectors<myMath::QSegment3D>>> QSegment3Dss_Render;
// 	typedef template_Render<spvector<myMath::baseTriangle>> Triangles_Render;
// 	typedef template_Render<spvectors<myMath::baseTriangle>> Triangless_Render;
// 	typedef template_Render<myModel::BoundBox> BoundBox_Render;
// 	typedef template_Render<std::vector<myModel::BoundBox*>> BoundBoxs_Render;
// 	typedef template_Render<myModel::OctCloud> OctCloud_Render;
// 	typedef template_Render<myModel::AABBtree> AABBtree_Render;
// 	typedef template_Render<myModel::Model> Model_Render;
// 	typedef template_Render<myModel::ModelManger> Manger_Render;
// 	typedef template_Render<myModel::ColoredModel> CModel_Render;
// 	typedef template_Render<myModel::ModelSelect> Select_Render;
// 	typedef template_Render<myModel::ModelPlane> Plane_Render;
// 	
// 
// 	typedef template_Render<myModel::ModelDigHole> Hole_Render;
// 	typedef template_Render<myModel::RegionGrowing> Growing_Render;
// 	typedef template_Render<myModel::VaildRegionGrowing> VGrowing_Render;
// 	typedef uiRender<myModel::ModelCutter> Cutter_Render;
// 	typedef uiRender<myModel::ColorHistogram2> CH2_Render;
// 	typedef uiRender<myModel::ColorHistogram> CH_Render;
// 	typedef uiRender<myModel::ModelSelect2> Select2_Render;
// 	typedef uiRender<myModel::ModelEdge> Edge_Render;

//////////////////////////////////////////////////////////////////////////
	
// 	typedef template_Render<std::string> string_Render;
// 	typedef template_Render<myMath::Point3D> Point3D_Render;
// 	typedef test::test_render<myMath::baseTriangle, wolrd_layer> Triangle_Render;
// 	typedef test::test_render<myMath::Segment3D, wolrd_layer> Segment3D_Render;
// 	typedef test::test_render<spvector<myMath::Segment3D>, wolrd_layer> Segment3Ds_Render;
// 	typedef test::test_render<spvectors<myMath::Segment3D>, wolrd_layer> Segment3Dss_Render;
// 	typedef test::test_render<std::pair<sp<myModel::Model>, spvectors<myMath::QSegment3D>>, wolrd_layer> QSegment3Dss_Render;
// 	typedef test::test_render<spvector<myMath::baseTriangle>, wolrd_layer> Triangles_Render;
// 	typedef test::test_render<spvectors<myMath::baseTriangle>, wolrd_layer> Triangless_Render;
// 	typedef test::test_render<myModel::BoundBox, wolrd_layer> BoundBox_Render;
// 	typedef test::test_render<std::vector<myModel::BoundBox*>, wolrd_layer> BoundBoxs_Render;
// 	typedef test::test_render<myModel::OctCloud, wolrd_layer> OctCloud_Render;
// 	typedef test::test_render<myModel::AABBtree, wolrd_layer> AABBtree_Render;
// 	typedef test::test_render<myModel::Model, wolrd_layer> Model_Render;
// 	typedef test::test_render<myModel::ModelManger, wolrd_layer> Manger_Render;
// 	typedef test::test_render<myModel::ColoredModel, wolrd_layer> CModel_Render;
// 	typedef test::test_render<myModel::ModelSelect, wolrd_layer> Select_Render;
// 	typedef test::test_render<myModel::ModelPlane, wolrd_layer> Plane_Render;
// 
// 
// 	typedef test::test_render<myModel::ModelDigHole, wolrd_layer> Hole_Render;
// 	typedef test::test_render<myModel::RegionGrowing, wolrd_layer> Growing_Render;
// 	typedef test::test_render<myModel::VaildRegionGrowing, wolrd_layer> VGrowing_Render;
// 	typedef test::test_render<myModel::ModelCutter, wolrd_layer, ui_layer> Cutter_Render;
// 	typedef test::test_render<myModel::ColorHistogram2, ui_layer> CH2_Render;
// 	typedef test::test_render<myModel::ColorHistogram, wolrd_layer> CH_Render;
// 	typedef test::test_render<myModel::ModelSelect2, wolrd_layer, ui_layer> Select2_Render;
// 	typedef test::test_render<myModel::ModelEdge, wolrd_layer, ui_layer> Edge_Render;

//////////////////////////////////////////////////////////////////////////

	typedef test2::layerrender<std::string> string_Render;
	typedef test2::layerrender<myMath::Point3D> Point3D_Render;
	typedef test2::layerrender<std::vector<myMath::Point3D>, world_layer> Point3Ds_Render;
	typedef test2::layerrender<myMath::baseTriangle, world_layer> Triangle_Render;
	typedef test2::layerrender<myMath::Segment3D, world_layer> Segment3D_Render;
	typedef test2::layerrender<spvector<myMath::Segment3D>, world_layer> Segment3Ds_Render;
	typedef test2::layerrender<spvectors<myMath::Segment3D>, world_layer> Segment3Dss_Render;
	typedef test2::layerrender<std::pair<sp<myModel::Model>, spvectors<myMath::QSegment3D>>, world_layer> QSegment3Dss_Render;
	typedef test2::layerrender<spvector<myMath::baseTriangle>, world_layer> Triangles_Render;
	typedef test2::layerrender<spvectors<myMath::baseTriangle>, world_layer> Triangless_Render;
	typedef test2::layerrender<myModel::BoundBox, world_layer> BoundBox_Render;
	typedef test2::layerrender<std::vector<myModel::BoundBox*>, world_layer> BoundBoxs_Render;
	typedef test2::layerrender<myModel::OctCloud, world_layer> OctCloud_Render;
	typedef test2::layerrender<myModel::AABBtree, world_layer> AABBtree_Render;
	typedef test2::layerrender<myModel::ModelPCA, world_layer> PCAAxis_Render;
	//typedef test2::layerrender<myModel::LDI, world_layer> LDI_Render;
	typedef test2::layerrender<myModel::Model, world_layer> Model_Render;
	//typedef test2::layerrender<myModel::ModelManger, world_layer> Manger_Render;
	typedef test2::multilayerrender<myModel::ModelManger, world_layer, world_transparent_layer> Manger_Render;
	typedef test2::layerrender<myModel::ColoredModel, world_layer> CModel_Render;
	typedef test2::layerrender<myModel::ColoredModel2, world_layer> Curvature_Render;
	typedef test2::layerrender<myModel::ModelSelect, world_layer> Select_Render;
	typedef test2::layerrender<myModel::ModelPlane, world_layer> Plane_Render;
	typedef test2::layerrender<myModel::ModelPlane2, world_transparent_layer> Plane2_Render;	

	typedef test2::layerrender<myModel::ModelSplitPlane, world_transparent_layer> SPlane_Render;
	typedef test2::layerrender<myModel::ModelDigHole, world_layer> Hole_Render;
	typedef test2::layerrender<myModel::RegionGrowing, world_layer> Growing_Render;
	typedef test2::layerrender<myModel::RegionGrowingControl, world_layer> GrowingC_Render;
	
	typedef test2::layerrender<myModel::VaildRegionGrowing, world_layer> VGrowing_Render;
	typedef test2::layerrender<myModel::ModelCutter, ui_layer> Cutter_Render;
	typedef test2::layerrender<myModel::ColorHistogram2, ui_layer> CH2_Render;
	typedef test2::multilayerrender<myModel::ColorHistogram, world_layer> CH_Render;
	typedef test2::multilayerrender<myModel::ModelSelect2, world_layer, ui_layer> Select2_Render;
	typedef test2::multilayerrender<myModel::ModelSelect3, world_layer, ui_layer> Select3_Render;
	typedef test2::multilayerrender<myModel::ModelEdge, world_layer, ui_layer> Edge_Render;
	typedef test2::multilayerrender<myModel::ModelLaserCutter, world_layer, world_transparent_layer> LaserCutter_Render;
}

#endif