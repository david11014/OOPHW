#include "stdafx.h"

#include <regex>
#include <windows.h>

#include "STLmodel.h"

//#include "Delaunay.h"
//#include "cgal.h"
//#include "alglib/dataanalysis.h"

#ifdef useParallel
#include <ppl.h>
#include <concurrent_vector.h>
#include <concurrent_unordered_map.h>
using namespace concurrency;
#endif

using namespace std;
using namespace myMath;

//////////////////////////////////////////////////////////////////////////
void myModel::LoadStlBinary(const string& fileName, tri_ind_t& trianglecount, float*& dest)
{
	using namespace std;
	mtime3.clear();
	mtime3.start();
	try
	{
		ifstream fin(fileName, ios::in | ios::binary);

		if (!fin.good())
		{
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
		}

		fin.seekg(0, ios::end);
		string contents((size_t)fin.tellg(), '0');
		fin.seekg(0, ios::beg);
		fin.read(&contents[0], contents.size());
		fin.close();

		memcpy(&trianglecount, &contents[80], 4);

		char solid[6] = {'\0'};
		memcpy(solid, &contents[0], 5);
		if (strcmp(solid, "solid") == 0 && contents.size() != trianglecount * 50 + 84)
		{
			throw new domain_error("asci file.");
		}
		if (trianglecount < 0 || trianglecount > 1e8)
		{
			throw new out_of_range("triangle size incorrect.");
		}

		if (dest)
		{
			delete[] dest;
		}
		dest = new float[trianglecount * 3 * 3];

		float* tempp = dest;
		char* tempb = &contents[84] + 12;


		auto check_hash = [](float* f)
		{
			auto hashfunc = [](float* f)
			{
				return _Hash_seq((unsigned char*)f, 12);
			};
			auto hh1 = hashfunc(f);
			auto hh2 = hashfunc(f + 3);
			auto hh3 = hashfunc(f + 6);
			return hh1 == hh2 || hh2 == hh3 || hh3 == hh1;
		};
		auto check_inf_nan = [](float* f)
		{
			for (unsigned int j = 0; j < 9; j++)
			{
				if (isinf(*(f + j)) || isnan(*(f + j)))
				{
					return true;
				}
			}
			return false;
		};

		for (tri_ind_t i = 0; i < trianglecount; i++, tempb += 50)
		{
			memcpy(tempp, tempb, 36);

			if (check_inf_nan(tempp) || check_hash(tempp))
			{
				continue;
			}

			tempp += 9;
		}
		trianglecount = (tempp - dest) / 9;
	}
	catch (exception*)
	{
		trianglecount = 0;
		if (dest)
		{
			delete[] dest;
			dest = nullptr;
		}
	}
	mtime3.end();
	mtime3.print("Loadstl: ");
	return;
}

float* myModel::LoadStlBinary(const string& fileName, tri_ind_t& trianglecount)
{
	float* dest = nullptr;
	LoadStlBinary(fileName, trianglecount, dest);
	return dest;
}

void myModel::SaveStlBinary(const string& fileName, const spvector<baseTriangle>& V)
{
	using namespace std;
	try
	{
		ofstream fout(fileName.c_str(), ios::out | ios::binary);

		if (!fout.good())
		{
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
		}

		const unsigned int tmpLen = 80;
		char tmpBuff[tmpLen] = {'\0'};

		tri_ind_t numofTri = V.size();
		fout.write(tmpBuff, tmpLen);
		fout.write((char *)&numofTri, 4);

		char tempchar[2] = {'\0', '\0'};
		float tempf;
		char* tempc = (char*)&tempf;
		for (auto&& i : V)
		{
			tempf = (float)i->n.x; fout.write(tempc, 4);
			tempf = (float)i->n.y; fout.write(tempc, 4);
			tempf = (float)i->n.z; fout.write(tempc, 4);

			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					tempf = (float)i->at(j).p[k];
					fout.write(tempc, 4);
				}
			}
			fout.write(tempchar, 2);
		}

		fout.clear();
		fout.close();
	}
	catch (exception*)
	{
		return;
	}
}
//////////////////////////////////////////////////////////////////////////
//myModel::SubTriangles::SubTriangles(const ModelTriangle* T, const spvector<QSegment3D>& S)
//{
//	vector<tuple<double, double, double>> dd; dd.reserve(3);
//	vector<tuple<double, double, double>> ss; ss.reserve(S.size() * 2);
//
//	dd.emplace_back(T->a->x, T->a->y, T->a->z);
//	dd.emplace_back(T->b->x, T->b->y, T->b->z);
//	dd.emplace_back(T->c->x, T->c->y, T->c->z);
//
//	for (auto&& qs : S)
//	{
//		ss.emplace_back(qs->pstart.x, qs->pstart.y, qs->pstart.z);
//		ss.emplace_back(qs->pend.x, qs->pend.y, qs->pend.z);
//	}
//
//	dtri_ind = testtest(dd, ss);
//}
//
//myModel::SubTriangles::SubTriangles(const ModelTriangle* T, const spvector<ModelTriangle>& S)
//{
//	vector<tuple<double, double, double>> dd; dd.reserve(3);
//	vector<tuple<double, double, double>> ss; ss.reserve(S.size() * 3);
//	dd.emplace_back(T->a->x, T->a->y, T->a->z);
//	dd.emplace_back(T->b->x, T->b->y, T->b->z);
//	dd.emplace_back(T->c->x, T->c->y, T->c->z);
//	for (size_t i = 0; i < S.size(); i++)
//	{
//		ss.emplace_back(S[i]->a->x, S[i]->a->y, S[i]->a->z);
//		ss.emplace_back(S[i]->b->x, S[i]->b->y, S[i]->b->z);
//		ss.emplace_back(S[i]->c->x, S[i]->c->y, S[i]->c->z);
//	}
//	//cout << ss.size() << endl;
//	dtri_ind = testtest2(dd, ss);
//}
//
//myModel::SubTriangles::SubTriangles(const ModelTriangle* T, const double edgelength)
//{
//	vector<tuple<double, double, double>> dd;
//
//	dd.emplace_back(T->a->x, T->a->y, T->a->z);
//	dd.emplace_back(T->b->x, T->b->y, T->b->z);
//	dd.emplace_back(T->c->x, T->c->y, T->c->z);
//
//	if (T->lab() > edgelength)
//	{
//		size_t lcount = lround(T->lab() / edgelength);
//
//		Vector3D ba = (*T->b - *T->a)/double(lcount);
//		for (size_t i = 1; i < lcount; i++)
//		{
//			auto pp = *T->a + ba*double(i);
//			dd.emplace_back(pp.x, pp.y, pp.z);
//		}
//	}
//
//	if (T->lbc() > edgelength)
//	{
//		size_t lcount = lround(T->lbc() / edgelength);
//
//		Vector3D cb = (*T->c - *T->b) / double(lcount);
//		for (size_t i = 1; i < lcount; i++)
//		{
//			auto pp = *T->b + cb*double(i);
//			dd.emplace_back(pp.x, pp.y, pp.z);
//		}
//	}
//
//	if (T->lca() > edgelength)
//	{
//		size_t lcount = lround(T->lca() / edgelength);
//
//		Vector3D ac = (*T->a - *T->c) / double(lcount);
//		for (size_t i = 1; i < lcount; i++)
//		{
//			auto pp = *T->c + ac*double(i);
//			dd.emplace_back(pp.x, pp.y, pp.z);
//		}
//	}
//
//	dtri_ind = subedge(dd);
//}

spvector<baseTriangle> myModel::SubTriangles::Triangulation() const
{
	spvector<baseTriangle> triarr;

	if (dtri_ind.size() < 18)
	{
		return triarr;
	}

	for (size_t i = 0; i < dtri_ind.size(); i += 9)
	{
		Point3D aa = {dtri_ind[i], dtri_ind[i + 1], dtri_ind[i + 2]};
		Point3D bb = {dtri_ind[i + 3], dtri_ind[i + 4], dtri_ind[i + 5]};
		Point3D cc = {dtri_ind[i + 6], dtri_ind[i + 7], dtri_ind[i + 8]};
		Point3D::erasezero(aa);
		Point3D::erasezero(bb);
		Point3D::erasezero(cc);

		hash<Point3D> hp;
		size_t hh[3] = {hp(aa), hp(bb), hp(cc)};

		if (hh[0] == hh[1] || hh[0] == hh[2] || hh[1] == hh[2])
		{
			cout << "0";
			continue;
		}

		auto temptri = make_shared<Triangle>(aa, bb, cc);

		if (temptri->area() != 0 && !temptri->isslim())
		{
			triarr.push_back(move(temptri));
		}

	}

	return triarr;
}
//////////////////////////////////////////////////////////////////////////
//myModel::ModelPCA::ModelPCA(spvector<ModelTriangle> *tarr, BoundBox *bbox)
//{
//	CalModelPCA(tarr, bbox);
//}
//
//void myModel::ModelPCA::CalModelPCA(spvector<ModelTriangle> *tarr, BoundBox *bbox)
//{
//	double *P = new double[tarr->size() * 9];
//
//	for (tri_ind_t i = 0; i < tarr->size(); i++)
//	{
//		P[i * 9] = tarr->at(i)->a->x;
//		P[i * 9 + 1] = tarr->at(i)->a->y;
//		P[i * 9 + 2] = tarr->at(i)->a->z;
//
//		P[i * 9 + 3] = tarr->at(i)->b->x;
//		P[i * 9 + 4] = tarr->at(i)->b->y;
//		P[i * 9 + 5] = tarr->at(i)->b->z;
//
//		P[i * 9 + 6] = tarr->at(i)->c->x;
//		P[i * 9 + 7] = tarr->at(i)->c->y;
//		P[i * 9 + 8] = tarr->at(i)->c->z;
//	}
//
//	alglib::real_2d_array ptInput;
//	//ptInput.setcontent(sourceP.size(), 3,sourceP.data()->p);
//	ptInput.setcontent(tarr->size() * 3, 3, P);
//
//	// this guy gets passed in and is filled in with an integer status code
//	alglib::ae_int_t info;
//
//	// scalar values that describe variances along each eigenvector
//	alglib::real_1d_array eigValues;
//
//	// unit eigenvectors which are the orthogonal basis that we want
//	alglib::real_2d_array eigVectors;
//
//	// perform the analysis
//	//pcabuildbasis(ptInput, sourceP.size(), 3, info, eigValues, eigVectors);
//	pcabuildbasis(ptInput, tarr->size() * 3, 3, info, eigValues, eigVectors);
//
//
//
//	basis0.x = eigVectors[0][0];
//	basis0.y = eigVectors[1][0];
//	basis0.z = eigVectors[2][0];
//
//	basis1.x = eigVectors[0][1];
//	basis1.y = eigVectors[1][1];
//	basis1.z = eigVectors[2][1];
//
//	basis2.x = eigVectors[0][2];
//	basis2.y = eigVectors[1][2];
//	basis2.z = eigVectors[2][2];
//
//	this->eigValues[0] = eigValues[0];
//	this->eigValues[1] = eigValues[1];
//	this->eigValues[2] = eigValues[2];
//
//	O.x = (bbox->M.x + bbox->m.x) / 2;
//	O.y = (bbox->M.y + bbox->m.y) / 2;
//	O.z = (bbox->M.z + bbox->m.z) / 2;
//}
////////////////////////////////////////////////////////////////////////////
//
//void myModel::LDI::UpdateLDI()
//{
//	cout << "D:" << D << endl;
//
//	//設定兩平面軸 i j 與射線軸k
//	const int i[] = { 1,2,0 };
//	const int j[] = { 2,0,1 };
//	const int k[] = { 0,1,2 };
//
//	Point3D M, m;
//	M = bbbox->bboxcorner(0);
//	BoxMax[0] = (int(M.x / D) + 1);
//	BoxMax[1] = (int(M.y / D) + 1);
//	BoxMax[2] = (int(M.z / D) + 1);
//
//	m = bbbox->bboxcorner(7);
//	BoxMin[0] = (int(m.x / D) - 1);
//	BoxMin[1] = (int(m.y / D) - 1);
//	BoxMin[2] = (int(m.z / D) - 1);
//
//	/*std::cout << M.ToString() << "\n";
//	std::cout << m.ToString() << "\n";*/
//
//	//cal tri intersection point and push in map
//	Point3D interP;
//	map<pair<int, int>, sp<Cells>>::iterator MapIt;
//	map<sp<ModelTriangle>, bool> uncatchTri;
//	for (auto&& tri : *ModelTr)
//	{
//		uncatchTri.insert(pair<sp<ModelTriangle>, bool>(tri, true));
//	}
//
//	for (auto&& tri : *ModelTr)
//	{		
//			vector<Point3D> triPs;
//			
//			
//			triPs.push_back(*(tri->a));
//			triPs.push_back(*(tri->b));
//			triPs.push_back(*(tri->c));
//
//			Point3D trim = { DBL_MAX, DBL_MAX, DBL_MAX };
//			Point3D triM = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
//
//			int TriMin[] = { DBL_MAX, DBL_MAX, DBL_MAX };
//			int TriMax[] = { DBL_MAX, DBL_MAX, DBL_MAX };
//			for (auto p : triPs)
//			{
//				trim = Point3D::min(trim, p);
//				triM = Point3D::max(triM, p);
//			}
//
//			TriMax[0] = (int(triM.x / D) + 1);
//			TriMax[1] = (int(triM.y / D) + 1);
//			TriMax[2] = (int(triM.z / D) + 1);
//
//			TriMin[0] = (int(trim.x / D) - 1);
//			TriMin[1] = (int(trim.y / D) - 1);
//			TriMin[2] = (int(trim.z / D) - 1);
//
//			int x[] = { 0,0,0 };
//			for (int mode = 0; mode < 3; mode++)
//			{
//				for (x[i[mode]] = TriMin[i[mode]]; x[i[mode]] <= TriMax[i[mode]]; x[i[mode]]++)
//				{
//					for (x[j[mode]] = TriMin[j[mode]]; x[j[mode]] <= TriMax[j[mode]]; x[j[mode]]++)
//					{
//						uVector3D Dir;
//						Dir[i[mode]] = 0;
//						Dir[j[mode]] = 0;
//						Dir[k[mode]] = 1;
//
//						Point3D TriO;
//						TriO[i[mode]] = x[i[mode]] * D;
//						TriO[j[mode]] = x[j[mode]] * D;
//						TriO[k[mode]] = BoxMin[k[mode]] * D;
//						if (RayIntersectsTriangle(TriO, Dir, tri.get(), interP))
//						{
//							if (uncatchTri.find(tri) != uncatchTri.end())
//							{
//								uncatchTri.find(tri)->second = false;
//							}
//
//							if (LDI_Data[mode].find(make_pair(x[i[mode]], x[j[mode]])) == LDI_Data[mode].end())
//							{
//								Point3D O;
//								O[i[mode]] = x[i[mode]] * D;
//								O[j[mode]] = x[j[mode]] * D;
//								O[k[mode]] = 0;
//
//								sp<Cells> T = make_shared<Cells>(O, Dir);
//								LDI_Data[mode].insert(pair<pair<double, double>, sp<Cells>>(make_pair(x[i[mode]], x[j[mode]]), T));
//							}
//
//							MapIt = LDI_Data[mode].find(make_pair(x[i[mode]], x[j[mode]]));
//							if (MapIt != LDI_Data[mode].end())
//							{
//								Cells::data temp = { interP.p[mode] ,tri };
//								MapIt->second->IntsDatas.push_back(temp);
//							}
//							else
//								cout << "LDI_Date intsert has some error" << endl;
//
//						}
//					}
//				}
//			}
//
//	}
//
//
//	//sort map data
//	for (int mode = 0; mode < 3; mode++)
//	{
//		for (auto CELL : LDI_Data[mode])
//		{
//			if (CELL.second->IntsDatas.size() > 1)
//				CELL.second->SortCell();
//			/*cout << CELL.first.first << " " << CELL.first.second << "->";
//			cout << CELL.second->ToString() << endl;*/
//
//		}
//	}
//
//	for (auto a : uncatchTri)
//	{
//		if (a.second)
//		{
//			cout << a.first << endl;
//		}
//
//	}
//}
//
////////////////////////////////////////////////////////////////////////////
myModel::Model::Model() :VBOobj(), showmesh(true), dynamicshowmesh(true)
{
	colorR = 1;//R255
	colorB = 0.537;//B137
	colorG = 0.804;//G205
}

myModel::Model::Model(const myModel::Model& M) :Model()
{
	colorR = M.colorR;
	colorB = M.colorB;
	colorG = M.colorG;

	tm = M.tm;

	float* tempp = new float[M.tarr.size() * 9];
	for (tri_ind_t i = 0; i < M.tarr.size(); i++)
	{
		auto& t = *M.tarr[i];
		for (vertex_ind_t j = 0; j < 3; j++)
		{
			tempp[i * 9 + j * 3] = (float)t[j][0];
			tempp[i * 9 + j * 3 + 1] = (float)t[j][1];
			tempp[i * 9 + j * 3 + 2] = (float)t[j][2];
		}
	}
	SetModel(tempp, M.tarr.size());
	delete[] tempp;
}

myModel::Model::Model(myModel::Model&& M)
{
	colorR = M.colorR, colorG = M.colorG, colorB = M.colorB;
	
	tm = M.tm;

	MinTriSize = M.MinTriSize;
	bbox = M.bbox;
	swap(ocd, M.ocd);
	swap(ABt, M.ABt);

	showmesh = M.showmesh;
	dynamicshowmesh = M.dynamicshowmesh;

	sourceP.swap(M.sourceP);
	tarr.swap(M.tarr);
	tarrinfo.swap(M.tarrinfo);
}

myModel::Model::~Model()
{
}

myModel::Model::Model(const spvector<baseTriangle>& V) :Model()
{
	float* tempp = new float[V.size() * 9];
	for (tri_ind_t i = 0; i < V.size(); i++)
	{
		auto& t = *V.at(i);
		for (vertex_ind_t j = 0; j < 3; j++)
		{
			tempp[i * 9 + j * 3] = (float)t[j][0];
			tempp[i * 9 + j * 3 + 1] = (float)t[j][1];
			tempp[i * 9 + j * 3 + 2] = (float)t[j][2];
		}
	}

	SetModel(tempp, V.size());
	delete[] tempp;
}

myModel::Model::Model(const string& fileName) :Model()
{
	tri_ind_t numofTri = 0;
	float* tempp = LoadStlBinary(fileName, numofTri);

	if (numofTri == 0)
	{
		return;
	}

	SetModel(tempp, numofTri);
	delete[] tempp;
}

myModel::Model myModel::Model::LoadObjASCI(const string& fileName)
{
	struct objp
	{
		size_t vind;
		size_t ttind;
		size_t nind;
	};
	struct face
	{
		objp p[3];
	};

	Model result;
	vector<Point3D>& p = result.sourceP;
	vector<Vector3D> v;
	vector<pair<float, float>> tt;
	vector<Point3DNode*> pn;
	vector<face> f;

	try
	{
		ifstream fin(fileName.c_str(), ios::in | ios::binary);
		if (!fin.good())
		{
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
		}
		mtime2.clear();
		mtime2.start();

		regex vertexrx(R"(^v (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+))");
		regex facerx(R"(^f (\d+)\/(\d+)\/(\d+) (\d+)\/(\d+)\/(\d+) (\d+)\/(\d+)\/(\d+))");

		smatch sm;
		while (!fin.eof())
		{
			string temp;
			getline(fin, temp, '\n');
			if (temp[0] == 'v' && temp[1] == ' ')
			{
				regex_search(temp, sm, vertexrx);
				Point3D pp = {stod(sm[1]), stod(sm[2]), stod(sm[3])};
				Point3D::erasezero(pp);
				p.push_back(pp);
			}
// 			else if (temp[0] == 'v' && temp[1] == 't')
// 			{
// 				regex_search(temp, sm, texturerx);
// 				tt.push_back({stof(sm[1]), stof(sm[2])});
// 			}
// 			else if (temp[0] == 'v' && temp[1] == 'n')
// 			{
// 				regex_search(temp, sm, normalrx);
// 				v.push_back({stod(sm[1]), stod(sm[2]), stod(sm[3])});
// 			}
			else if (temp[0] == 'f' && temp[1] == ' ')
			{
				regex_search(temp, sm, facerx);
				face tf;
				for (int i = 0; i < 3; i++)
				{
#ifdef _WIN64
					tf.p[i] = {stoull(sm[i * 3 + 1]), stoull(sm[i * 3 + 2]), stoull(sm[i * 3 + 3])};
#else
					tf.p[i] = {stoul(sm[i * 3 + 1]), stoul(sm[i * 3 + 2]), stoul(sm[i * 3 + 3])};
#endif
				}
				f.push_back(tf);
			}
		}
		fin.close();

		pn.resize(f.size() * 3);

#ifndef useParallel
		result.tarr.reserve(f.size());
		result.tarrinfo.reserve(f.size());
		auto n = pn.begin();
		for (tri_ind_t i = 0; i < f.size(); i++)
		{
			auto t = make_shared<ModelTriangle>();

			for (vertex_ind_t j = 0; j < 3; j++)
			{
				*n++ = new Point3DNode(t->p[j] = &result.sourceP[f[i].p[j].vind - 1], i * 4 + j);
			}

			t->caln();

			result.tarr.push_back(t);
			result.tarrinfo.push_back(make_shared<triangleinfo>());
		}
#else
		result.tarr.resize(f.size());
		result.tarrinfo.resize(f.size());

		parallel_for((tri_ind_t)0, result.tarr.size(), [&](tri_ind_t i)
		{
			result.tarr[i] = make_shared<ModelTriangle>();
			result.tarrinfo[i] = make_shared<triangleinfo>();

			for (vertex_ind_t j = 0; j < 3; j++)
			{
				pn[i * 3 + j] = new Point3DNode(result.tarr[i]->p[j] = &result.sourceP[f[i].p[j].vind - 1], i * 4 + j);
			}
			result.tarr[i]->caln();
		});
#endif
		mtime2.end();
		mtime2.print("pn: ");

		result.ocd = make_shared<OctCloud>(pn);

		result.ABt.reset();

		mtime2.clear();
		mtime2.start();
		result.Triangle_update();
		mtime2.end();
		mtime2.print("tri_update: ");

		result.SetMinTriSize();

		result.bbox = *result.ocd->root;
	}
	catch (exception*)
	{
		result.tarr.clear();
		result.tarrinfo.clear();
	}

	return result;
}

void myModel::Model::SetModel(float* tempp, tri_ind_t numofTri)
{
	if (numofTri == 0)
	{
		return;
	}

	sourceP.resize(numofTri * 3);
	vector<Point3DNode*> pn(numofTri * 3);

	mtime2.clear();
	mtime2.start();
#ifndef useParallel
	tarr.reserve(numofTri);
	tarrinfo.reserve(numofTri);
	auto p = sourceP.data();
	auto n = pn.begin();
	for (tri_ind_t i = 0; i < numofTri; i++, p += 3)
	{
		auto t = make_shared<ModelTriangle>();

		for (vertex_ind_t j = 0; j < 3; j++, tempp += 3)
		{
			copy(tempp, tempp + 3, p[j].p);
			Point3D::erasezero(p[j]);
			*n++ = new Point3DNode(t->p[j] = &p[j], triid_pid(i, j));
		}

		t->caln();

		tarr.push_back(t);
		tarrinfo.push_back(make_shared<triangleinfo>());
	}
#else
	tarr.resize(numofTri);
	tarrinfo.resize(numofTri);

	copy(tempp, tempp + sourceP.size()*3, sourceP.data()->p);

	parallel_for((tri_ind_t)0, numofTri, [&](tri_ind_t i)
	{
		tarr[i] = make_shared<ModelTriangle>();
		tarrinfo[i] = make_shared<triangleinfo>();

		for (vertex_ind_t j = 0; j < 3; j++)
		{
			auto ind1 = i * 3 + j;
			auto ind2 = tempp + ind1 * 3;
			copy(ind2, ind2 + 3, sourceP[ind1].p);
			Point3D::erasezero(sourceP[ind1]);
			pn[ind1] = new Point3DNode(tarr[i]->p[j] = &sourceP[ind1], triid_pid(i, j));
		}
		tarr[i]->caln();
	});
#endif
	mtime2.end();
	mtime2.print("pn: ");

	mtime2.clear();
	mtime2.start();
	ocd = make_shared<OctCloud>(pn);
	mtime2.end();
	mtime2.print("ocd: ");

	ABt.reset();

	mtime2.clear();
	mtime2.start();
	Triangle_update();
	mtime2.end();
	mtime2.print("tri_update: ");

	SetMinTriSize();
	bbox = *ocd->root;

	/*mtime2.clear();
	mtime2.start();
	ldi = new LDI(bbox,tarr);
	mtime2.end();
	mtime2.print("LDI: ");*/
	

	mtime2.clear();
	mtime2.start();
	PAxis.CalModelPCA(&tarr,&bbox);	
	mtime2.end();
	mtime2.print("PCA: ");
	std::cout << PAxis.ToString();

	
}

void myModel::Model::SaveStlBinary(const string& fileName) const
{
	tri_ind_t numofTri = tarr.size();
	using namespace std;
	try
	{
		ofstream fout(fileName.c_str(), ios::out | ios::binary);

		if (!fout.good())
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");

		const unsigned int tmpLen = 80;
		char tmpBuff[tmpLen] = {'\0'};

		fout.write(tmpBuff, tmpLen);
		fout.write((char *)&numofTri, 4);

		char tempchar[2] = {'\0', '\0'};
		float tempf;
		char* tempc = (char*)&tempf;
		for (auto&& i : tarr)
		{
			tempf = (float)i->n.x; fout.write(tempc, 4);
			tempf = (float)i->n.y; fout.write(tempc, 4);
			tempf = (float)i->n.z; fout.write(tempc, 4);

			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					tempf = (float)i->p[j]->p[k];
					fout.write(tempc, 4);
				}
			}
			fout.write(tempchar, 2);
		}

		fout.clear();
		fout.close();
	}
	catch (exception*)
	{
		return;
	}
}

void myModel::Model::SavePointCloudASCII(const string& fileName) const
{
	using namespace std;
	try
	{
		ofstream fout(fileName.c_str(), ios::out);

		if (!fout.good())
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
		fout << fixed << setprecision(17);
		fout << ocd->pnarr.size() << "\n";
		for (auto&& i : ocd->pnarr)
		{
			for (int j = 0; j < 3; j++)
			{
				fout << i->pdata->p[j] << " ";
			}
			fout << "\n";
		}

		fout.clear();
		fout.close();
	}
	catch (exception*)
	{
		return;
	}
}

void myModel::Model::SavePointCloudBinary(const string& fileName) const
{
	using namespace std;
	try
	{
		ofstream fout(fileName.c_str(), ios::out | ios::binary);

		if (!fout.good())
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
		size_t pcount = ocd->pnarr.size();
		fout.write((char *)&pcount, 4);

		float tempf;
		char* tempc = (char*)&tempf;
		for (auto&& i : ocd->pnarr)
		{
			for (int j = 0; j < 3; j++)
			{
				tempf = (float)i->pdata->p[j];
				fout.write(tempc, 4);
			}
		}

		fout.clear();
		fout.close();
	}
	catch (exception*)
	{
		return;
	}
}

void myModel::Model::Add_tri(const spvector<baseTriangle>& bt)
{
	auto spdata = sourceP.data();

	tarr.reserve(tarr.size() + bt.size());
	sourceP.reserve(sourceP.size() + bt.size() * 9);

	for (auto&& t : tarr)
	{
		for (vertex_ind_t j = 0; j < 3; j++)
		{
			t->p[j] = sourceP.data() + (t->p[j] - spdata);
		}
	}

	for (auto&& t : bt)
	{
		auto mt = make_shared<ModelTriangle>();

		for (vertex_ind_t j = 0; j < 3; j++)
		{
			sourceP.push_back(t->at(j));
			mt->p[j] = &sourceP.back();
		}
		mt->caln();
		tarr.push_back(mt);
	}

	RebuildOctCloud();

	ABt.reset();	

	SetMinTriSize();

	bbox = *ocd->root;
}

void myModel::Model::ApplyTMat()
{
	auto cen = bbox.bboxcenter();

	mtime2.clear();
	mtime2.start();
#ifndef useParallel
	for (auto&& pn : ocd->pnarr)
	{
		auto temp = (tm.applyMat(*pn->pdata - Vector3D(cen)) + cen).tofarr();
		*pn->pdata = {temp[0], temp[1], temp[2]};
	}
#else
	parallel_for_each(ocd->pnarr.begin(), ocd->pnarr.end(), [&](Point3DNode* pn)
	{
		auto temp = Point3D(tm.applyMat(*pn->pdata - Vector3D(cen)) + cen).tofarr();
		*pn->pdata = {temp[0], temp[1], temp[2]};
	});
#endif
	mtime2.end();
	mtime2.print("amat: ");

	for (auto&& t : tarr)
	{
		t->caln();
	}
	ocd->resegmentation();
	BuildAABB();

	tm.ResetMat();
	bbox = *ocd->root;

	//mtime2.clear();
	//mtime2.start();
	//ldi = new LDI(bbox,tarr);
	//mtime2.end();
	//mtime2.print("LDI: ");
	

	mtime2.clear();
	mtime2.start();
	PAxis.CalModelPCA(&tarr, &bbox);
	mtime2.end();
	mtime2.print("PCA: ");
	std::cout << PAxis.ToString();
}

void myModel::Model::RebuildOctCloud()
{
	tarrinfo.clear();
	tarrinfo.reserve(tarr.size());

	vector<Point3DNode*> pn(tarr.size() * 3);

	auto n = pn.begin();

	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		*n++ = new Point3DNode(tarr[i]->a, triid_pid(i, 0));
		*n++ = new Point3DNode(tarr[i]->b, triid_pid(i, 1));
		*n++ = new Point3DNode(tarr[i]->c, triid_pid(i, 2));

		tarrinfo.push_back(make_shared<triangleinfo>());
	}

	ocd = make_shared<OctCloud>(pn);
	Triangle_update();
}

void myModel::Model::Triangle_update()
{

#ifndef useParallel
	for (auto&& i : ocd->pnarr)
#else
	parallel_for_each(ocd->pnarr.begin(), ocd->pnarr.end(), [&](Point3DNode* i)
#endif
	{
		auto& adjt = i->adj_tri;

		for (auto&& j : adjt)
		{
			auto& t1 = tarr[j.tid()];
			t1->p[j.pid()] = i->pdata;
			t1->np[j.pid()] = i;
		}
	}
#ifdef useParallel
	);
#endif

#ifndef useParallel
	for (tri_ind_t i = 0; i < tarr.size(); i++)
#else
	parallel_for(tri_ind_t(0), tarr.size(), [&](tri_ind_t i)
#endif
	{
		auto& adjt = tarrinfo[i]->a;
		adjt.reserve(tarr[i]->na->adj_tri.size() + tarr[i]->nb->adj_tri.size() + tarr[i]->nc->adj_tri.size() - 3);
		for (vertex_ind_t j = 0; j < 3; j++)
		{
			for (auto&& jj : tarr[i]->np[j]->adj_tri)
			{
				auto jjj = jj.tid();
				if (jjj != i)
				{
					adjt.push_back(jjj * 4 + j);
				}
			}
		}

		auto& adj2t = tarrinfo[i]->b;

		sort(adjt.begin(), adjt.end());
		adjt.erase(unique(adjt.begin(), adjt.end(), [&](const triid_pid& a, const triid_pid& b)
		{
			if (a == b)
			{
				adj2t.emplace_back(a.tid(), 3 - a.pid() - b.pid());
				return true;
			}
			return false;
		}), adjt.end());

		unsigned char checkedge = 0;
		for (auto&& j : adj2t)
		{
			checkedge |= (1 << (2 - j.pid()));
		}
		tarrinfo[i]->c = 7 - checkedge;

		tarr[i]->hasEdge = (checkedge != 7);
	}
#ifdef useParallel
	);
#endif

#ifdef DEBUG
	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (tarrinfo[i]->c != 0)
		{
			cout << i << " ";
			for (auto&& ii : tarrinfo[i]->b)
			{
				cout << "(" << ii.first << " " << ii.second << ") ";
			}
			cout << endl;
		}
	}
#endif // DEBUG

#ifndef useParallel
	for (auto&& i : ocd->pnarr)
#else
	parallel_for_each(ocd->pnarr.begin(), ocd->pnarr.end(), [&](Point3DNode* i)
#endif
	{
		auto& adjp = i->adj_Point;
		adjp.reserve(i->adj_tri.size() * 2);
		for (auto&& j : i->adj_tri)
		{
			auto& t1 = tarr[j.tid()];

			adjp.push_back(t1->np[(j.pid() + 1) % 3]);
			adjp.push_back(t1->np[(j.pid() + 2) % 3]);
		}
		sort(adjp.begin(), adjp.end());
		adjp.erase(unique(adjp.begin(), adjp.end()), adjp.end());
		adjp.shrink_to_fit();
	}
#ifdef useParallel
	);
#endif
}

void myModel::Model::SetMinTriSize()
{
	mtime3.clear();
	mtime3.start();
	MinTriSize = 10000;
	AvgTriSize = 0;
	MaxTriSize = 0;

	double esize = 1/(3.0 * tarr.size());
#ifndef useParallel
	for (auto&& t : tarr)
	{
		MinTriSize = min({(*t->b - *t->a).lengthsq(), (*t->c - *t->b).lengthsq(), (*t->a - *t->c).lengthsq(), MinTriSize});
		MaxTriSize = max({(*t->b - *t->a).lengthsq(), (*t->c - *t->b).lengthsq(), (*t->a - *t->c).lengthsq(), MaxTriSize});

		AvgTriSize += (t->lab() + t->lbc() + t->lca()) * esize;
	}
#else
	combinable<double> m([]() {return 10000.0; });
	combinable<double> M([]() {return 0.0; });
	combinable<double> A([]() {return 0.0; });
	parallel_for_each(tarr.begin(), tarr.end(), [&](const sp<ModelTriangle>& t)
	{
		m.local() = min({(*t->b - *t->a).lengthsq(), (*t->c - *t->b).lengthsq(), (*t->a - *t->c).lengthsq(), m.local()});
		M.local() = max({(*t->b - *t->a).lengthsq(), (*t->c - *t->b).lengthsq(), (*t->a - *t->c).lengthsq(), M.local()});

		A.local() += (t->lab() + t->lbc() + t->lca()) * esize;
	});

	MinTriSize = m.combine([](double l, double r) {return std::min(l, r); });
	MaxTriSize = M.combine([](double l, double r) {return std::max(l, r); });
	AvgTriSize = A.combine([](double l, double r) {return l+r; });
#endif

	MinTriSize = max(sqrt(MinTriSize), 0.05);
	MaxTriSize = sqrt(MaxTriSize);
	mtime3.end();
	mtime3.print("find tri size: ");
}

void myModel::Model::BuildAABB()
{
	auto tarrsize = tarr.size();
	vector<nodedata> AABBnodes(tarrsize);

#ifndef useParallel
#pragma omp parallel for
	for (tri_ind_t i = 0; i < tarrsize; i++)
	{
		AABBnodes[i] = BoundBox::genbbox(*tarr[i]);
	}
#else
	parallel_for((tri_ind_t)0, tarrsize, [&](tri_ind_t i)
	{
		AABBnodes[i] = BoundBox::genbbox(*tarr[i]);
	});
#endif

	ABt = make_shared<AABBtree>(move(AABBnodes));
}

tuple<bool, char, char, double> myModel::Model::checkDegenerationTriangle(const myModel::ModelTriangle& T, double llimit) const
{
	bitset<3> lflag;
	bitset<3> dflag;

	//unsigned int lflag = 0;
	//unsigned int dflag = 0;
	double lab = T.lab();
	double lbc = T.lbc();
	double lca = T.lca();

	lflag[0] = lab < llimit;
	lflag[1] = lbc < llimit;
	lflag[2] = lca < llimit;

	// 	if (lab < llimit)
	// 	{
	// 		//lflag += 1;
	// 		lflag[0] = true;
	// 	}
	// 	if (lbc < llimit)
	// 	{
	// 		//lflag += 2;
	// 		lflag[1] = true;
	// 	}
	// 	if (lca < llimit)
	// 	{
	// 		//lflag += 4;
	// 		lflag[2] = true;
	// 	}

	if (lflag.any()/*lflag>0*/)//lflag == 1 || lflag == 2 || lflag == 4 || lflag == 7)
	{
		double minlength = min({lab, lbc, lca});//(lab < lbc) ? ((lab < lca) ? lab : lca) : ((lbc < lca) ? lbc : lca);
												//cout << "ml, " << minlength << endl;
		return make_tuple(true, 0, char(lflag.to_ulong()), minlength);
	}
	else
	{
		double dab = T.dab();
		double dbc = T.dbc();
		double dca = T.dca();

		dflag[0] = dab < llimit * 10;
		dflag[1] = dbc < llimit * 10;
		dflag[2] = dca < llimit * 10;

		// 		if (dab < llimit*10)
		// 		{
		// 			//dflag += 1;
		// 			dflag[0] = true;
		// 		}
		// 		if (dbc < llimit*10)
		// 		{
		// 			//dflag += 2;
		// 			dflag[1] = true;
		// 		}
		// 		if (dca < llimit*10)
		// 		{
		// 			//dflag += 4;
		// 			dflag[2] = true;
		// 		}

		if (dflag.count() == 2/*dflag == 3 || dflag == 5 || dflag == 6*/)//dabc,dcab,dbca
		{
			double mindeg = min({dab, dbc, dca});//;(dab < dbc) ? ((dab < dca) ? dab : dca) : ((dbc < dca) ? dbc : dca);
												 //cout << "md, " << mindeg << endl;
			return make_tuple(true, 1, char(dflag.to_ulong()), mindeg);
		}
	}

	//double rlimit = 150;
	//double rab = lbc / lab + lca / lab;
	//double rbc = lab / lbc + lca / lbc;
	//double rca = lab / lca + lbc / lca;
	//if (rab > rlimit)
	//{
	//	//cout << "rab, " << rab << endl;
	//	return make_tuple(true, 0, 1, lab);
	//}
	//if (rbc > rlimit)
	//{
	//	//cout << "rbc, " << rbc << endl;
	//	return make_tuple(true, 0, 2, lbc);
	//}
	//if (rca > rlimit)
	//{
	//	//cout << "rca, " << rca << endl;
	//	return make_tuple(true, 0, 4, lca);
	//}
	return make_tuple(false, 0, 0, 0);
}

spvector<baseTriangle> myModel::Model::searchDegenerationTriangle(double llimit) const
{
	spvector<baseTriangle> v;

	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (get<0>(checkDegenerationTriangle(*tarr[i], llimit)))
		{
			v.push_back(tarr[i]);
		}
	}
	return v;
}

void myModel::Model::removeDegenerationTriangle(double llimit)
{
	char* isDTri = new char[tarr.size()];
	size_t DtriCount = 0;

	using ty = tuple<tri_ind_t, char, double>;
	vector<ty> edge_collapse;
	vector<ty> edge_swapping;

	reference_wrapper<vector<ty>> edged[] = {edge_collapse, edge_swapping};
	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		isDTri[i] = 0;
		bool checkt;
		char dmethod;
		char dflag;
		double comp;
		tie(checkt, dmethod, dflag, comp) = checkDegenerationTriangle(*tarr[i], llimit);
		if (checkt /*&& dmethod != 1*/)
		{
			DtriCount++;
			edged[dmethod].get().emplace_back(i, dflag, comp);
			isDTri[i] = dmethod + 1;
		}
	}

	auto ls = [](const ty& d1, const ty& d2) {return get<2>(d1) < get<2>(d2); };

	std::sort(edge_collapse.begin(), edge_collapse.end(), ls);
	std::sort(edge_swapping.begin(), edge_swapping.end(), ls);

	for (size_t i = 0; i < edge_swapping.size(); i++)
	{
		tri_ind_t ind;
		int mflag;
		tie(ind, mflag, ignore) = edge_swapping[i];

		if (isDTri[ind] != 2)
		{
			continue;
		}

		if (mflag == 3 || mflag == 5 || mflag == 6)
		{
			auto t = tarr[ind];
			auto adjt = &(tarrinfo[ind]->b);

			vertex_ind_t startpid = -1;
			if (mflag == 3)
			{
				startpid = 2;
			}
			else if (mflag == 5)
			{
				startpid = 1;
			}
			else if (mflag == 6)
			{
				startpid = 0;
			}
			Point3D startp = *t->p[startpid];
			for (auto&& it : *adjt)
			{
				if (it.pid() == startpid)
				{
					cout << "oh ";
				}
			}
			cout << "\n";
			auto antri = find_if(adjt->begin(), adjt->end(), [&](const triid_pid& t) {return t.pid() == startpid; });
			auto antriind = antri->tid();

			auto t2 = tarr[antriind];
			auto adjt2 = &tarrinfo[antriind]->b;
			auto endpn = find_if(adjt2->begin(), adjt2->end(), [&](const triid_pid& t) {return t.tid() == ind; });
			auto endpid = endpn->pid();

			if (t->p[startpid] == t2->p[endpid])
			{
				cout << "start same with end\n";
				continue;
			}

			Point3D endp = *t2->p[endpid];

			if (isDTri[antriind] == 2)
			{
				if ((endp - startp).length() < llimit)
				{
					Point3D cen = (endp + startp) / 2;
					auto& t2edpadjtri = t2->np[endpid]->adj_tri;
					for (auto&& tind : t2edpadjtri)
					{
						if (ind != tind.tid())
						{
							tarr[tind.tid()]->p[tind.pid()] = t->p[startpid];
							tarr[tind.tid()]->np[tind.pid()] = t->np[startpid];
							t->np[startpid]->adj_tri.push_back(tind);
						}
					}
					*t->np[startpid]->pdata = cen;
				}
				isDTri[ind] = 1;
				isDTri[antriind] = 1;
				continue;
			}
			auto& v1 = t->np[(startpid + 1) % 3]->adj_tri;
			auto& v2 = t->np[(startpid + 2) % 3]->adj_tri;
			v1.erase(remove_if(v1.begin(), v1.end(), [&](const triid_pid& t) {return t.tid() == antriind; }), v1.end());
			v2.erase(remove_if(v2.begin(), v2.end(), [&](const triid_pid& t) {return t.tid() == ind; }), v2.end());

			t->np[startpid]->adj_tri.emplace_back(antriind, (endpid + 2) % 3);
			t2->np[endpid]->adj_tri.emplace_back(ind, (startpid + 2) % 3);

			t->p[(startpid + 2) % 3] = t2->p[endpid];
			t->np[(startpid + 2) % 3] = t2->np[endpid];
			t2->p[(endpid + 2) % 3] = t->p[startpid];
			t2->np[(endpid + 2) % 3] = t->np[startpid];

			isDTri[ind] = 0;
			DtriCount--;
		}
	}

	for (size_t i = 0; i < edge_collapse.size(); i++)
	{
		tri_ind_t ind;
		int mflag;

		tie(ind, mflag, ignore) = edge_collapse[i];

		ModelTriangle* t = tarr[ind].get();
		if (tarr[ind]->area() < 1e-15)
		{
			if (isDTri[ind] == 0)
			{
				isDTri[ind] = 1;
				DtriCount++;
			}
			continue;
		}

		if (mflag == 1)
		{
			Point3D cen = (*t->a + *t->b) / 2;
			for (auto&& ind : t->nb->adj_tri)
			{
				// 				(&tarr[ind.first]->a)[ind.second] = t->a;
				// 				(&tarr[ind.first]->na)[ind.second] = t->na;
				tarr[ind.tid()]->p[ind.pid()] = t->a;
				tarr[ind.tid()]->np[ind.pid()] = t->na;
				t->na->adj_tri.push_back(ind);
			}
			*t->na->pdata = cen;
		}
		else if (mflag == 2)
		{
			Point3D cen = (*t->b + *t->c) / 2;
			for (auto&& ind : t->nc->adj_tri)
			{
				// 				(&tarr[ind.first]->a)[ind.second] = t->b;
				// 				(&tarr[ind.first]->na)[ind.second] = t->nb;
				tarr[ind.tid()]->p[ind.pid()] = t->b;
				tarr[ind.tid()]->np[ind.pid()] = t->nb;
				t->nb->adj_tri.push_back(ind);
			}
			*t->nb->pdata = cen;
		}
		else if (mflag == 4)
		{
			Point3D cen = (*t->a + *t->c) / 2;
			for (auto&& ind : t->nc->adj_tri)
			{
				// 				(&tarr[ind.first]->a)[ind.second] = t->a;
				// 				(&tarr[ind.first]->na)[ind.second] = t->na;
				tarr[ind.tid()]->p[ind.pid()] = t->a;
				tarr[ind.tid()]->np[ind.pid()] = t->na;
				t->na->adj_tri.push_back(ind);
			}
			*t->na->pdata = cen;
		}
		else// if (mflag == 7)
		{
			Point3D cen = (*t->a + *t->b + *t->c) / 3;
			for (auto&& ind : t->nb->adj_tri)
			{
				// 				(&tarr[ind.first]->a)[ind.second] = t->a;
				// 				(&tarr[ind.first]->na)[ind.second] = t->na;
				tarr[ind.tid()]->p[ind.pid()] = t->a;
				tarr[ind.tid()]->np[ind.pid()] = t->na;
				t->na->adj_tri.push_back(ind);
			}
			for (auto&& ind : t->nc->adj_tri)
			{
				// 				(&tarr[ind.first]->a)[ind.second] = t->a;
				// 				(&tarr[ind.first]->na)[ind.second] = t->na;
				tarr[ind.tid()]->p[ind.pid()] = t->a;
				tarr[ind.tid()]->np[ind.pid()] = t->na;
				t->na->adj_tri.push_back(ind);
			}
			*t->na->pdata = cen;
		}
	}

	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (isDTri[i] == 0)
		{
			if (tarr[i]->a == tarr[i]->b || tarr[i]->c == tarr[i]->b || tarr[i]->a == tarr[i]->c)
			{
				isDTri[i] = 1;
				DtriCount++;
			}
			//  			if (tarr[i]->area() < 1e-15)
			//  			{			
			//  				isDTri[i] = 1;
			//  				DtriCount++;
			//  			}
		}
	}

	size_t numoftri = tarr.size() - DtriCount;
	//numofTri = numofTri - DtriCount;
	float* tempp = new float[numoftri * 9];
	float* p = tempp;
	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (isDTri[i] != 0)
		{
			continue;
		}

		p[0] = (float)tarr[i]->a->x;
		p[1] = (float)tarr[i]->a->y;
		p[2] = (float)tarr[i]->a->z;
		p[3] = (float)tarr[i]->b->x;
		p[4] = (float)tarr[i]->b->y;
		p[5] = (float)tarr[i]->b->z;
		p[6] = (float)tarr[i]->c->x;
		p[7] = (float)tarr[i]->c->y;
		p[8] = (float)tarr[i]->c->z;
		p += 9;
	}

	delete[] isDTri;

	ocd.reset();
	sourceP.clear();
	tarr.clear();
	tarrinfo.clear();

	SetModel(tempp, numoftri);
	delete[] tempp;
}

pair<bitset<3>, bitset<3>> Determiniant6(const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D, const Point3D& E, const Point3D& F)
{
	bitset<3> flag, flagzero(0);

	double dd = Vector3D::dot(A - D, Vector3D::cross(B - D, C - D));
	double ee = Vector3D::dot(A - E, Vector3D::cross(B - E, C - E));
	double ff = Vector3D::dot(A - F, Vector3D::cross(B - F, C - F));
	flag[0] = dd <= 0;
	flag[1] = ee <= 0;
	flag[2] = ff <= 0;


	if (flag[2])
	{
		flag.flip();
	}

	if (flag.any())
	{
		flagzero[0] = abs(dd) < 1e-15;//FLT_EPSILON;
		flagzero[1] = abs(ee) < 1e-15;//FLT_EPSILON;
		flagzero[2] = abs(ff) < 1e-15;//FLT_EPSILON;
		if (flagzero.count() > 1)
		{
			//cout << __FUNCTION__ << __LINE__ << " zero~~~ " << dd << " " << ee << " " << ff << endl;
			flag = flagzero;
			if (flag[2])
			{
				flag.flip();
			}
		}
	}

	return make_pair(flag, flagzero);
}

pair<bitset<2>, bitset<2>> Determiniant4(const Point3D& A, const Point3D& B, const Point3D& C/*, const Point3D& D*/, const Point3D& E, const Point3D& F)
{
	bitset<2> flag, flagzero(0);

	double ee = Vector3D::dot(A - E, Vector3D::cross(B - E, C - E));
	double ff = Vector3D::dot(A - F, Vector3D::cross(B - F, C - F));

	flag[0] = ee <= 0;
	flag[1] = ff <= 0;

	if (flag[1])
	{
		flag.flip();
	}

	if (flag.any())
	{
		flagzero[0] = abs(ee) < 1e-15;//FLT_EPSILON;
		flagzero[1] = abs(ff) < 1e-15;//FLT_EPSILON;
									  // 		if (flagzero.count() > 0)
									  // 		{
									  // 			cout << __FUNCTION__ << __LINE__ << " zero~~~ " << ee << " " << ff << endl;
									  // 		}
	}

	return make_pair(flag, flagzero);
}

Point3D TxL(const myModel::ModelTriangle& T, const myMath::Line3D& L)
{
	auto E1 = *T.b - *T.a;
	auto E2 = *T.c - *T.a;
	auto TT = L.p - *T.a;
	auto PP = Vector3D::cross(L.n, E2);
	auto QQ = Vector3D::cross(TT, E1);
	auto tt = 1 / Vector3D::dot(PP, E1);

	auto result = L.p + L.n*(Vector3D::dot(QQ, E2)*tt);
	Point3D::erasezero(result);

	return result;
}

//pair<bool, QSegment3D*> myModel::Model::checkTriangleCross(const myModel::ModelTriangle& T1, const myModel::ModelTriangle& T2)
//{
//	auto nocross = make_pair(false, nullptr);
//
//	Point3D& T1a = *T1.a;
//	Point3D& T1b = *T1.b;
//	Point3D& T1c = *T1.c;
//	Point3D& T2a = *T2.a;
//	Point3D& T2b = *T2.b;
//	Point3D& T2c = *T2.c;
//
//	auto flag_sameside_T2 = Determiniant6(T1a, T1b, T1c, T2a, T2b, T2c);
//
//	if (flag_sameside_T2.first.none())
//	{
//		//std::cout << "nocross\n";
//		return nocross;
//	}
//
//	auto flag_sameside_T1 = Determiniant6(T2a, T2b, T2c, T1a, T1b, T1c);
//
//	if (flag_sameside_T1.first.none())
//	{
//		//std::cout << "nocross\n";
//		return nocross;
//	}
//
//	int start1 = flag_sameside_T1.first.to_ulong() - 1;
//	int	start2 = flag_sameside_T2.first.to_ulong() - 1;
//
//	Point3D P1[3];
//	Point3D P2[3];
//
//	for (int i = 0; i < 3; i++)
//	{
//		P1[i] = *T1.p[(start1 + i) % 3];
//		P2[i] = *T2.p[(start2 + i) % 3];
//	}
//
//	auto temp1 = P1[0] - P2[0];
//	auto temp2 = -temp1;
//
//	if (Vector3D::dot(temp1, T1.n) < 0)
//	{
//		swap(P1[1], P1[2]);
//	}
//	if (Vector3D::dot(temp2, T2.n) < 0)
//	{
//		swap(P2[1], P2[2]);
//	}
//
//	bool resultb = Determinant(P1[0], P1[1], P2[0], P2[1]) <= 0 && Determinant(P1[0], P1[2], P2[2], P2[0]) <= 0;
//	if (!resultb)
//	{
//		return nocross;
//	}
//
//	auto b1 = Determinant(temp1, P1[2] - P2[0], P2[1] - P2[0]) <= 0 ? TxL(T1, Line3D(P2[1], P2[0])) : TxL(T2, Line3D(P1[2], P1[0]));
//	auto b2 = Determinant(temp1, P1[1] - P2[0], P2[2] - P2[0]) <= 0 ? TxL(T2, Line3D(P1[1], P1[0])) : TxL(T1, Line3D(P2[2], P2[0]));
//
//	if ((b1 - b2).length() < FLT_EPSILON)
//	{
//		//cout << __FUNCTION__ << __LINE__ << " " << (b1 - b2).length() << " b1b2tooshort\n";
//		return nocross;
//	}
//
//	return make_pair(true, new QSegment3D(b1, b2));
//}
//
//pair<bool, QSegment3D*> myModel::Model::checkTriangleCross_onesame(const myModel::ModelTriangle& T1, const myModel::ModelTriangle& T2)
//{
//	auto nocross = make_pair(false, nullptr);
//
//	int pt1 = -1;
//	int pt2 = -1;
//	for (int i = 0; i < 3; i++)
//	{
//		bool flag = false;
//		for (int j = 0; j < 3; j++)
//		{
//			if (T1.np[i] == T2.np[j])
//			{
//				pt1 = i + 1;
//				pt2 = j + 1;
//				flag = true;
//				break;
//			}
//		}
//		if (flag)
//		{
//			break;
//		}
//	}
//
//	Point3D P1[3];
//	Point3D P2[3];
//
//	for (int i = 0; i < 3; i++)
//	{
//		P1[i] = *T1.p[(pt1 + i) % 3];
//		P2[i] = *T2.p[(pt2 + i) % 3];
//	}
//
//	auto flag_sameside_T2 = Determiniant4(P1[2], P1[0], P1[1]/*, P2[2]*/, P2[0], P2[1]);
//
//	if (flag_sameside_T2.first.none() || flag_sameside_T2.second.all())
//	{
//		return nocross;
//	}
//
//	auto flag_sameside_T1 = Determiniant4(P2[2], P2[0], P2[1]/*, P1[2]*/, P1[0], P1[1]);
//
//	if (flag_sameside_T1.first.none() || flag_sameside_T1.second.all())
//	{
//		return nocross;
//	}
//
//	auto n1 = T1.n;
//	auto n2 = T2.n;
//
//	if (flag_sameside_T2.second[0])
//	{
//		swap(P2[0], P2[1]);
//		n2 = -n2;
//	}
//	if (flag_sameside_T1.second[0])
//	{
//		swap(P1[0], P1[1]);
//		n1 = -n1;
//	}
//
//	auto temp1 = P1[0] - P2[0];
//	auto temp2 = -temp1;
//
//	if (Vector3D::dot(temp1, n1) < 0)
//	{
//		swap(P1[1], P1[2]);
//	}
//	if (Vector3D::dot(temp2, n2) < 0)
//	{
//		swap(P2[1], P2[2]);
//	}
//
//	bool resultb = Determinant(P1[0], P1[1], P2[0], P2[1]) < 0 && Determinant(P1[0], P1[2], P2[2], P2[0]) < 0;
//	if (!resultb)
//	{
//		return nocross;
//	}
//
//	auto b1 = Determinant(temp1, P1[2] - P2[0], P2[1] - P2[0]) <= 0 ? TxL(T1, Line3D(P2[1], P2[0])) : TxL(T2, Line3D(P1[2], P1[0]));
//	auto b2 = Determinant(temp1, P1[1] - P2[0], P2[2] - P2[0]) <= 0 ? TxL(T2, Line3D(P1[1], P1[0])) : TxL(T1, Line3D(P2[2], P2[0]));
//
//	if ((b1 - b2).length() < FLT_EPSILON)
//	{
//		//cout << __FUNCTION__ << __LINE__ << " " << (b1 - b2).length() << " b1b2tooshort\n";
//		return nocross;
//	}
//
//	return make_pair(true, new QSegment3D(b1, b2));
//}

myModel::Model::TriCrossResult myModel::Model::searchTriangleCross()
{
	mtime3.clear();
	mtime3.start();
	if (ABt == nullptr)
	{
		BuildAABB();
	}
	mtime3.end();
	mtime3.print("BuildAABB: ");
	mtime3.clear();
	mtime3.start();
	auto testind = ABt->returnind();
	mtime3.end();
	mtime3.print("returnind: ");

	return searchTriangleCross(this, testind);
}

myModel::Model::TriCrossResult myModel::Model::searchTriangleCross(const Model* rhs, const inddatas& rhsind) const
{
	const bool onemodel = (rhs == this);
	TriCrossResult sresult;

	struct ttest
	{
		tri_ind_t k;
		bool operator()(const triid_pid& aaa) const
		{
			return k == aaa.tid();
		}
	} /*tk*/;

	vector<tri_ind_t> crossind;

#ifndef useParallel
	spvector<QSegment3D> qs;
	auto& crossind2 = crossind;
	for (tri_ind_t j = 0; j < tarr.size(); j++)
#else
	concurrent_vector<sp<QSegment3D>> qs;
	concurrent_vector<tri_ind_t> crossind2;
	parallel_for(tri_ind_t(0), tarr.size(), [&](tri_ind_t j)
#endif
	{
		ModelTriangle& t1 = *tarr[j];
		auto& adjtrim = tarrinfo[j]->a;
		auto& adj2trim = tarrinfo[j]->b;
		ttest tk;
		bool flag = false;
		for (auto&& k : rhsind[j])
		{
			auto cTC = &checkTriangleCross;

			if (onemodel)
			{
				tk.k = k;

				if (find_if(adjtrim.begin(), adjtrim.end(), tk) != adjtrim.end())
				{
					if (find_if(adj2trim.begin(), adj2trim.end(), tk) == adj2trim.end())
					{
						cTC = &checkTriangleCross_onesame;
					}
					else
					{
						continue;
					}
				}
			}

			auto testresult = cTC(t1, *rhs->tarr[k]);

			if (testresult.first)
			{
				shared_ptr<QSegment3D> ts(testresult.second);
				ts->tind1 = j;
				ts->tind2 = k;
				qs.push_back(ts);
				flag = true;
				if (onemodel)
				{
					crossind2.push_back(k);
				}
			}
			else
			{
				delete testresult.second;
			}
		}
		if (flag)
		{
			crossind2.push_back(j);
		}
	}
#ifdef useParallel
	);

	crossind.resize(crossind2.size());
	copy(crossind2.begin(), crossind2.end(), crossind.begin());
#endif
	sort(crossind.begin(), crossind.end());
	crossind.erase(unique(crossind.begin(), crossind.end()), crossind.end());

	for (auto&& i : crossind)
	{
		sresult.emplace(i, spvector<QSegment3D>());
	}
	for (auto&& ss : qs)
	{
		sresult.at(ss->tind1).push_back(ss);
		if (onemodel)
		{
			sresult.at(ss->tind2).push_back(ss);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	hash<Point3D> hp;
	typedef pair<size_t, reference_wrapper<Point3D>> phv;
	vector<phv> phvarr;
	for (auto&& sit : sresult)
	{
		for (auto&& j : sit.second)
		{
			phvarr.emplace_back(hp(j->pstart), ref(j->pstart));
			phvarr.emplace_back(hp(j->pend), ref(j->pend));
		}
	}
	sort(phvarr.begin(), phvarr.end(), [](const phv& p1, const phv& p2) {return p1.first < p2.first; });
	unique(phvarr.begin(), phvarr.end(), [](const phv& p1, const phv& p2)
	{
		if (p1.first == p2.first && p1.second.get().equal_tof(p2.second.get()))
		{
			p2.second.get() = p1.second.get();
			return true;
		}
		return false;
	});

	return sresult;
}

myModel::Model::TriSeparateResult myModel::Model::separateTriangle(const TriCrossResult& crossresult)
{
	TriSeparateResult mresult;
	if (crossresult.size() == 0)
	{
		return mresult;
	}

	for (auto&& i : crossresult)
	{
		auto stt = SubTriangles(tarr[i.first].get(), i.second).Triangulation();
		if (!stt.empty())
		{
			mresult.emplace(i.first, move(stt));
		}
	}

	return mresult;

	//  	spvector<baseTriangle>* V = new spvector<baseTriangle>();
	//  	V->reserve(tarr.size());
	//  	for (unsigned int i = 0; i < tarr.size(); i++)
	//  	{
	//  		auto pa = mresult.find(i);
	//  		if (pa != mresult.end())
	//  		{
	//  			for (auto&& j : *pa->second)
	//  			{
	//  				V->push_back(j);
	//  			}
	//  		}
	//  		//  		else
	//  		//  		{
	//  		//  			V->push_back(tarr[i]);
	//  		//  		}
	//  	}
	//  	myModel::SaveStlBinary("..\\..\\testsample\\testsssoh3.stl", V);
	//  	delete V;


	/*vector<reference_wrapper<spvector<ModelTriangle>>> testtemp1;
	testtemp1.push_back(ref(tarr));
	vector<TriSeparateResult> testtemp;
	testtemp.push_back(move(mresult));
	repairtsr(testtemp1, testtemp);
	return testtemp[0];*/
}

spvector<Segment3D> myModel::Model::searchTriangleCross_fordisplay(const TriCrossResult& x)
{
	spvector<Segment3D> sresult;

	for (auto&& i : x)
	{
		for (auto&& j : i.second)
		{
			if (i.first == j->tind1)
			{
				sresult.push_back(j);
			}
		}
	}

	return sresult;
}

spvectors<Segment3D> myModel::Model::searchTriangleCross_fordisplay2(const TriCrossResult& x)
{
	spvectors<Segment3D> sresult;

	spvector<QSegment3D> xx;

	for (auto&& i : x)
	{
		for (auto&& j : i.second)
		{
			if (i.first == j->tind1)
			{
				xx.push_back(j);
			}
		}
	}
	hash<Point3D> hp;
	while (!xx.empty())
	{
		spvector<Segment3D> sc;
		spvector<QSegment3D> tc;
		auto starts = xx.back();
		xx.pop_back();

		tc.push_back(starts);

		while (!tc.empty())
		{
			auto currs = tc.back();
			tc.pop_back();
			sc.push_back(currs);

			auto hp1 = hp(currs->pstart);
			auto hp2 = hp(currs->pend);

			auto check = [&](const shared_ptr<QSegment3D>& s)
			{
				auto hps = hp(s->pstart);
				
				if (hps == hp1)
				{
					return s->pstart.equal_tof(currs->pstart);
				}
				if (hps == hp2)
				{
					return s->pstart.equal_tof(currs->pend);
				}

				auto hpe = hp(s->pend);

				if (hpe == hp1)
				{
					return s->pend.equal_tof(currs->pstart);
				}
				if (hpe == hp2)
				{
					return s->pend.equal_tof(currs->pend);
				}
				return false;
			};

			auto it = find_if(xx.begin(), xx.end(), check);

			while (it != xx.end())
			{
				tc.push_back(*it);
				if (it != --xx.end())
				{
					*it = xx.back();
					xx.pop_back();
					it = find_if(it, xx.end(), check);
				}
				else
				{
					xx.pop_back();
					break;
				}
			}
		}

		sresult.push_back(sc);
	}

	return sresult;
}

spvector<baseTriangle> myModel::Model::subtri_fordisplay(const TriSeparateResult& st)
{
	spvector<baseTriangle> tresult;
	size_t tsize = 0;

	for (auto&& i : st)
	{
		tsize += i.second.size();
	}

	tresult.reserve(tsize);

	for (auto&& i : st)
	{
		tresult.insert(tresult.end(), i.second.begin(), i.second.end());
	}

	return tresult;
}

spvectors<Segment3D> myModel::Model::searchEdgecontour_fordisplay() const
{
	return searchEdgecontour_fordisplay(*this, searchEdgecontour());
}

spvectors<Segment3D> myModel::Model::searchEdgecontour_fordisplay(const Model& m, const vectors<Edgeind>& vei)
{
	spvectors<Segment3D> sed;
	sed.reserve(vei.size());
	for (auto&& i : vei)
	{
		spvector<Segment3D> tempss;
		tempss.reserve(i.size());
		for (auto&& j : i)
		{
			tempss.push_back(make_shared<Segment3D>(m.tarr[j.a]->at(j.b), m.tarr[j.a]->at(j.c)));
		}
		sed.push_back(move(tempss));
	}
	return sed;
}

struct ethv
{
	size_t hv;
	float lensq;
	shared_ptr<baseTriangle> t;
	bool isslim;

	ethv() {};
	ethv(const size_t hv_, const float lensq_, const shared_ptr<baseTriangle>& t_ = nullptr, const bool isslim_ = false) :hv(hv_), lensq(lensq_), t(t_), isslim(isslim_) {};
	bool operator<(const ethv& e2) const
	{
		return hv < e2.hv || (!(e2.hv < hv) && lensq < e2.lensq);
	}
	bool check_slim() const
	{
		return t != nullptr && isslim;
	}
};

void myModel::Model::repairtsr(vector<reference_wrapper<spvector<ModelTriangle>>>& modelarr, vector<Model::TriSeparateResult>& tsr)
{
	auto counttsr = [&]()
	{
		size_t cc = 0;

		for (auto&& ts : tsr)
		{
			for (auto&& t : ts)
			{
				cc += t.second.size();
			}
		}

		//cout << "cc: " << cc << endl;
		return cc;
	};

	size_t cc1 = counttsr();
	//////////////////////////////////////////////////////////////////////////
	typedef pair<size_t, sp<baseTriangle>> tttt;
	vector<tttt> temparr;

	hash<sp<baseTriangle>> ht;
	for (auto&& mm : tsr)
	{
		for (auto&& pa : mm)
		{
			for (auto&& j : pa.second)
			{
				temparr.push_back({ht(j), j});
			}
		}
	}

	sort(temparr.begin(), temparr.end(), [](const tttt& a, const tttt& b) { return a.first < b.first; });
	unique(temparr.begin(), temparr.end(), [](const tttt& a, const tttt& b)
	{
		if (a.first == b.first)
		{
			a.second->at(0).x = nan("");
			b.second->at(0).x = nan("");
			return true;
		}
		return false;
	});
	//////////////////////////////////////////////////////////////////////////
	size_t tsrsize = 0;
	for (auto&& mm : tsr)
	{
		tsrsize += mm.size();
	}

	vector<ethv> ethvarrtemp;
	vector<ethv> normalethvarr;
	normalethvarr.reserve(tsrsize * 3);

	hash<myMath::Point3D> hp;
	for (size_t i = 0; i < tsr.size(); i++)
	{
		for (auto&& pa : tsr[i])
		{
			if (pa.second.size() > 0)
			{
				auto& t = modelarr[i].get().at(pa.first);
				size_t hpp[3] = {hp(t->at(0)), hp(t->at(1)), hp(t->at(2))};
				normalethvarr.emplace_back(hpp[1] ^ hpp[2], (float)(t->at(1) - t->at(2)).lengthsq());
				normalethvarr.emplace_back(hpp[0] ^ hpp[2], (float)(t->at(0) - t->at(2)).lengthsq());
				normalethvarr.emplace_back(hpp[0] ^ hpp[1], (float)(t->at(0) - t->at(1)).lengthsq());

				for (auto&& j : pa.second)
				{
					if (isnan(j->at(0).x))
					{
						continue;
					}
					size_t hpp[3] = {hp(j->at(0)), hp(j->at(1)), hp(j->at(2))};
					bool slimflag = j->isslim();
					ethvarrtemp.emplace_back(hpp[1] ^ hpp[2], (float)(j->at(1) - j->at(2)).lengthsq(), j, slimflag);
					ethvarrtemp.emplace_back(hpp[0] ^ hpp[2], (float)(j->at(0) - j->at(2)).lengthsq(), j, slimflag);
					ethvarrtemp.emplace_back(hpp[0] ^ hpp[1], (float)(j->at(0) - j->at(1)).lengthsq(), j, slimflag);
				}
			}
		}
	}

	sort(ethvarrtemp.begin(), ethvarrtemp.end());
	sort(normalethvarr.begin(), normalethvarr.end());
	normalethvarr.erase(unique(normalethvarr.begin(), normalethvarr.end(), [](ethv& e1, ethv& e2)
	{
		if (e1.hv == e2.hv && e1.lensq == e2.lensq)
		{
			e1.isslim = e2.isslim = true;
			return true;
		}
		return false;
	}), normalethvarr.end());
	normalethvarr.erase(remove_if(normalethvarr.begin(), normalethvarr.end(), [](const ethv& e1) {return e1.isslim; }), normalethvarr.end());
	//////////////////////////////////////////////////////////////////////////
	vector<ethv> ethvarr(normalethvarr.size() + ethvarrtemp.size());

	merge(normalethvarr.begin(), normalethvarr.end(), ethvarrtemp.begin(), ethvarrtemp.end(), ethvarr.begin());

	typedef vector<ethv>::iterator ethvit;
	auto Iterative_for_each = [&ethvarr](size_t(*fn)(ethvit, ethvit))
	{
		size_t innercount = 0;
		for (auto it = ethvarr.begin(); it != ethvarr.end();)
		{
			auto nit = upper_bound(it, ethvarr.end(), *it);
			innercount += fn(it, nit);
			it = nit;
		}
		return innercount;
	};
	auto count_unnormalt = [](ethvit it, ethvit nit)->size_t
	{
		if (std::distance(it, nit) == 1)
		{
			if (it->t != nullptr)
			{
				it->t->at(0).x = nan("");
				return 1;
			}
			else
			{
				cout << __FILE__ << __LINE__ << "wtf\n";
			}
		}

		return 0;
	};
	auto count_slim = [](ethvit it, ethvit nit)->size_t
	{
		auto itdis = std::distance(it, nit);
		if (itdis >= 4 && itdis % 2 == 0)
		{
			vector<ethv> count_temp;
			auto ii = back_inserter(count_temp);
			copy_if(it, nit, ii, [](const ethv& e1) {return e1.check_slim(); });
			if (count_temp.size() > 1 && count_temp.size() % 2 == 0)
			{
				for (auto&& it2 : count_temp)
				{
					it2.t->at(0).x = nan("");
				}
				return count_temp.size();
			}
		}
		else if (itdis == 3)
		{
			for (auto i = it; i != nit; i++)
			{
				if (i->check_slim())
				{
					i->t->at(0).x = nan("");
				}
			}
			return 1;
		}
		return 0;
	};
	auto ethvarrremover = [&ethvarr]()
	{
		ethvarr.erase(remove_if(ethvarr.begin(), ethvarr.end(), [](const ethv& e1)
		{
			return e1.t != nullptr && isnan(e1.t->at(0).x);
		}), ethvarr.end());
	};

	size_t slimcounttemp = 0;
	while (true)
	{
		size_t unnormaltcount = Iterative_for_each(count_unnormalt);

		if (unnormaltcount == 0)
		{
			size_t slimcount = Iterative_for_each(count_slim);

			if (slimcount == slimcounttemp)
			{
				break;
			}
			else
			{
				slimcounttemp = slimcount;
			}

			//cout << "slimtri: " << slimcount << endl;
		}

		//cout << "removed: " << unnormaltcount << endl;

		ethvarrremover();
	}

	for (auto&& mm : tsr)
	{
		for (auto&& nn : mm)
		{
			auto& vv = nn.second;

			vv.erase(remove_if(vv.begin(), vv.end(), [](const shared_ptr<baseTriangle>& t)
			{
				if (isnan(t->at(0).x))
				{
					return true;
				}
				if (t->area() == 0)
				{
					return true;
				}
				return false;
			}), vv.end());
		}
	}
	//////////////////////////////////////////////////////////////////////////
	size_t cc2 = counttsr();
	if (cc1 != cc2)
	{
		MessageBox(
			GetActiveWindow(),
			(LPCWSTR)L"這還有用",
			(LPCWSTR)L"repairtsr!",
			MB_ICONWARNING | MB_OK | MB_DEFBUTTON1
		);
	}
}

spvector<baseTriangle> myModel::Model::searchEdgeTriangle() const
{
	spvector<baseTriangle> v;// = new spvector<baseTriangle>();

	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (tarr[i]->hasEdge)
		{
			v.push_back(tarr[i]);
		}
	}

	return v;
}

spvector<Segment3D> myModel::Model::searchEdge() const
{
	spvector<Segment3D> v;// = new spvector<Segment3D>();
	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (tarr[i]->hasEdge)
		{
			auto & tc = tarrinfo[i]->c;
			for (int j = 0; j < 3; j++)
			{
				char k = 1 << j;
				if (tc & k)
				{
					int inda = j * 2 % 3;
					int indb = (j * 2 + 1) % 3;
					v.push_back(make_shared<Segment3D>(tarr[i]->at(inda), tarr[i]->at(indb)));
				}
			}
		}
	}
	return v;
}

vector<myModel::Model::Edgeind> myModel::Model::searchEdgeind() const
{
	vector<Edgeind> v;// = new vector<Edgeind>();
	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		if (tarr[i]->hasEdge)
		{
			auto & tc = tarrinfo[i]->c;
			for (int j = 0; j < 3; j++)
			{
				char k = 1 << j;
				if (tc & k)
				{
					unsigned char inda = j * 2 % 3;
					unsigned char indb = (j * 2 + 1) % 3;
					v.push_back({i, inda, indb});
				}
			}
		}
	}
	return v;
}

vectors<myModel::Model::Edgeind> myModel::Model::searchEdgecontour() const
{
	// 	vector<Edgeind> sourceEdge;
	// 	auto tempv = searchEdgeind();
	// 	sourceEdge.swap(*tempv);
	// 	delete tempv;
	// 	return searchEdgecontour(sourceEdge);
	return searchEdgecontour(searchEdgeind());
}

vectors<myModel::Model::Edgeind> myModel::Model::searchEdgecontour(vector<Edgeind>& sourceEdge) const
{
	vectors<Edgeind> v;

	while (!sourceEdge.empty())
	{
		vector<Edgeind> edgec;
		vector<Edgeind> tempedge;

		tempedge.push_back(sourceEdge.back());
		sourceEdge.pop_back();

		while (!tempedge.empty())
		{
			auto currEdge = tempedge.back();
			tempedge.pop_back();
			edgec.push_back(currEdge);

			auto np1 = tarr[currEdge.a]->np[currEdge.b];
			auto np2 = tarr[currEdge.a]->np[currEdge.c];

			auto check = [&](const Edgeind& e)
			{
				auto np = tarr[e.a]->np;
				return np[e.b] == np1 || np[e.b] == np2 || np[e.c] == np1 || np[e.c] == np2;
			};

			auto it = find_if(sourceEdge.begin(), sourceEdge.end(), check);

			while (it != sourceEdge.end())
			{
				tempedge.push_back(*it);
				if (it != --sourceEdge.end())
				{
					*it = sourceEdge.back();
					sourceEdge.pop_back();

					it = find_if(it, sourceEdge.end(), check);
				}
				else
				{
					sourceEdge.pop_back();
					break;
				}
			}
		}

		v.push_back(edgec);
	}
	return v;
}

sp<baseTriangle> myModel::Model::searchTriangle(const myMath::Point3D& P) const
{
	/*Point3DNode* pn = ocd->searchNode(P);*/
	shared_ptr<ModelTriangle> t(nullptr);
// 	if (pn == nullptr)
// 	{
		double dis = DBL_MAX , tempdis;
		for (auto&& tt : tarr)
		{
			tempdis = tt->distance(P);
			if (tempdis < dis)
			{
				dis = tempdis;
				t = tt;
			}
		}		
// 	}
// 	else
// 	{
// 		double dis = DBL_MAX, tempdis;
// 		for (auto&& i : pn->adj_tri)
// 		{
// 			auto tid = i.tid();
// 
// 			tempdis = tarr[tid]->distance(P);
// 			if (tempdis < dis)
// 			{
// 				dis = tempdis;
// 				t = tarr[tid];
// 			}
// 		}
// 
// 	}
	return t;
}

myModel::tri_ind_t myModel::Model::searchTriangleind(const myMath::Point3D& P) const
{
	Point3DNode* pn = ocd->searchNode(P);

	tri_ind_t ind = -1;
	if (pn != nullptr)
	{
		double dis = DBL_MAX, tempdis;
		for (auto&& i : pn->adj_tri)
		{
			tempdis = tarr[i.tid()]->distance(P);
			if (tempdis < dis)
			{
				dis = tempdis;
				ind = i.tid();
			}
		}
	}

	if (ind == -1)
	{
		cout << "wtf searchTriangleind\n";
	}

	return ind;
}

spvector<baseTriangle> myModel::Model::getTriangles(const inddata& ind) const
{
	spvector<myMath::baseTriangle> tas;
	tas.reserve(ind.size());
	for (auto&& tind : ind)
	{
		tas.push_back(tarr[tind]);
	}
	return tas;
}

spvectors<baseTriangle> myModel::Model::getTriangles(const inddatas& inds) const
{
	spvectors<myMath::baseTriangle> tass;
	tass.reserve(inds.size());
	for (auto&& ind : inds)
	{
		spvector<myMath::baseTriangle> tas;
		tas.reserve(ind.size());
		for (auto&& tind : ind)
		{
			tas.push_back(tarr[tind]);
		}
		tass.push_back(tas);
	}
	return tass;
}

double myModel::Model::calcMeanCurvature(const Point3DNode* pn) const
{
	double aa = 0;
	for (auto& adjt : pn->adj_tri)
	{
		aa += tarr[adjt.tid()]->area();
	}
	aa /= 3;

	auto cott = [&](const Point3DNode* pn2)
	{
		auto a = find_if(pn->adj_tri.begin(), pn->adj_tri.end(), [&](const triid_pid& ind)
		{
			return tarr[ind.tid()]->np[(ind.pid() + 2) % 3] == pn2;
		});
		auto b = find_if(pn->adj_tri.begin(), pn->adj_tri.end(), [&](const triid_pid& ind)
		{
			return tarr[ind.tid()]->np[(ind.pid() + 1) % 3] == pn2;
		});

		auto& t1 = tarr[a->tid()];
		auto& t2 = tarr[b->tid()];

		auto v11 = t1->at((a->pid() + 1) % 3) - *pn->pdata;
		auto v12 = t1->at((a->pid() + 1) % 3) - *pn2->pdata;

		auto v21 = t2->at((b->pid() + 2) % 3) - *pn2->pdata;
		auto v22 = t2->at((b->pid() + 2) % 3) - *pn->pdata;

		return Vector3D::dot(v11, v12) / Vector3D::cross(v12, v11).length() + Vector3D::dot(v21, v22) / Vector3D::cross(v22, v21).length();
	};

	Vector3D vv {0,0,0};
	for (auto& adjpn : pn->adj_Point)
	{
		vv += (*adjpn->pdata - *pn->pdata)*cott(adjpn);
	}

	return signCurvature(pn) ? 0.5*(vv / aa / 2).length() : -0.5*(vv / aa / 2).length();
}

double myModel::Model::calcGaussCurvature(const Point3DNode* pn) const
{
	double aa = 0;
	double tt = 0;
	for (auto& adjt : pn->adj_tri)
	{
		auto t = tarr[adjt.tid()];
		aa += t->area();

		uVector3D v11 = t->at((adjt.pid() + 1) % 3) - *pn->pdata;
		uVector3D v12 = t->at((adjt.pid() + 2) % 3) - *pn->pdata;

		auto dd = myMath::Vector3D::dot(v11, v12);
		tt += std::acos(dd > 0 ? std::min(dd, 1.0) : std::max(dd, -1.0));
	}

	aa /= 3;

	return signCurvature(pn) ? (2 * M_PI - tt) / aa : -(2 * M_PI - tt) / aa;
}

std::pair<double, double> myModel::Model::calcCurvature(const Point3DNode* pn) const
{
	double aa = 0;
	double tt = 0;
	for (auto& adjt : pn->adj_tri)
	{
		auto t = tarr[adjt.tid()];
		aa += t->area();

		uVector3D v11 = t->at((adjt.pid() + 1) % 3) - *pn->pdata;
		uVector3D v12 = t->at((adjt.pid() + 2) % 3) - *pn->pdata;

		auto dd = myMath::Vector3D::dot(v11, v12);
		tt += std::acos(dd > 0 ? std::min(dd, 1.0) : std::max(dd, -1.0));
	}

	aa /= 3;

	auto K = (2 * M_PI - tt) / aa;

	auto cott = [&](const Point3DNode* pn2)
	{
		auto a = find_if(pn->adj_tri.begin(), pn->adj_tri.end(), [&](const triid_pid& ind)
		{
			return tarr[ind.tid()]->np[(ind.pid() + 2) % 3] == pn2;
		});
		auto b = find_if(pn->adj_tri.begin(), pn->adj_tri.end(), [&](const triid_pid& ind)
		{
			return tarr[ind.tid()]->np[(ind.pid() + 1) % 3] == pn2;
		});

		auto& t1 = tarr[a->tid()];
		auto& t2 = tarr[b->tid()];

		auto v11 = t1->at((a->pid() + 1) % 3) - *pn->pdata;
		auto v12 = t1->at((a->pid() + 1) % 3) - *pn2->pdata;

		auto v21 = t2->at((b->pid() + 2) % 3) - *pn2->pdata;
		auto v22 = t2->at((b->pid() + 2) % 3) - *pn->pdata;

		return Vector3D::dot(v11, v12) / Vector3D::cross(v12, v11).length() + Vector3D::dot(v21, v22) / Vector3D::cross(v22, v21).length();
	};


	Vector3D vv {0,0,0};
	for (auto& adjpn : pn->adj_Point)
	{
		vv += (*adjpn->pdata - *pn->pdata)*cott(adjpn);
	}

	auto H = 0.5*(vv / aa / 2).length();

	auto temp = sqrt(H*H - K);
	cout << H << " " << K << endl;
	return make_pair(H + temp, H - temp);
};

bool myModel::Model::signCurvature(const Point3DNode* pn) const
{
	Vector3D vn {0,0,0};
	Vector3D vp {0,0,0};
	for (auto& adjt : pn->adj_tri)
	{
		vn += tarr[adjt.tid()]->n;
	}
	for (auto& adjpn : pn->adj_Point)
	{
		vp += (*adjpn->pdata - *pn->pdata);
	}

	return Vector3D::dot(vn, vp) >= 0;
}

double myModel::Model::calcCurvature(const myMath::Point3D& P) const
{
	auto t = searchTriangle(P);

	if (t != nullptr)
	{
		auto mt = dynamic_cast<ModelTriangle*>(t.get());

		//return (calcMeanCurvature(mt->na) + calcMeanCurvature(mt->nb) + calcMeanCurvature(mt->nc)) / 3;
		return (calcMeanCurvature(mt->na) + calcMeanCurvature(mt->nb) + calcMeanCurvature(mt->nc)) / 3;
	}
	cout << "fuck\n";
	return 0;
}

myMath::Point3D myModel::Model::PointinModelCoor(const Point3D& P) const
{
	//auto& mdvm = ModelMatrixr;
	auto c = bbox.bboxcenter();
	// 	Vector3D vx {mdvm[0], mdvm[1], mdvm[2]};
	// 	Vector3D vy {mdvm[4], mdvm[5], mdvm[6]};
	// 	Vector3D vz {mdvm[8], mdvm[9], mdvm[10]};
	// 	Point3D pmc = P + Mtt*-1 + c*-1;
	// 	Point3D npmc;
	// 	npmc.x = Vector3D::dot(vx, pmc) + c.x;
	// 	npmc.y = Vector3D::dot(vy, pmc) + c.y;
	// 	npmc.z = Vector3D::dot(vz, pmc) + c.z;
	// 
	// 	return npmc;

	return c + tm.removeMat(P + c*-1);
}

myMath::Vector3D myModel::Model::VectorinModelCoor(const Vector3D& V) const
{
	// 	auto& mdvm = ModelMatrixr;
	// 	Vector3D vx {mdvm[0], mdvm[1], mdvm[2]};
	// 	Vector3D vy {mdvm[4], mdvm[5], mdvm[6]};
	// 	Vector3D vz {mdvm[8], mdvm[9], mdvm[10]};
	// 	Vector3D nvmc;
	// 	nvmc.x = Vector3D::dot(vx, V);
	// 	nvmc.y = Vector3D::dot(vy, V);
	// 	nvmc.z = Vector3D::dot(vz, V);
	// 
	// 	return nvmc;
	return tm.removeMat(V);
}

myMath::Point3D myModel::Model::PointinWorldCoor(const Point3D& P) const
{
	auto c = bbox.bboxcenter();

	return c + tm.applyMat(P + c*-1);
}

myMath::Vector3D myModel::Model::VectorinWorldCoor(const Vector3D& V) const
{
	return tm.applyMat(V);
}

spvector<baseTriangle> myModel::Model::shiftTriangle(double D, unsigned int dir) const
{
	spvector<baseTriangle> v;
	v.reserve(tarr.size());
	Point3D A, B, C;

	auto V = [&](const Point3DNode* const& PN)
	{
		auto& k = PN->adj_tri;
		Vector3D vv = {0.0, 0.0, 0.0};

		unsigned int s = k.size() < 3 ? 2 : 1;

		auto vec = PN->find_adjtri(s);

		for (auto&& j : vec)
		{
			vv += tarr[j]->n;
		}

		vv.normalize();

		return vv;
	};

	map<const Point3DNode* const, Vector3D> mV;
	for (auto&& i : ocd->pnarr)
	{
		mV.emplace(i, V(i));
	}

	for (auto&& i : tarr)
	{
		if (i->hasEdge)
		{
			A = *i->a + mV[i->na] * D;
			B = *i->b + mV[i->nb] * D;
			C = *i->c + mV[i->nc] * D;
		}
		else
		{
			A = *i->a + mV[i->na] * D;
			B = *i->b + mV[i->nb] * D;
			C = *i->c + mV[i->nc] * D;
		}

		if (dir == 0)
		{
			v.push_back(make_shared<Triangle>(A, C, B));
		}
		else
		{
			v.push_back(make_shared<Triangle>(A, B, C));
		}
	}
	return v;
}

spvector<baseTriangle> myModel::Model::shiftTriangle2(double D, unsigned int dir, Vector3D n) const
{
	spvector<baseTriangle> v;
	v.reserve(tarr.size());
	Point3D A, B, C;

	auto nD = n*D;

	for (auto&& i : tarr)
	{
		A = *i->a + nD;
		B = *i->b + nD;
		C = *i->c + nD;

		if (dir == 0)
		{
			v.push_back(make_shared<Triangle>(A, C, B));
		}
		else
		{
			v.push_back(make_shared<Triangle>(A, B, C));
		}
	}
	return v;
}

spvector<baseTriangle> myModel::Model::shiftTriangle3(const map<Point3DNode*, double>& mpnd, double shiftd, unsigned int dir) const
{
	spvector<baseTriangle> v;
	v.reserve(tarr.size());
	Point3D A, B, C;

	auto V = [&](const Point3DNode* const& PN)
	{
		auto& k = PN->adj_tri;
		Vector3D vv = {0.0, 0.0, 0.0};

		unsigned int s = k.size() < 3 ? 2 : 1;

		auto vec = PN->find_adjtri(s);

		for (auto&& j : vec)
		{
			vv += tarr[j]->n;
			// 			vv.x += tarr[j]->n.x;
			// 			vv.y += tarr[j]->n.y;
			// 			vv.z += tarr[j]->n.z;
		}

		vv.normalize();

		return vv;
	};

	map<const Point3DNode* const, Vector3D> mV;
	for (auto&& i : ocd->pnarr)
	{
		mV.emplace(i, V(i));
	}

	for (auto&& i : tarr)
	{
		if (i->hasEdge)
		{
			A = *i->a + mV[i->na] * (mpnd.at(i->na) + shiftd);
			B = *i->b + mV[i->nb] * (mpnd.at(i->nb) + shiftd);
			C = *i->c + mV[i->nc] * (mpnd.at(i->nc) + shiftd);
		}
		else
		{
			A = *i->a + mV[i->na] * (mpnd.at(i->na) + shiftd);
			B = *i->b + mV[i->nb] * (mpnd.at(i->nb) + shiftd);
			C = *i->c + mV[i->nc] * (mpnd.at(i->nc) + shiftd);
		}

		if (dir == 0)
		{
			v.push_back(make_shared<Triangle>(A, C, B));
		}
		else
		{
			v.push_back(make_shared<Triangle>(A, B, C));
		}
	}
	return v;
}

spvector<baseTriangle> myModel::Model::shiftTriangle4(const map<Point3DNode*, double>& mpnd, double shiftd, unsigned int dir, Vector3D n) const
{
	spvector<baseTriangle> v;
	v.reserve(tarr.size());
	Point3D A, B, C;

	for (auto&& i : tarr)
	{
		A = *i->a + n*(mpnd.at(i->na) + shiftd);
		B = *i->b + n*(mpnd.at(i->nb) + shiftd);
		C = *i->c + n*(mpnd.at(i->nc) + shiftd);

		if (dir == 0)
		{
			v.push_back(make_shared<Triangle>(A, C, B));
		}
		else
		{
			v.push_back(make_shared<Triangle>(A, B, C));
		}
	}
	return v;
}

pair<Point3D, Point3D> myModel::Model::Intersection(const Line3D& i) const
{
	auto pp = bbox.Intersection(i);
	auto ind = ABt->returnind(i);
	Point3D pm;
	Point3D pM;

	auto nn = pp.first - pp.second;
	nn.normalize();
	auto pp2 = make_pair(pp.first + nn, pp.second - nn);
	pp = pp2;

	double dm = DBL_MAX;
	double dM = DBL_MAX;
	for (auto ii : ind)
	{
		auto p = Plane(*tarr[ii]).Intersection(i);
		if (tarr[ii]->isintri(p))
		{
			auto temp1 = p.distancesq(pp.first);
			auto temp2 = p.distancesq(pp.second);
			if (temp1 < dm)
			{
				dm = temp1;
				pm = p;
			}
			if (temp2 < dM)
			{
				dM = temp2;
				pM = p;
			}
		}
	}

	return make_pair(pm, pM);
}

sp<myModel::Model> myModel::Model::filpmodel() const
{
	Model M;

	M.tm = tm;

	float* tempp = new float[tarr.size() * 9];
	for (tri_ind_t i = 0; i < tarr.size(); i++)
	{
		auto& t = *tarr[i];

		tempp[i * 9] = (float)t[0][0];
		tempp[i * 9 + 1] = (float)t[0][1];
		tempp[i * 9 + 2] = (float)t[0][2];

		tempp[i * 9 + 3] = (float)t[2][0];
		tempp[i * 9 + 4] = (float)t[2][1];
		tempp[i * 9 + 5] = (float)t[2][2];

		tempp[i * 9 + 6] = (float)t[1][0];
		tempp[i * 9 + 7] = (float)t[1][1];
		tempp[i * 9 + 8] = (float)t[1][2];
	}
	M.SetModel(tempp, tarr.size());
	delete[] tempp;

	return make_shared<Model>(move(M));
}

sp<myModel::Model> myModel::Model::refine() const
{
	return refine(AvgTriSize);
}

sp<myModel::Model> myModel::Model::refine(double edgelength) const
{
	spvector<baseTriangle> bt;

	for (auto&& t : tarr)
	{
		auto sts = SubTriangles(t.get(), edgelength).Triangulation();
		if (sts.size() > 1)
		{
			bt.insert(bt.end(), sts.begin(), sts.end());
		}
		else
		{
			bt.push_back(t);
		}		
	}

	return make_shared<Model>(bt);
}

myMath::Point3D myModel::Model::find_biggest_tri_cen_without_self_cross()
{
	auto tcr = searchTriangleCross();

	if (tcr.empty())
	{
		return tarr[0]->getcen();
	}

	vector<tri_ind_t> tind(tarr.size());
	tri_ind_t ind = 0;
	generate(tind.begin(), tind.end(), [&]() {return ind++; });

	vector<tri_ind_t> resultind;

	for (auto&& i : tind)
	{
		if (tcr.find(i) == tcr.end())
		{
			resultind.push_back(i);
		}
	}

	vector<pair<double, tri_ind_t>> area_ind(resultind.size());

	transform(resultind.begin(), resultind.end(), area_ind.begin(), [&](tri_ind_t i)
	{
		return make_pair(tarr[i]->area(), i);
	});

	auto it = max_element(area_ind.begin(), area_ind.end(), [](auto& a, auto& b)->bool {return a.first < b.first; });

	if (it->second < tarr.size())
	{
		return tarr[it->second]->getcen();
	}
	else
	{
		return tarr[0]->getcen();
	}
}
//////////////////////////////////////////////////////////////////////////
map<myModel::Point3DNode*, double> myModel::ModelSelect2::get_mpnd(double shiftd, unsigned int level) const
{
	mtime3.clear();
	mtime3.start();
	map<Point3DNode*, double> mpnd;

	vector<double> tcidvalue(mmodel->tarr.size(), shiftd);
	for (auto&& ind : tcid)
	{
		tcidvalue[ind.first] = ind.second;
	}
	
#ifndef useParallel
	for (auto&& pn : mmodel->ocd->pnarr)
#else
	concurrent_vector<pair<Point3DNode*, double>> mpndtemp;
	parallel_for_each(mmodel->ocd->pnarr.begin(), mmodel->ocd->pnarr.end(), [&](Point3DNode* pn)
#endif	
	{
		double d = 0;

		auto ind2 = pn->find_adjtri(level);

		for (auto&& ind : ind2)
		{
// 			auto it = tcid.find(ind);
// 			d += it != tcid.end() ? it->second : shiftd;

			d += tcidvalue[ind];
		}
#ifndef useParallel
		mpnd[pn] = d / ind2.size();
	}
#else
		mpndtemp.push_back(make_pair(pn, d / ind2.size()));
	});

	mpnd.insert(mpndtemp.begin(), mpndtemp.end());
#endif
	mtime3.end();
	mtime3.print("mpnd: ");
	return mpnd;
}

map<myModel::Point3DNode*, double> myModel::ModelSelect3::get_mpnd(double shiftd, unsigned int level) const
{
	mtime3.clear();
	mtime3.start();
	map<Point3DNode*, double> mpnd;

	vector<double> tcidvalue(mmodel->tarr.size(), shiftd);
	for (auto&& ind : Tr2ID)
	{
		tcidvalue[ind.first] = ind.second;
	}

#ifndef useParallel
	for (auto&& pn : mmodel->ocd->pnarr)
#else
	concurrent_vector<pair<Point3DNode*, double>> mpndtemp;
	parallel_for_each(mmodel->ocd->pnarr.begin(), mmodel->ocd->pnarr.end(), [&](Point3DNode* pn)
#endif	
	{
		double d = 0;

		auto ind2 = pn->find_adjtri(level);

		for (auto&& ind : ind2)
		{
			// 			auto it = tcid.find(ind);
			// 			d += it != tcid.end() ? it->second : shiftd;

			d += tcidvalue[ind];
		}
#ifndef useParallel
		mpnd[pn] = d / ind2.size();
	}
#else
		mpndtemp.push_back(make_pair(pn, d / ind2.size()));
	});

	mpnd.insert(mpndtemp.begin(), mpndtemp.end());
#endif
	mtime3.end();
	mtime3.print("mpnd: ");
	return mpnd;
}

map<myModel::Point3DNode*, double> myModel::ModelSelect3::get_PointLoc(double shiftd, int mode, unsigned int level) const
{
	mtime3.clear();
	mtime3.start();
	map<Point3DNode*, double> mpnd;

	vector<double> tcidvalue(mmodel->tarr.size(), shiftd);

	switch (mode) {
	case 0:
		for (auto&& ind : Tr2ID)
			tcidvalue[ind.first] = ID2Pr.at(ind.second).allowance; //allowance
			//tcidvalue[ind.first] = atid.at(ind.second).first; //allowance
		break;
	case 1:
		for (auto&& ind : Tr2ID)
			tcidvalue[ind.first] = ID2Pr.at(ind.second).allowance + ID2Pr.at(ind.second).thinkness; //thinkness
			//tcidvalue[ind.first] = atid.at(ind.second).second + atid.at(ind.second).first; //thinkness
		break;
	default:
		for (auto&& ind : Tr2ID)
			tcidvalue[ind.first] = ind.second;
	}



#ifndef useParallel
	for (auto&& pn : mmodel->ocd->pnarr)
#else
	concurrent_vector<pair<Point3DNode*, double>> mpndtemp;
	parallel_for_each(mmodel->ocd->pnarr.begin(), mmodel->ocd->pnarr.end(), [&](Point3DNode* pn)
#endif	
	{
		double d = 0;

		auto ind2 = pn->find_adjtri(level);

		for (auto&& ind : ind2)
		{
			// 			auto it = tcid.find(ind);
			// 			d += it != tcid.end() ? it->second : shiftd;

			d += tcidvalue[ind];
		}
#ifndef useParallel
		mpnd[pn] = d / ind2.size();
	}
#else
		mpndtemp.push_back(make_pair(pn, d / ind2.size()));
});

	mpnd.insert(mpndtemp.begin(), mpndtemp.end());
#endif
	mtime3.end();
	mtime3.print("mpnd: ");
	return mpnd;
}