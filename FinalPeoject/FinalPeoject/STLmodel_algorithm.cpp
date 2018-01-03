#include "stdafx.h"
#include <future>
#include "STLmodel.h"
#include "cgal.h"

using namespace std;
using namespace myMath;
//////////////////////////////////////////////////////////////////////////
sp<myModel::Model> myModel::ModelCutter::cutModel(const Model& M, double(&mdvm)[16])
{
	size_t ps = sparr.size();
	if (ps < 3)
	{
		return nullptr;
	}
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
	Plane planenn({mdvm[2], mdvm[6], mdvm[10]}, transp[0]);
	vector<Segment3D> seg;
	for (size_t i = 0; i < sparr.size(); i++)
	{
		seg.emplace_back(transp[i], transp[i + 1]);
	}

	auto isinpolygon = [&](const Point3D& P)->bool
	{
		Point3D testp = planenn.Projection(P);

		double sumangle = 0;
		for (auto&& s : seg)
		{
			Vector3D aaa = s.pstart - testp;
			Vector3D bbb = s.pend - testp;
			aaa.normalize();
			bbb.normalize();

			double d1 = Vector3D::dot(aaa, bbb);
			d1 = d1 > 0 ? std::min(d1, 1.0) : std::max(d1, -1.0);
			sumangle += std::copysign(acos(d1), Vector3D::dot(Vector3D::cross(aaa, bbb), planenn.n));
		}

		return abs(sumangle / 
			/ 2) > 0.99;
	};

	spvector<baseTriangle> tarr;
	for (auto&& t : M.tarr)
	{
		if (isinpolygon(t->getcen()) && isinpolygon(*t->a) && isinpolygon(*t->b) && isinpolygon(*t->c))
		{
			tarr.push_back(t);
		}
	}
	if (tarr.size() == 0)
	{
		cout << "\nshit\n";
		return nullptr;
	}

	return make_shared<Model>(tarr);
}
//////////////////////////////////////////////////////////////////////////
myModel::ModelManger::ModelManger(const string& filename, bool compucross) :ModelManger(filename, make_shared<Model>(filename), compucross)
{
}

myModel::ModelManger::ModelManger(const string& filename, shared_ptr<Model> Mmodel, bool compucross) : mmodel(Mmodel), modelfilefullname(filename), isshown(true), isMatchange(false), transparent_value(10)
{
	if (!mmodel->tarr.empty())
	{
		auto s1 = modelfilefullname.find_last_of('\\') + 1;
		auto slen = modelfilefullname.length() - s1;
		modelfilepath = modelfilefullname.substr(0, s1);
		modelname = modelfilefullname.substr(s1, slen);
		auto p = modelname.find_last_of('.');
		modelnamewtoext = (p == string::npos ? modelname : modelname.substr(0, p));

		renderflags.reset();
		renderflags[5] = true;
		
		for (auto&& i : renderarr)
		{
			i = nullptr;
		}
		initialthis(compucross);
	}
}

myModel::ModelManger::~ModelManger()
{
	for (auto&& i : renderarr)
	{
		if (i)
		{
			delete i;
		}
	}
}

void myModel::ModelManger::initialthis(bool compucross)
{
	auto secc = std::async([this]() {return mmodel->searchEdgecontour_fordisplay(); });
	auto tr1c = std::async([this]() {return mmodel->searchEdgeTriangle(); });
	auto tr2c = std::async([this]() {return mmodel->searchDegenerationTriangle(3e-3); });

	mtime2.clear();
	mtime2.start();
	auto tcr = compucross ? mmodel->searchTriangleCross() : Model::TriCrossResult();
	mtime2.end();
	mtime2.print("AABB: ");

	mtime2.clear();
	mtime2.start();
	auto tsr = compucross ? mmodel->separateTriangle(tcr) : Model::TriSeparateResult();
	mtime2.end();
	mtime2.print("ST: ");

	auto sec = secc.get();
	auto tr1 = tr1c.get();
	auto tr2 = tr2c.get();

	visibleflags[0] = !sec.empty();
	visibleflags[1] = !tr1.empty();
	visibleflags[2] = !tr2.empty();
	visibleflags[3] = true;
	visibleflags[4] = !tcr.empty();
	visibleflags[5] = true;
	visibleflags[6] = !tsr.empty();
	visibleflags[7] = true;
	visibleflags[8] = true;
	//visibleflags[9] = true;

	mtime2.clear();
	mtime2.start();
	renderarr[0] = new OBJrender::Segment3Dss_Render(move(sec));
	renderarr[1] = new OBJrender::Triangles_Render(move(tr1), mycolor4 {0, 0, 1, 1});
	renderarr[2] = new OBJrender::Triangles_Render(move(tr2), mycolor4 {1, 0, 0, 1});
	renderarr[3] = new OBJrender::BoundBox_Render(mmodel->bbox);
	renderarr[4] = new OBJrender::Segment3Dss_Render(Model::searchTriangleCross_fordisplay2(tcr));
	renderarr[5] = new OBJrender::Model_Render(mmodel);
	renderarr[6] = new OBJrender::Triangles_Render(Model::subtri_fordisplay(tsr), mycolor4 {0, 0, 1, 1});
	renderarr[7] = new OBJrender::OctCloud_Render(mmodel->ocd, mycolor4 {0, 1, 0, 1});
	renderarr[8] = new OBJrender::PCAAxis_Render(mmodel->PAxis, mycolor4{ 0, 1, 0, 1 });
	//renderarr[9] = new OBJrender::LDI_Render(*(mmodel->ldi), mycolor4{ 0, 1, 0, 1 });

	mtime2.end();
	mtime2.print("Render: ");

	cout << "tcr size: " << tcr.size() << endl;
}

void myModel::ModelManger::updateMat()
{
	if (isMatchange)
	{
		mmodel->ApplyTMat();
		reinitial();
		isMatchange = false;
	}
}

void myModel::ModelManger::modifyname(const string& name)
{	
	modelname = name;
	auto p = modelname.find_last_of('.');
	modelnamewtoext = (p == string::npos ? modelname : modelname.substr(0, p));
	modelfilefullname = modelfilepath + modelname;
}

void myModel::ModelManger::reinitial(bool compucross)
{
	for (auto&& i : renderarr)
	{
		if (i)
		{
			delete i;
		}
	}
	initialthis(compucross);
}

void myModel::ModelManger::updatevbo(const inddata& ind)
{
	if (ind.empty())
	{
		return;
	}

	for (auto&& i : ind)
	{
		mmodel->tarr[i]->caln();
	}

	auto it = remove_if(mmodel->tarr.begin(), mmodel->tarr.end(), [](sp<ModelTriangle>& t)
	{
		return *t->a == *t->b || *t->a == *t->c || *t->b == *t->c;
	});

	if (it != mmodel->tarr.end())
	{
		cout << "Remove triangle: " << distance(it, mmodel->tarr.end()) << endl;
		mmodel->tarr.erase(it, mmodel->tarr.end());
		
		std::thread th1([this]()
		{
			mmodel->RebuildOctCloud();
		});
		std::thread th2([this]()
		{
			mmodel->BuildAABB();
		});
		mmodel->PAxis.CalModelPCA(&(mmodel->tarr),&(mmodel->bbox));
		
		delete renderarr[5];
		renderarr[5] = new OBJrender::Model_Render(mmodel);
		
		th1.join();
		th2.join();
	}
	else
	{	
		std::thread th1([this]()
		{
			mmodel->ocd->resegmentation();
		});//23~25ms
		std::thread th2([this, &ind]()
		{
			for (auto&& i : ind)
			{
				mmodel->ABt->sourcendata[i] = BoundBox::genbbox(*mmodel->tarr[i]);
			}

			mmodel->ABt->reBuildtree();
		});//18ms

		mtime3.clear();
		mtime3.start();
// 		delete renderarr[5];
// 		renderarr[5] = new OBJrender::Model_Render(mmodel);
		mmodel->modifyvbo(ind);
		mtime3.end();
		mtime3.print("vbo: ");
		
		th1.join();
		th2.join();
	}
}

void myModel::ModelManger::releaseVBO()
{
	for (auto&& r : renderarr)
	{
		r->releaseVBO();
	}
}
//////////////////////////////////////////////////////////////////////////
void myModel::ColorHistogram::fontall(float x, float y) const
{
	OBJrender::displayfunc<std::string> strdf;
  	std::stringstream ss;
 	std::string str;
 	ss << center - halfrange;
 	ss >> str;
 
 	strdf(str, {0, 0, 0, 1}, {x, y - 15, 0});
 
 	ss.clear();
 	str.clear();
 	ss << center;
 	ss >> str;
 	strdf(str, {0, 0, 0, 1}, {x, y + 135, 0});
 
 	ss.clear();
 	str.clear();
 	ss << center + halfrange;
 	ss >> str;
 	strdf(str, {0, 0, 0, 1}, {x, y + 285, 0});
}
//////////////////////////////////////////////////////////////////////////
void myModel::ColorHistogram2::fontall(float x, float y) const
{
	OBJrender::displayfunc<std::string> strdf;
	auto str = my_to_string(lowlim, 2);
  	strdf(str, {0, 0, 0, 1}, {x - str.size() * 9, y, 0});
	str = my_to_string((hilim + lowlim) / 2, 2);
  	strdf(str, {0, 0, 0, 1}, {x - str.size() * 9, y + 150, 0});
	str = my_to_string(hilim, 2);
  	strdf(str, {0, 0, 0, 1}, {x - str.size() * 9, y + 300, 0});
}

void myModel::ColorHistogram2::fontcurrt(float x, float y) const
{
	OBJrender::displayfunc<std::string> strdf;
   
  	strdf(my_to_string(getcurrtd(), 2), {0, 0, 0, 1}, {x, y + currtc, 0});
}
//////////////////////////////////////////////////////////////////////////
myModel::ColoredModel::ColoredModel(const shared_ptr<Model>& ma, const shared_ptr<Model>& mb, const shared_ptr<ColorHistogram2>& ct) :VBOobj(1), mmodel(ma), mmodel2(mb), colorht(ct)
{
}

myModel::ColoredModel2::ColoredModel2(const shared_ptr<Model>& m, const shared_ptr<ColorHistogram2>& ct) : VBOobj(1), mmodel(m), colorht(ct)
{
}
//////////////////////////////////////////////////////////////////////////
myMath::Point3D myModel::RegionGrowing::tricen(tri_ind_t tind) const
{
	auto& pp = mmodel->tarr.at(tind)->p;
	return (*pp[0] + *pp[1] + *pp[2]) / 3;
};

double myModel::RegionGrowing::tricur(tri_ind_t tind1, tri_ind_t tind2, vertex_ind_t pind) const
{
	auto& n2 = mmodel->tarr.at(tind2)->n;

	auto& t = mmodel->tarr.at(tind1);
	auto& p = t->p;
	Point3D* pa = p[(pind + 1) % 3];
	Point3D* pb = p[pind];

	double nn = (Vector3D::dot(t->n, n2) - 1) / 2;

	return copysign(nn, Vector3D::dot(n2, *pa - *pb));
}

void myModel::RegionGrowing::addseed(const Point3D& p, bool newseed)
{
	if (growstep > 0)
	{
		growstep = 0;
		groupcount = 0;
		submodelind.clear();
		renderind.clear();
	}
	tri_ind_t ind = mmodel->searchTriangleind(mmodel->PointinModelCoor(p));

	if (ind > mmodel->tarr.size())
	{
		return;
	}

	if (find_if(seed0.begin(), seed0.end(), [ind](const pair<groupind, tri_ind_t>& p) {return p.second == ind; }) != seed0.end())
	{
		return;
	}

	if (newseed && seed0.size() != 0)
	{
		groupcount++;
	}

	if (newseed || seed0.size() == 0)
	{
		renderind.push_back(vector<gl_tri_ind_t>());
	}

	seed0.emplace_back(groupcount, ind);

	renderind.back().push_back(gl_tri_ind_t(ind * 3 + 0));
	renderind.back().push_back(gl_tri_ind_t(ind * 3 + 1));
	renderind.back().push_back(gl_tri_ind_t(ind * 3 + 2));
}

void myModel::RegionGrowing::Growing1()
{
	growstep = 1;

	vector<double> seeddisarr;
	vector<Point3D> tricenarr;
	seeddisarr.resize(seed0.size(), DBL_MAX);
	tricenarr.resize(seed0.size());
	transform(seed0.begin(), seed0.end(), tricenarr.begin(), [&](const pair<groupind, tri_ind_t>& p) {return tricen(p.second); });

	vector<seedtemp1> seed1;
	for (size_t i = 0; i < seed0.size(); i++)
	{
		auto& p1 = tricenarr[i];
		for (size_t j = i + 1; j < seed0.size(); j++)
		{
			if (seed0[i].first == seed0[j].first)
			{
				continue;
			}
			seeddisarr[i] = min(seeddisarr[i], p1.distancesq(tricenarr[j]));
			seeddisarr[j] = min(seeddisarr[j], seeddisarr[i]);
		}

		seed1.push_back({seed0[i].second, seed0[i].first, p1, 0.1225*(seeddisarr[i])});
	}

	seed0.clear();
	vector<multiset<seedtemp2, greater<seedtemp2>>> seed2;
	seed2.resize(groupcount + 1);

	while (!seed1.empty())
	{
		auto currseed = seed1.back();
		seed1.pop_back();

		if (!growcheck[currseed.tind])
		{
			submodelind[currseed.group].push_back(currseed.tind);
			growcheck[currseed.tind] = true;

			auto& tt = mmodel->tarrinfo[currseed.tind]->b;
			for (auto&& t : tt)
			{
				if (!growcheck[t.tid()])
				{
					if (currseed.cen.distancesq(tricen(t.tid())) < currseed.dis)
					{
						seed1.push_back({t.tid(), currseed.group, currseed.cen, currseed.dis});
					}
					else
					{
						seed2[currseed.group].insert({tricur(currseed.tind, t.tid(), t.pid()), currseed.group, t.tid()});
					}
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	size_t seed2count = 0;
	vector<seedtemp2> maxcur;
	for (auto&& group : seed2)
	{
		seed2count += group.size();
		maxcur.push_back(*group.begin());
	}
	while (seed2count != 0)
	{
		auto& currseed = *max_element(maxcur.begin(), maxcur.end());
		auto& currgroup = seed2[currseed.group];
		if (currgroup.empty())
		{
			continue;
		}
		currgroup.erase(currgroup.begin());
		seed2count--;

		if (!growcheck[currseed.tind])
		{
			submodelind[currseed.group].push_back(currseed.tind);
			growcheck[currseed.tind] = true;
			auto& tt = mmodel->tarrinfo[currseed.tind]->b;

			for (auto&& t : tt)
			{
				if (!growcheck[t.tid()])
				{
					currgroup.insert({tricur(currseed.tind, t.tid(), t.pid()), currseed.group, t.tid()});
					seed2count++;
				}
			}
		}

		if (currgroup.empty())
		{
			maxcur[currseed.group].cur = -DBL_MAX;
		}
		else
		{
			maxcur[currseed.group] = *currgroup.begin();
		}
	}
}

void myModel::RegionGrowing::Growing()
{
	growcheck = new bool[mmodel->tarr.size()];
	fill(growcheck, growcheck + mmodel->tarr.size(), false);

	submodelind.resize(groupcount + 1);

	Growing1();

	delete[] growcheck;

	renderind.clear();
	renderind.resize(groupcount + 1);
	for (size_t i = 0; i < renderind.size(); i++)
	{
		for (auto&& j : submodelind[i])
		{
			renderind[i].push_back(gl_tri_ind_t(j * 3 + 0));
			renderind[i].push_back(gl_tri_ind_t(j * 3 + 1));
			renderind[i].push_back(gl_tri_ind_t(j * 3 + 2));
		}
	}
}

spvector<myModel::Model> myModel::RegionGrowing::resultmodel()
{
	spvector<Model> result;
	for (auto&& sm : submodelind)
	{
		spvector<baseTriangle> temp;
		temp.reserve(sm.size());
		for (auto&& ind : sm)
		{
			temp.push_back(mmodel->tarr[ind]);
		}
		result.push_back(std::make_shared<Model>(std::move(temp)));
	}
	return result;
}
//////////////////////////////////////////////////////////////////////////
void myModel::VaildRegionGrowing::addsubtri(shared_ptr<baseTriangle>& t)
{
	auto checkit = [&](const eihv& e)
	{
		return e.t == t && e.b != true;
	};

	auto it2 = find_if(earr.begin(), earr.end(), checkit);

	if (it2 != earr.end())
	{
		subtri.push_back(t);
	}
	while (it2 != earr.end())
	{
		estack.push_back(*it2);
		it2 = find_if(++it2, earr.end(), checkit);
	}
};

void myModel::VaildRegionGrowing::subgrow(triid_pid t, tri_ind_t currseed)
{
	growcheck[t.tid()] = 1;
	auto& currt = mmodel->tarr[currseed];

	size_t currthv[3] = {hp(currt->at(0)), hp(currt->at(1)), hp(currt->at(2))};
	size_t hh = currthv[0] ^ currthv[1] ^ currthv[2] ^ currthv[t.pid()];
	float lensq = (float)(currt->at((t.pid() + 1) % 3) - currt->at((t.pid() + 2) % 3)).lengthsq();
	estack.push_back({hh, lensq, currt, currseed, t.pid(), false});

	while (!estack.empty())
	{
		auto tempe = estack.back();
		estack.pop_back();

		if (tempe.b)
		{
			continue;
		}

		auto it = equal_range(earr.begin(), earr.end(), tempe);
		auto itl = it.second - it.first;

		for (auto iii = it.first; iii != it.second; iii++)
		{
			iii->b = true;
		}

		if (itl == 0)
		{
			cout << "wtf vg meet hole\n";
		}
		else if (itl == 1)
		{
			auto& iii = it.first;
			if (iii->t != tempe.t)
			{
				addsubtri(iii->t);
			}
			else
			{
				auto& exitt = mmodel->tarr[tempe.tid];

				size_t thv[3] = {hp(exitt->at(0)), hp(exitt->at(1)), hp(exitt->at(2))};
				size_t exite[3] = {thv[1] ^ thv[2], thv[0] ^ thv[2], thv[0] ^ thv[1]};

				vertex_ind_t ee = -1;
				for (vertex_ind_t i = 0; i < 3; i++)
				{
					if (tempe.hv == exite[i])
					{
						ee = i;
						break;
					}
				}
				if (ee > 2)
				{
					cout << "wtf exit?\n";
					continue;
				}

				auto& adjt = mmodel->tarrinfo[tempe.tid]->b;
				auto eit = find_if(adjt.begin(), adjt.end(), [&](const triid_pid& ti) {return ti.pid() == ee; });
				if (eit != adjt.end() && growcheck[eit->tid()] == 0)
				{
					growstack.push_back(eit->tid());
					//cout << "exit!!\n";
				}
			}
		}
		else if (itl == 2)
		{
			for (auto iii = it.first; iii != it.second; iii++)
			{
				if (iii->t != tempe.t)
				{
					addsubtri(iii->t);
					break;
				}
			}
		}
		else if (itl >= 4)
		{
			Point3D p1 = tempe.t->at(tempe.eid);
			Point3D p2 = tempe.t->at((tempe.eid + 1) % 3);

			uVector3D p2p1 = p2 - p1;
			double p1p2p = -2;
			auto it4 = it.first;
			for (auto iii = it.first; iii != it.second; iii++)
			{
				if (iii->tid == tempe.tid)
				{
					continue;
				}

				auto& p = iii->t->at(iii->eid);
				auto v1 = p - p1;

				if (signbit(Vector3D::dot(v1, tempe.t->n)) != signbit(Vector3D::dot(v1, iii->t->n)))
				{
					uVector3D p2p = p2 - p;
					double p1p2ptemp = Vector3D::dot(p2p, p2p1);
					if (p1p2ptemp > p1p2p)
					{
						p1p2p = p1p2ptemp;
						it4 = iii;
					}
				}
			}

			if (p1p2p != -2)
			{
				growcheck[it4->tid] = 1;
				addsubtri(it4->t);
			}
			else
			{
				cout << __FUNCTION__ << __LINE__ << " " << itl << endl;
			}
		}
		else
		{
			cout << setprecision(20);
			for (auto iii = it.first; iii != it.second; iii++)
			{
				cout << iii->hv << " " << iii->lensq << " " << iii->t << " " << iii->tid << " " << iii->eid << endl;
				for (edge_ind_t ijk = 0; ijk < 3; ijk++)
				{
					if (ijk == iii->eid)
					{
						continue;
					}
					cout << iii->t->at(ijk).x << " " << iii->t->at(ijk).y << " " << iii->t->at(ijk).z << endl;
				}
			}
			cout << itl << " www\n";

			Plane p1(*it.first->t);
			Plane p2(*(it.first + 1)->t);
			Plane p3(*(it.first + 2)->t);

			if (it.first->tid != (it.first + 1)->tid)
			{
				auto& a = it.first->t;
				auto& b = (it.first + 1)->t;

				int aa = p2.onplane(a->at(0)) + p2.onplane(a->at(1)) + p2.onplane(a->at(2));
				int bb = p1.onplane(b->at(0)) + p1.onplane(b->at(1)) + p1.onplane(b->at(2));
				cout << aa << bb << endl;
			}
			else if (it.first->tid != (it.first + 2)->tid)
			{
				auto& a = it.first->t;
				auto& b = (it.first + 2)->t;

				int aa = p3.onplane(a->at(0)) + p3.onplane(a->at(1)) + p3.onplane(a->at(2));
				int bb = p1.onplane(b->at(0)) + p1.onplane(b->at(1)) + p1.onplane(b->at(2));
				cout << aa << bb << endl;
			}
			else if ((it.first + 2)->tid != (it.first + 1)->tid)
			{
				auto& a = (it.first + 2)->t;
				auto& b = (it.first + 1)->t;

				int aa = p2.onplane(a->at(0)) + p2.onplane(a->at(1)) + p2.onplane(a->at(2));
				int bb = p3.onplane(b->at(0)) + p3.onplane(b->at(1)) + p3.onplane(b->at(2));
				cout << aa << bb << endl;
			}
		}
	}
}

void myModel::VaildRegionGrowing::setseedtri(const Point3D& p)
{
	if (growstep > 0)
	{
		growstep = 0;
		renderind.clear();
	}

	tri_ind_t ind = mmodel->searchTriangleind(mmodel->PointinModelCoor(p));

	if (ind > mmodel->tarr.size())
	{
		return;
	}
	seedtri.push_back(ind);
	renderind.push_back(gl_tri_ind_t(ind * 3 + 0));
	renderind.push_back(gl_tri_ind_t(ind * 3 + 1));
	renderind.push_back(gl_tri_ind_t(ind * 3 + 2));
}

void myModel::VaildRegionGrowing::Growing()
{
	growstep = 1;

	//////////////////////////////////////////////////////////////////////////
	tcr = mmodel->searchTriangleCross();
	tsr = mmodel->separateTriangle(tcr);
	//////////////////////////////////////////////////////////////////////////

	growcheck = new int[mmodel->tarr.size()];
	int ind = 0;
	generate(growcheck, growcheck + mmodel->tarr.size(), [&]() {return (tsr.find(ind++) == tsr.end()) - 1; });

	for (auto&& i : tsr)
	{
		for (auto&& j : i.second)
		{
			size_t hh[3] = {hp(j->at(0)), hp(j->at(1)), hp(j->at(2))};
			earr.push_back({hh[0] ^ hh[1], (float)(j->at(0) - j->at(1)).lengthsq(), j, i.first, 2, false});
			earr.push_back({hh[0] ^ hh[2], (float)(j->at(0) - j->at(2)).lengthsq(), j, i.first, 1, false});
			earr.push_back({hh[2] ^ hh[1], (float)(j->at(2) - j->at(1)).lengthsq(), j, i.first, 0, false});
		}
	}
	sort(earr.begin(), earr.end());

	growstack.swap(seedtri);

	while (!growstack.empty())
	{
		auto currseed = growstack.back();
		growstack.pop_back();

		if (growcheck[currseed] == 0)
		{
			resultind.push_back(currseed);
			growcheck[currseed] = 1;

			auto& tt = mmodel->tarrinfo[currseed]->b;
			for (auto&& t : tt)
			{
				if (growcheck[t.tid()] == 0)
				{
					growstack.push_back(t.tid());
				}
				else if (growcheck[t.tid()] == -1)
				{
					subgrow(t, currseed);
				}
			}
		}
	}

	delete[] growcheck;
	tcr.clear();
	tsr.clear();

	sort(resultind.begin(), resultind.end());
	renderind.clear();
	renderind.reserve(resultind.size() * 3);
	for (auto&& i : resultind)
	{
		renderind.push_back(gl_tri_ind_t(i * 3 + 0));
		renderind.push_back(gl_tri_ind_t(i * 3 + 1));
		renderind.push_back(gl_tri_ind_t(i * 3 + 2));
	}
}

sp<myModel::Model> myModel::VaildRegionGrowing::resultmodel()
{
	spvector<baseTriangle> subtri2;
	subtri2.reserve(mmodel->tarr.size() + subtri.size());
	subtri2 = subtri;
	for (auto&& i : resultind)
	{
		subtri2.push_back(mmodel->tarr[i]);
	}

	return make_shared<Model>(subtri2);
}
//////////////////////////////////////////////////////////////////////////
void myModel::ModelDigHole::addCylinder(const Cylinder& cy_)
{
	cyarr.emplace_back(cy_);
}

sp<myModel::Model> myModel::ModelDigHole::Triangulation(bool twoside) const
{
	spvector<baseTriangle> resultarr;

	resultarr.reserve(Cylinder::nedge * 4 * cyarr.size());
	for (auto&& i : cyarr)
	{
		auto cy = i.cy.Triangulation(twoside);
		resultarr.insert(resultarr.end(), cy.begin(), cy.end());
	}

	return make_shared<Model>(resultarr);
}

sp<myModel::Model> myModel::ModelDigHole::Triangulation2(bool twoside) const
{
	spvector<baseTriangle> resultarr;

	std::vector<std::pair<double, double>> tocgal;
	for (auto& p : brush->pm.poly.parr)
	{
		tocgal.push_back(make_pair(p.x, p.y));
	}

	auto fromcgal = testtest3(tocgal);

	spvector<baseTriangle> face;
	auto h1 = twoside ? -1.0 : 0.0;
	auto h2 = 1.0;
	for (size_t i = 0; i < fromcgal.size(); i+= 6)
	{
		Point3D a = {fromcgal[i], fromcgal[i + 1], h2};
		Point3D b = {fromcgal[i + 2], fromcgal[i + 3], h2};
		Point3D c = {fromcgal[i + 4], fromcgal[i + 5], h2};
		Point3D d = {fromcgal[i], fromcgal[i + 1], h1};
		Point3D f = {fromcgal[i + 2], fromcgal[i + 3], h1};
		Point3D e = {fromcgal[i + 4], fromcgal[i + 5], h1};

		face.push_back(make_shared<myMath::Triangle>(a, b, c));
		face.push_back(make_shared<myMath::Triangle>(d, e, f));
	}

	auto n = brush->pm.poly.parr.size();
	auto& parr = brush->pm.poly.parr;
	for (size_t i = 0; i < n; i++)
	{
		auto ind = (i + 1) % n;
		Point3D a = {parr[i].x, parr[i].y, h2};
		Point3D b = {parr[ind].x, parr[ind].y, h2};
		Point3D c = {parr[ind].x, parr[ind].y, h1};
		Point3D d = {parr[i].x, parr[i].y, h1};

		face.push_back(make_shared<myMath::Triangle>(a, c, b));
		face.push_back(make_shared<myMath::Triangle>(c, a, d));
	}

	for (auto&& i : cyarr)
	{
		auto tm2 = i.tm;
		tm2.t = i.cy.cen;

		TransformationMatrix tm;
		auto zx = i.cy.n;
		uVector3D xx = Vector3D::genVVec(zx);
		uVector3D yx = Vector3D::cross(zx, xx);

		tm.r[0] = xx.x*i.cy.radius;
		tm.r[1] = xx.y*i.cy.radius;
		tm.r[2] = xx.z*i.cy.radius;
		tm.r[4] = yx.x*i.cy.radius;
		tm.r[5] = yx.y*i.cy.radius;
		tm.r[6] = yx.z*i.cy.radius;
		tm.r[8] = zx.x*i.cy.height;
		tm.r[9] = zx.y*i.cy.height;
		tm.r[10] = zx.z*i.cy.height;
		

		for (auto& t : face)
		{
			resultarr.push_back(make_shared<myMath::Triangle>(tm2.applyMat(tm.applyMat(*t))));
		}
	}

	return make_shared<Model>(resultarr);
}

sp<myModel::Model> myModel::ModelDigHole::Triangulation3() const
{
	spvector<baseTriangle> resultarr;

	std::vector<std::pair<double, double>> tocgal;
	for (auto& p : brush->pm.poly.parr)
	{
		tocgal.push_back(make_pair(p.x, p.y));
	}

	auto fromcgal = testtest3(tocgal);

	spvector<baseTriangle> face;
	auto h1 = -1.0;
	auto h2 = 1.0;
	for (size_t i = 0; i < fromcgal.size(); i += 6)
	{
		Point3D a = {fromcgal[i], fromcgal[i + 1], h2};
		Point3D b = {fromcgal[i + 2], fromcgal[i + 3], h2};
		Point3D c = {fromcgal[i + 4], fromcgal[i + 5], h2};

		face.push_back(make_shared<myMath::Triangle>(a, b, c));
	}

	auto n = brush->pm.poly.parr.size();
	auto& parr = brush->pm.poly.parr;
	for (size_t i = 0; i < n; i++)
	{
		auto ind = (i + 1) % n;
		Point3D a = {parr[i].x, parr[i].y, h2};
		Point3D b = {parr[ind].x, parr[ind].y, h2};
		Point3D c = {0, 0, h1};

		face.push_back(make_shared<myMath::Triangle>(a, c, b));
	}

	for (auto&& i : cyarr)
	{
		auto tm2 = i.tm;
		tm2.t = i.cy.cen;

		TransformationMatrix tm;
		auto zx = i.cy.n;
		uVector3D xx = Vector3D::genVVec(zx);
		uVector3D yx = Vector3D::cross(zx, xx);

		auto r = 2 * (i.cy.radius + 1);

		tm.r[0] = xx.x*r;
		tm.r[1] = xx.y*r;
		tm.r[2] = xx.z*r;
		tm.r[4] = yx.x*r;
		tm.r[5] = yx.y*r;
		tm.r[6] = yx.z*r;
		tm.r[8] = zx.x*(i.cy.radius + 1);
		tm.r[9] = zx.y*(i.cy.radius + 1);
		tm.r[10] = zx.z*(i.cy.radius + 1);

		for (auto& t : face)
		{
			resultarr.push_back(make_shared<myMath::Triangle>(tm2.applyMat(tm.applyMat(*t))));
		}
	}

	return make_shared<Model>(resultarr);
}

sp<myModel::Model> myModel::ModelDigHole::Triangulation4() const
{
	spvector<baseTriangle> resultarr;

	std::vector<std::pair<double, double>> tocgal;
	for (auto& p : brush->pm.poly.parr)
	{
		tocgal.push_back(make_pair(p.x, p.y));
	}

	auto fromcgal = testtest3(tocgal);

	spvector<baseTriangle> faceu;
	spvector<baseTriangle> facel;
	auto h1 = -1.0;
	auto h2 = 1.0;
	for (size_t i = 0; i < fromcgal.size(); i += 6)
	{
		Point3D a = {fromcgal[i], fromcgal[i + 1], h2};
		Point3D b = {fromcgal[i + 2], fromcgal[i + 3], h2};
		Point3D c = {fromcgal[i + 4], fromcgal[i + 5], h2};
		Point3D d = {fromcgal[i], fromcgal[i + 1], h1};
		Point3D f = {fromcgal[i + 2], fromcgal[i + 3], h1};
		Point3D e = {fromcgal[i + 4], fromcgal[i + 5], h1};

		faceu.push_back(make_shared<myMath::Triangle>(a, b, c));
		facel.push_back(make_shared<myMath::Triangle>(d, e, f));
	}

	auto n = brush->pm.poly.parr.size();
	auto& parr = brush->pm.poly.parr;
	
	double r_ratio = 1.25;
	double r_shift = 0;
	double r_hratio = 0.25;
	for (auto&& i : cyarr)
	{
		auto tm2 = i.tm;
		tm2.t = i.cy.cen;

		TransformationMatrix tm, tm3, tmm;
		auto zx = i.cy.n;
		uVector3D xx = Vector3D::genVVec(zx);
		uVector3D yx = Vector3D::cross(zx, xx);

		tm.r[0] = xx.x*i.cy.radius;
		tm.r[1] = xx.y*i.cy.radius;
		tm.r[2] = xx.z*i.cy.radius;
		tm.r[4] = yx.x*i.cy.radius;
		tm.r[5] = yx.y*i.cy.radius;
		tm.r[6] = yx.z*i.cy.radius;
		tm.r[8] = zx.x*i.cy.height;
		tm.r[9] = zx.y*i.cy.height;
		tm.r[10] = zx.z*i.cy.height;

		tm3.r[0] = xx.x * (i.cy.radius*r_ratio + r_shift + i.cy.height*r_hratio);
		tm3.r[1] = xx.y * (i.cy.radius*r_ratio + r_shift + i.cy.height*r_hratio);
		tm3.r[2] = xx.z * (i.cy.radius*r_ratio + r_shift + i.cy.height*r_hratio);
		tm3.r[4] = yx.x * (i.cy.radius*r_ratio + r_shift + i.cy.height*r_hratio);
		tm3.r[5] = yx.y * (i.cy.radius*r_ratio + r_shift + i.cy.height*r_hratio);
		tm3.r[6] = yx.z * (i.cy.radius*r_ratio + r_shift + i.cy.height*r_hratio);
		tm3.r[8] = zx.x * i.cy.height*r_hratio;
		tm3.r[9] = zx.y * i.cy.height*r_hratio;
		tm3.r[10] = zx.z * i.cy.height*r_hratio;

		tmm.r[0] = xx.x*i.cy.radius;
		tmm.r[1] = xx.y*i.cy.radius;
		tmm.r[2] = xx.z*i.cy.radius;
		tmm.r[4] = yx.x*i.cy.radius;
		tmm.r[5] = yx.y*i.cy.radius;
		tmm.r[6] = yx.z*i.cy.radius;
		tmm.r[8] = zx.x;
		tmm.r[9] = zx.y;
		tmm.r[10] = zx.z;

		for (auto& t : facel)
		{
			resultarr.push_back(make_shared<myMath::Triangle>(tm2.applyMat(tm.applyMat(*t))));
		}

		for (auto& t : faceu)
		{
			resultarr.push_back(make_shared<myMath::Triangle>(tm2.applyMat(tm3.applyMat(*t))));
		}

		for (size_t ii = 0; ii < n; ii++)
		{
			auto ind = (ii + 1) % n;
			Point3D a = tm2.applyMat(tm3.applyMat(Point3D {parr[ii].x, parr[ii].y, h2}));
			Point3D b = tm2.applyMat(tm3.applyMat(Point3D {parr[ind].x, parr[ind].y, h2}));
			
			Point3D c = tm2.applyMat(tmm.applyMat(Point3D {parr[ii].x, parr[ii].y, -(i.cy.radius*(r_ratio - 1) + r_shift)}));
			Point3D d = tm2.applyMat(tmm.applyMat(Point3D {parr[ind].x, parr[ind].y, -(i.cy.radius*(r_ratio - 1) + r_shift)}));

			Point3D e = tm2.applyMat(tm.applyMat(Point3D {parr[ii].x, parr[ii].y, h1}));
			Point3D f = tm2.applyMat(tm.applyMat(Point3D {parr[ind].x, parr[ind].y, h1}));


			resultarr.push_back(make_shared<myMath::Triangle>(a, d, b));
			resultarr.push_back(make_shared<myMath::Triangle>(d, a, c));
			resultarr.push_back(make_shared<myMath::Triangle>(c, f, d));
			resultarr.push_back(make_shared<myMath::Triangle>(f, c, e));
		}
	}

	return make_shared<Model>(resultarr);
}

sp<myModel::Model> myModel::ModelDigHole::Triangulation5() const
{
	spvector<baseTriangle> resultarr;

	std::vector<std::pair<double, double>> tocgal;
	for (auto& p : brush->pm.poly.parr)
	{
		tocgal.push_back(make_pair(p.x, p.y));
	}

	auto fromcgal = testtest3(tocgal);

	spvector<baseTriangle> faceu;
	spvector<baseTriangle> facel;
	auto h1 = -1.0;
	auto h2 = 1.0;
	for (size_t i = 0; i < fromcgal.size(); i += 6)
	{
		Point3D a = {fromcgal[i], fromcgal[i + 1], h2};
		Point3D b = {fromcgal[i + 2], fromcgal[i + 3], h2};
		Point3D c = {fromcgal[i + 4], fromcgal[i + 5], h2};
		Point3D d = {fromcgal[i], fromcgal[i + 1], h1};
		Point3D f = {fromcgal[i + 2], fromcgal[i + 3], h1};
		Point3D e = {fromcgal[i + 4], fromcgal[i + 5], h1};

		faceu.push_back(make_shared<myMath::Triangle>(a, b, c));
		facel.push_back(make_shared<myMath::Triangle>(d, e, f));
	}

	auto n = brush->pm.poly.parr.size();
	auto& parr = brush->pm.poly.parr;
	
	double r_ratio = 1.25;
	double r_shift = 0;
	double r_hratio = 0.25;
	for (auto&& i : cyarr)
	{
		auto tm2 = i.tm;
		tm2.t = i.cy.cen;

		TransformationMatrix tm, tm3, tmm;
		auto zx = i.cy.n;
		uVector3D xx = Vector3D::genVVec(zx);
		uVector3D yx = Vector3D::cross(zx, xx);

		double rrr = i.cy.radius*r_ratio + r_shift;
		double h_shift = (1 - sqrt(1 - rrr*rrr*i.k*i.k)) / i.k;

		if (isnan(h_shift))
		{
			cout << "radius too big." << endl;
			continue;
		}

		tm.r[0] = xx.x*i.cy.radius;
		tm.r[1] = xx.y*i.cy.radius;
		tm.r[2] = xx.z*i.cy.radius;
		tm.r[4] = yx.x*i.cy.radius;
		tm.r[5] = yx.y*i.cy.radius;
		tm.r[6] = yx.z*i.cy.radius;
		tm.r[8] = zx.x*i.cy.height;
		tm.r[9] = zx.y*i.cy.height;
		tm.r[10] = zx.z*i.cy.height;		

		tm3.r[0] = xx.x * (rrr + i.cy.height*r_hratio);
		tm3.r[1] = xx.y * (rrr + i.cy.height*r_hratio);
		tm3.r[2] = xx.z * (rrr + i.cy.height*r_hratio);
		tm3.r[4] = yx.x * (rrr + i.cy.height*r_hratio);
		tm3.r[5] = yx.y * (rrr + i.cy.height*r_hratio);
		tm3.r[6] = yx.z * (rrr + i.cy.height*r_hratio);
		tm3.r[8] = zx.x * (i.cy.height*r_hratio + h_shift);
		tm3.r[9] = zx.y * (i.cy.height*r_hratio + h_shift);
		tm3.r[10] = zx.z * (i.cy.height*r_hratio + h_shift);

		tmm.r[0] = xx.x*i.cy.radius;
		tmm.r[1] = xx.y*i.cy.radius;
		tmm.r[2] = xx.z*i.cy.radius;
		tmm.r[4] = yx.x*i.cy.radius;
		tmm.r[5] = yx.y*i.cy.radius;
		tmm.r[6] = yx.z*i.cy.radius;
		tmm.r[8] = zx.x;
		tmm.r[9] = zx.y;
		tmm.r[10] = zx.z;

		for (auto& t : facel)
		{
			resultarr.push_back(make_shared<myMath::Triangle>(tm2.applyMat(tm.applyMat(*t))));
		}

		for (auto& t : faceu)
		{
			resultarr.push_back(make_shared<myMath::Triangle>(tm2.applyMat(tm3.applyMat(*t))));
		}

		for (size_t ii = 0; ii < n; ii++)
		{
			auto ind = (ii + 1) % n;
			Point3D a = tm2.applyMat(tm3.applyMat(Point3D {parr[ii].x, parr[ii].y, h2}));
			Point3D b = tm2.applyMat(tm3.applyMat(Point3D {parr[ind].x, parr[ind].y, h2}));

			Point3D c = tm2.applyMat(tmm.applyMat(Point3D {parr[ii].x, parr[ii].y, -(rrr - i.cy.radius) + h_shift}));
			Point3D d = tm2.applyMat(tmm.applyMat(Point3D {parr[ind].x, parr[ind].y, -(rrr - i.cy.radius) + h_shift}));

			Point3D e = tm2.applyMat(tm.applyMat(Point3D {parr[ii].x, parr[ii].y, h1}));
			Point3D f = tm2.applyMat(tm.applyMat(Point3D {parr[ind].x, parr[ind].y, h1}));


			resultarr.push_back(make_shared<myMath::Triangle>(a, d, b));
			resultarr.push_back(make_shared<myMath::Triangle>(d, a, c));
			resultarr.push_back(make_shared<myMath::Triangle>(c, f, d));
			resultarr.push_back(make_shared<myMath::Triangle>(f, c, e));
		}
	}

	return make_shared<Model>(resultarr);
}
//////////////////////////////////////////////////////////////////////////
spvector<baseTriangle> myModel::ModelHoleFill::getresult(const spvector<Segment3D>& ss)
{
	spvector<baseTriangle> triarr;

	vector<Point3D> pp;
	for (auto&& s: ss)
	{
		pp.push_back(s->pstart);
		pp.push_back(s->pend);
	}
	
	bool revflag = true;

	for (size_t i = 0; i < pp.size()-2; i+=2)
	{
		if (pp[i] == pp[i+2])
		{
			swap(pp[i], pp[i + 1]);
			if (i == 0)
			{
				revflag = false;
			}
		}
		else if (pp[i] == pp[i+3])
		{
			swap(pp[i], pp[i + 1]);
			swap(pp[i + 2], pp[i + 3]);
			if (i == 0)
			{
				revflag = false;
			}
		}
		else if (pp[i + 1] == pp[i + 3])
		{
			swap(pp[i + 2], pp[i + 3]);
		}
	}

	pp.erase(unique(pp.begin(), pp.end()), pp.end());

	if (!(pp[0] == pp.back()))
	{
		cout << "?????????\n";

		return triarr;
	}

	vector<tuple<double, double, double>> sss;
	for (auto&& p:pp)
	{
		sss.emplace_back(p.x, p.y, p.z);
	}

	auto dtri_ind = holefilling(sss);

	for (size_t i = 0; i < dtri_ind.size(); i += 9)
	{
		Point3D aa = {dtri_ind[i], dtri_ind[i + 1], dtri_ind[i + 2]};
		Point3D bb = {dtri_ind[i + 3], dtri_ind[i + 4], dtri_ind[i + 5]};
		Point3D cc = {dtri_ind[i + 6], dtri_ind[i + 7], dtri_ind[i + 8]};
		auto temptri = revflag ? make_shared<Triangle>(aa, cc, bb) : make_shared<Triangle>(aa, bb, cc);
		triarr.push_back(move(temptri));
	}

	return triarr;
}

spvector<baseTriangle> myModel::ModelHoleFill::getresult(const sp<Model> mmodel, const vector<Model::Edgeind>& ss)
{
	spvector<baseTriangle> triarr;

	vector<Point3DNode*> pp;
	for (auto&& s : ss)
	{
		pp.push_back(mmodel->tarr[s.a]->np[s.b]);
		pp.push_back(mmodel->tarr[s.a]->np[s.c]);
	}

	bool revflag = true;

	for (size_t i = 0; i < pp.size() - 2; i += 2)
	{
		if (pp[i] == pp[i + 2])
		{
			swap(pp[i], pp[i + 1]);
			if (i == 0)
			{
				revflag = false;
			}
		}
		else if (pp[i] == pp[i + 3])
		{
			swap(pp[i], pp[i + 1]);
			swap(pp[i + 2], pp[i + 3]);
			if (i == 0)
			{
				revflag = false;
			}
		}
		else if (pp[i + 1] == pp[i + 3])
		{
			swap(pp[i + 2], pp[i + 3]);
		}
	}

	pp.erase(unique(pp.begin(), pp.end()), pp.end());

	if (!(pp[0] == pp.back()))
	{
		cout << "?????????\n";

		return triarr;
	}

	vector<tuple<double, double, double>> sss;
	for (auto&& p : pp)
	{
		sss.emplace_back(p->pdata->x, p->pdata->y, p->pdata->z);
	}

	auto dtri_ind = holefilling(sss);

	
	for (size_t i = 0; i < dtri_ind.size(); i += 9)
	{
		Point3D aa = {dtri_ind[i], dtri_ind[i + 1], dtri_ind[i + 2]};
		Point3D bb = {dtri_ind[i + 3], dtri_ind[i + 4], dtri_ind[i + 5]};
		Point3D cc = {dtri_ind[i + 6], dtri_ind[i + 7], dtri_ind[i + 8]};
		auto temptri = revflag ? make_shared<Triangle>(aa, cc, bb) : make_shared<Triangle>(aa, bb, cc);
		triarr.push_back(move(temptri));
	}

	return triarr;
}

spvector<baseTriangle> myModel::ModelHoleFill::getresult(const sp<Model> mmodel, const vectors<Model::Edgeind>& sss)
{
	spvector<baseTriangle> triarr;
	for (auto&& ss : sss)
	{
		auto t = getresult(mmodel, ss);
		triarr.insert(triarr.begin(), t.begin(), t.end());
	}
	return triarr;
}

spvector<baseTriangle> myModel::ModelHoleFill::getresult(const vector<Point3D>& pp)
{
	spvector<baseTriangle> triarr;

// 	if (!(pp[0] == pp.back()))
// 	{
// 		cout << "?????????\n";
// 
// 		return triarr;
// 	}

	vector<tuple<double, double, double>> sss;
	for (auto&& p : pp)
	{
		sss.emplace_back(p.x, p.y, p.z);
	}

	auto dtri_ind = holefilling(sss);

	for (size_t i = 0; i < dtri_ind.size(); i += 9)
	{
		Point3D aa = {dtri_ind[i], dtri_ind[i + 1], dtri_ind[i + 2]};
		Point3D bb = {dtri_ind[i + 3], dtri_ind[i + 4], dtri_ind[i + 5]};
		Point3D cc = {dtri_ind[i + 6], dtri_ind[i + 7], dtri_ind[i + 8]};
		auto temptri = make_shared<Triangle>(aa, bb, cc);
		triarr.push_back(move(temptri));
	}

	return triarr;
}
//////////////////////////////////////////////////////////////////////////
void myModel::RegionGrowingControl::setmodel(sp<Model>& m)
{
	mmodel = m;
	update();
}

void myModel::RegionGrowingControl::update()
{
	rg.reset();
	rgr.reset();

	if (!seed.empty())
	{
		rg = make_shared<RegionGrowing>(mmodel);
		rgr = make_shared<OBJrender::Growing_Render>(rg);
				
		for (auto&& ss : seed)
		{
			bool newflag = true;
			for (auto&& s : ss)
			{
				rg->addseed(s, newflag);
				newflag = false;
			}
		}
		
		rg->Growing();
	}	
}

void myModel::RegionGrowingControl::addseed(const Point3D& p, bool newflag)
{
	if (seed.empty() || newflag)
	{
		seed.push_back(std::vector<myMath::Point3D>());
	}
	seed.back().push_back(p);	

	update();
}