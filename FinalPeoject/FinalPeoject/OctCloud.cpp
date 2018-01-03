#include "stdafx.h"
#include <algorithm>
#include <bitset>
#include <limits>
#include <set>
#include <stack>
#include <unordered_map>
#include "OctCloud.h"

#ifdef useParallel
#include <ppl.h>
#include <concurrent_vector.h>
using namespace concurrency;
#endif

#if defined (__AVX__)
#include <intrin.h>
#endif

using namespace std;
using namespace myModel;

Point3DNode::Point3DNode(myMath::Point3D* Data, const triid_pid& tpid) :pdata(Data), triind_(tpid), pp(*Data)
{
}

inddata Point3DNode::find_adjtri(size_t turns) const
{	
	inddata tind;
	set<Point3DNode*> pn_set;
	set<Point3DNode*> pn1;
	set<Point3DNode*> pn2;
	
	pn1.insert(const_cast<Point3DNode*>(this));
	while (turns > 0 && !pn1.empty())
	{
		for (auto&& p : pn1)
		{
			for (auto&& ind : p->adj_tri)
			{
				tind.push_back(ind.tid());
			}
		}
		if (--turns)
		{
			pn_set.insert(pn1.begin(), pn1.end());
			for (auto&& p : pn1)
			{
				pn2.insert(p->adj_Point.begin(), p->adj_Point.end());
			}
			pn1.clear();
			set_difference(pn2.begin(), pn2.end(), pn_set.begin(), pn_set.end(), inserter(pn1, pn1.begin()));
		}
	}
	

	sort(tind.begin(), tind.end());
	tind.erase(unique(tind.begin(), tind.end()), tind.end());
	return tind;
}

inddatas Point3DNode::find_adjtris(size_t turns) const
{
	inddatas tinds;
	set<tri_ind_t> tindset;
	set<Point3DNode*> pn_set;
	set<Point3DNode*> pn1;
	set<Point3DNode*> pn2;

	pn1.insert(const_cast<Point3DNode*>(this));
	while (turns > 0 && !pn1.empty())
	{
		inddata tind;
		inddata tindtemp;
		for (auto&& p : pn1)
		{
			for (auto&& ind : p->adj_tri)
			{
				auto ii = ind.tid();
				if (tindset.insert(ii).second)
				{
					tind.push_back(ii);
				}				
			}
		}
		sort(tind.begin(), tind.end());
		unique_copy(tind.begin(), tind.end(), inserter(tindtemp, tindtemp.begin()));
		tinds.push_back(move(tindtemp));

		if (--turns)
		{
			pn_set.insert(pn1.begin(), pn1.end());
			for (auto&& p : pn1)
			{
				pn2.insert(p->adj_Point.begin(), p->adj_Point.end());
			}
			pn1.clear();
			set_difference(pn2.begin(), pn2.end(), pn_set.begin(), pn_set.end(), inserter(pn1, pn1.begin()));
		}		
	}
	
	return tinds;
}
//////////////////////////////////////////////////////////////////////////
OctCloudNode::OctCloudNode() :BoundBox({{DBL_MAX, DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX, -DBL_MAX}}), Nodeptr(nullptr)
{
}

OctCloudNode::~OctCloudNode()
{
	if (Nodeptr)
	{
		delete[] Nodeptr;
	}
}

OctCloudNode::OctCloudNode(vector<Point3DNode*>& d):OctCloudNode()
{
	p.swap(d);
}

void OctCloudNode::Add(Point3DNode* d)
{
	p.push_back(d);
}

void OctCloudNode::calpoint()
{
	Point3D averp {0, 0, 0};
	size_t count = 1;
 	double countt;
 	for (auto&& i : p)
 	{
		auto& tempp = i->pp;//*(i->pdata);//

		m = Point3D::min(tempp, m);
		M = Point3D::max(tempp, M);

		countt = 1.0 / count;
		averp.x += (tempp.x - averp.x) * countt;
		averp.y += (tempp.y - averp.y) * countt;
		averp.z += (tempp.z - averp.z) * countt;
 
 		count++;
 	}

	cen = averp;

	/*combinable<Point3D> ppaverp([&]() {return averp; });
	combinable<Point3D> mpp([&]() { return m; });
	combinable<Point3D> Mpp([&](){ return M; });

	parallel_for_each(p.begin(), p.end(), [&](Point3DNode* i)
	{
		myMath::Point3D& tempp = *(i->pdata);

		mpp.local() = Point3D::min(tempp, mpp.local());
		Mpp.local() = Point3D::max(tempp, Mpp.local());

		ppaverp.local().x += tempp.x;
		ppaverp.local().y += tempp.y;
		ppaverp.local().z += tempp.z;
	});

	m = mpp.combine([](const Point3D& a, const Point3D& b) {return Point3D::min(a, b); });
	M = Mpp.combine([](const Point3D& a, const Point3D& b) {return Point3D::max(a, b); });

	cen = ppaverp.combine([](const Point3D& a, const Point3D& b) {return a + b; })/p.size();*/
}

void OctCloudNode::calpoint2()
{
//  	__m256d mp = _mm256_set1_pd(DBL_MAX);
//  	__m256d Mp = _mm256_set1_pd(-DBL_MAX);
//  	__m256d cp = _mm256_setzero_pd();
//  	unsigned int ps = p.size();
//  
//   	combinable<__m256d> mpp([&](){ return mp; });
//   	combinable<__m256d> Mpp([&](){ return Mp; });
//   	combinable<__m256d> cpp([&](){ return cp; });
//   
//    	parallel_for(0u, ps, [&](unsigned int i)
//    	{
//    		__m256d tp = _mm256_loadu_pd(p.at(i)->data->p);
//    
//    		mpp.local() = _mm256_min_pd(mpp.local(), tp);
//    		Mpp.local() = _mm256_max_pd(Mpp.local(), tp);
//    
//    		cpp.local() = _mm256_add_pd(cpp.local(), tp);
//    	});
//    
//    	mp = mpp.combine([](const __m256d& a, const __m256d& b){return _mm256_min_pd(a, b); });
//    	Mp = Mpp.combine([](const __m256d& a, const __m256d& b){return _mm256_max_pd(a, b); });
//    	cp = cpp.combine([](const __m256d& a, const __m256d& b){return _mm256_add_pd(a, b); });
//  
// //  	for (unsigned int i = 0; i < ps; i++)
// //  	{
// // 		__m256d tp = _mm256_loadu_pd(p.at(i)->data->p);
// //  
// //  		mp = _mm256_min_pd(mp, tp);
// //  		Mp = _mm256_max_pd(Mp, tp);
// //  
// //  		cp = _mm256_add_pd(cp, tp);
// //  	}
//  
//  	for (unsigned int i = 0; i < 3; i++)
//  	{
//  		cen[i] = cp.m256d_f64[i] / ps;
//  		maxp[i] = Mp.m256d_f64[i];
//  		minp[i] = mp.m256d_f64[i];
//  	}
}

bool OctCloudNode::is_leaf() const
{
	return Nodeptr == nullptr;
}

void OctCloudNode::separate(size_t lim)
{	
	auto psize = p.size();
	if (psize < lim)
	{
		Nodeptr = nullptr;
		return;
	}

	Nodeptr = new OctCloudNode[8];

	for (int i = 0; i < 8; i++)
	{
		(Nodeptr + i)->p.reserve(psize / 8);
	}


	for (auto&& i : p)
	{
		(Nodeptr + (/**(i->pdata)*/i->pp < cen))->Add(i);

	}

	int count = 0;
	for (int i = 0; i < 8 && count < 2; i++)
	{
		count += (Nodeptr + i)->p.empty();
	}
	if (count == 2)
	{
		delete[] Nodeptr;
		Nodeptr = nullptr;
	}
	 
// 	for (int i = 0; i < 8; i++)
// 	{
// 		if ((Nodeptr + i)->p.empty())
// 		{
// 			delete[] Nodeptr;
// 			Nodeptr = nullptr;
// 			break;
// 		}
// 	}
}

bool OctCloudNode::isinbox(const myMath::Point3D& P, double lim) const
{
	myMath::Vector3D VM = P-M;
	myMath::Vector3D Vm = m-P;

	return (VM.x < lim && VM.y < lim && VM.z < lim) && (Vm.x < lim && Vm.y < lim && Vm.z < lim);
}

void OctCloudNode::remove_duplicate_point()
{
	unordered_map<myMath::Point3D, Point3DNode*> mmap(p.size());

	for (auto&& i : p)
	{
		auto& ii = mmap.emplace(*i->pdata, i);
		if (!ii.second)
		{
			ii.first->second->adj_tri.push_back(i->triind_);
			delete i;
			i = nullptr;
		}
		else
		{
			i->adj_tri.push_back(i->triind_);
		}
	}
	p.erase(remove_if(p.begin(), p.end(), [&](Point3DNode*& pp){return pp == nullptr; }), p.end());
}

void OctCloudNode::remove_duplicate_point3()
{
	hash<myMath::Point3D> hp;
	typedef pair<size_t, reference_wrapper<Point3DNode*>> phv;
	vector<phv> phvarr;
	phvarr.reserve(p.size());

	for (auto&& pit : p)
	{
		phvarr.emplace_back(hp(*pit->pdata), ref(pit));
	}
	
	sort(phvarr.begin(), phvarr.end(), [&](const phv& p1, const phv& p2){return p1.first < p2.first; });
	//////////////////////////////////////////////////////////////////////////
	auto checknullptr = [](const phv& p){return p.second.get() != nullptr; };
	auto addadj = [](Point3DNode*& p1, Point3DNode*& p2)
	{
		if (p1->pdata->equal_tof(*p2->pdata))
		{
			p1->adj_tri.push_back(p2->triind_);
			delete p2;
			p2 = nullptr;
			return false;
		}
		return true;
	};

	auto sit = phvarr.begin();
	while (sit != phvarr.end())
	{
		bool collision = false;
		auto eit = sit+1;
		while  (eit != phvarr.end() && eit->first == sit->first)
		{
			if (addadj(sit->second.get(), eit->second.get()))
			{
				collision = true;
			}
			eit++;
		}
		if (collision)
		{
			cout << "collision~~~\n";
			auto nit = find_if(sit + 1, eit, checknullptr);
			cout << sit->second.get()->pdata->at(0) << " " << sit->second.get()->pdata->at(1) << " " << sit->second.get()->pdata->at(2) << endl;
			cout << nit->second.get()->pdata->at(0) << " " << nit->second.get()->pdata->at(1) << " " << nit->second.get()->pdata->at(2) << endl;
			while (nit != eit)
			{
				auto nit2 = find_if(nit + 1, eit, checknullptr);
				while (nit2 != eit)
				{
					addadj(nit->second.get(), nit2->second.get());
					nit2 = find_if(nit2 + 1, eit, checknullptr);
				}
				nit = find_if(nit + 1, eit, checknullptr);
			}
		}
		sit = eit;
	}

	//////////////////////////////////////////////////////////////////////////
	p.erase(remove_if(p.begin(), p.end(), [](Point3DNode* p1){return p1 == nullptr; }), p.end());

	for_each(p.begin(), p.end(), [&](Point3DNode*& p1)
	{
		p1->adj_tri.push_back(p1->triind_);
		sort(p1->adj_tri.begin(), p1->adj_tri.end());
		
		auto it = adjacent_find(p1->adj_tri.begin(), p1->adj_tri.end());
		while (it != p1->adj_tri.end())
		{
			cout << it->tid() << " " << it->pid() << " ";
			it++;
			cout << it->tid() << " " << it->pid() << "\n";
			it = adjacent_find(it, p1->adj_tri.end());
		}
	});

	p.shrink_to_fit();
}
//////////////////////////////////////////////////////////////////////////
OctCloud::OctCloud(vector<Point3DNode*>& PN) :hasduplicate(true)
{
	root = new OctCloudNode(PN);

	segmentation();
}

OctCloud::~OctCloud()
{
	if (root)
	{
		delete root;
	}

	for (auto&& p : pnarr)
	{
		delete p;
	}
}

Point3DNode* OctCloud::searchNode(const myMath::Point3D& P, double lim)
{
	if (root == nullptr || LeafNode.empty() || lim < 0)
	{
		return nullptr;
	}

	Point3DNode* nearest = nullptr;
	double dis = lim*lim, tempdis;

	if (LeafNode.size() > 27)
	{
		vector<pair<double, OctCloudNode*>> ocnarr;
		vector<pair<double, OctCloudNode*>*> ocnparr;
		ocnarr.reserve(LeafNode.size());
		ocnparr.reserve(LeafNode.size());

		auto pin = ocnarr.data();
		for (auto&& i : LeafNode)
		{
			auto tdis = i->distancesq(P);
			if (tdis < dis)
			{
				ocnarr.emplace_back(tdis, i);
				ocnparr.push_back(pin++);
			}
		}

		sort(ocnparr.begin(), ocnparr.end(), [](pair<double, OctCloudNode*>* p1, pair<double, OctCloudNode*>* p2){return p1->first < p2->first; });

		size_t ocnsize = 27 < ocnparr.size() ? 27 : ocnparr.size();
		for (size_t i = 0; i < ocnsize; i++)
		{
			for (auto&& j : ocnparr[i]->second->p)
			{
				tempdis = P.distancesq(*(j->pdata));
				if (tempdis < dis)
				{
					nearest = j;
					dis = tempdis;
				}
			}
		}
	}
	else
	{
		for (auto&& i : LeafNode)
		{
			for (auto&& j : i->p)
			{
				tempdis = P.distancesq(*(j->pdata));
				if (tempdis < dis)
				{
					nearest = j;
					dis = tempdis;
				}
			}
		}
	}

	return nearest;
}

void OctCloud::segmentation()
{
	int ps = (int)sqrtf((float)root->p.size()) + 1;

	const unsigned int pointsize = max(512, (int)pow(2, ceil(log2(ps))));//(512 > ps ? 512 : ps);//  /** 6*/;

	LeafNode.clear();
	
	mtime3.clear();
	mtime3.start();
#ifndef useParallel
	stack<OctCloudNode*> nodestack;

	nodestack.push(root);
	while (nodestack.size() > 0)
	{
		OctCloudNode* ocn = nodestack.top();
		nodestack.pop();


		ocn->calpoint();
		ocn->separate(pointsize);

		if (ocn->is_leaf())
		{
			if (ocn->p.size() > 0)
			{
				LeafNode.push_back(ocn);
			}
		}
		else
		{
			ocn->p.swap(std::vector<myModel::Point3DNode*>());

			for (int i = 0; i < 8; i++)
			{
				if ((ocn->Nodeptr + i)->p.size() > 0)
				{
					nodestack.push(ocn->Nodeptr + i);
				}
			}
		}
	}

#else
	concurrent_vector<OctCloudNode*> conLeafNode;

	std::function<void(OctCloudNode*)> cal_separate = [&](OctCloudNode* ocn)->void
	{
		ocn->calpoint();
		//300,417
		//274,281
		ocn->separate(pointsize);
		//200,364
		//196,250
		if (ocn->is_leaf())
		{
			if (!ocn->p.empty())
			{
				conLeafNode.push_back(ocn);
			}
		}
		else
		{
			ocn->p.swap(std::vector<myModel::Point3DNode*>());

			parallel_for(0, 8, [&](int i)
			//for (int i = 0; i < 8; i++)
			{
				if (!(ocn->Nodeptr + i)->p.empty())
				{
					cal_separate(ocn->Nodeptr + i);
				}
			});
		}
	};

	cal_separate(root);
	//524,790//243,382
	//478,543//227,320

	LeafNode.resize(conLeafNode.size());
	copy(conLeafNode.begin(), conLeafNode.end(), LeafNode.begin());


#endif
	mtime3.end();
	mtime3.print("cal: ");
	//////////////////////////////////////////////////////////////////////////
	if (hasduplicate)
	{
		mtime3.clear();
		mtime3.start();
		size_t s = 0;
#ifndef useParallel
		for (auto&& i : LeafNode)
		{
			i->remove_duplicate_point3();
			s += i->p.size();
		}
#else
		combinable<size_t> ss;
		parallel_for_each(LeafNode.begin(), LeafNode.end(), [&](OctCloudNode* i)
		{
			i->remove_duplicate_point3();
			ss.local() += i->p.size();
		});
		s = ss.combine(plus<size_t>());
#endif
		mtime3.end();
		mtime3.print("duplicate: ");

		pnarr.reserve(s);

		for (auto&& i : LeafNode)
		{
			pnarr.insert(pnarr.end(), i->p.begin(), i->p.end());
// 			for (auto&& j : i->p)
// 			{
// 				pnarr.push_back(j);
// 			}
		}

		hasduplicate = false;
	}	
}

void OctCloud::resegmentation()
{
	if (root == nullptr || pnarr.size() == 0)
	{
		return;
	}

	delete root;
	root = new OctCloudNode();
	root->p = pnarr;
	
	for (auto&& pn:pnarr)
	{
		pn->pp = *pn->pdata;
	}

	segmentation();
}