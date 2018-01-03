#include "stdafx.h"

#include <stack>

#include "AABBtree.h"

#ifdef useParallel
#include <ppl.h>
using namespace concurrency;
#endif

#if defined (__AVX__)
#include <intrin.h>
#define __useAVX__ 1
#define mycmp2(A, B) (7 == (_mm256_movemask_pd(_mm256_cmp_pd((A), (B), _CMP_LT_OQ))))
#define mycmp3(A, B, C) (mycmp2((A), _mm256_loadu_pd((C)->M.p)) && mycmp2(_mm256_loadu_pd((C)->m.p), (B)))
#endif

using namespace std;
using namespace myMath;
using namespace myModel;

AABBtreenode::AABBtreenode() :BoundBox({{DBL_MAX, DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX, -DBL_MAX}}), child(nullptr)
{
}

AABBtreenode::~AABBtreenode()
{
	if (child)
	{
		delete[] child;
	}
}

AABBtreenode::AABBtreenode(vector<nodedata*>& d) :AABBtreenode()
{
	ndata.swap(d);
}

void AABBtreenode::Add(nodedata* d)
{
	ndata.push_back(d);
}

void AABBtreenode::calpoint()
{
	Point3D averp {0, 0, 0};
	size_t count = 1;
	double countt;
	for (auto&& i : ndata)
	{
		auto& tempm = i->m;
		auto& tempM = i->M;

		m = Point3D::min(tempm, m);
		M = Point3D::max(tempM, M);

		countt = 1.0 / count;
		averp.x += ((tempm.x + tempM.x)*0.5 - averp.x) * countt;
		averp.y += ((tempm.y + tempM.y)*0.5 - averp.y) * countt;
		averp.z += ((tempm.z + tempM.z)*0.5 - averp.z) * countt;

		count++;
	}

	cen = averp;
}

bool AABBtreenode::is_leaf() const
{
	return child == nullptr;
}

void AABBtreenode::separate(size_t lim)
{
	auto ndatasize = ndata.size();
	if (ndatasize < lim)
	{
		child = nullptr;
		return;
	}

	child = new AABBtreenode[8];

	for (int i = 0; i < 8; i++)
	{
		(child + i)->ndata.reserve(ndatasize / 8);
	}

	for (auto&& i : ndata)
	{
		auto ff = f0(*i, cen);
		auto ff1 = ff.first;
		
		if (ff1)
		{
			auto ff2 = ff.second | ff1;
			for (int j = 0; j < 8; j++)
			{
				if ((j | ff1) == ff2)
				{
					(child + j)->Add(i);
				}
			}
		}
		else
		{
			(child + ff.second)->Add(i);
		}
	}
	
	for (int j = 0; j < 8; j++)
	{
		if (ndatasize == (child + j)->ndata.size())
		{
			delete[] child;
			child = nullptr;
			break;
		}
	}
}

template<typename boolean_allind>
vector<nodedata*> AABBtreenode::returnind(const nodedata& d, const boolean_allind& allind) const
{
	vector<AABBtreenode*> nodestack;

	vector<nodedata*> t;
	nodestack.push_back(const_cast<AABBtreenode*>(this));

	while (!nodestack.empty())
	{
		auto b = nodestack.back();
		nodestack.pop_back();

		if (b->iscross(d))
		{
			if (b->is_leaf())
			{
				pushind(t, d, b->ndata, allind);
			}
			else
			{
				auto flagf = f0(d, b->cen);
 				auto f2 = flagf.first | flagf.second;
				for (int i = 0; i < 8; i++)
				{
 					if ((i | flagf.first) == f2)
 					{
						nodestack.push_back(b->child + i);
					}
				}
			}
		}
	}
	
	sort(t.begin(), t.end());
	t.erase(unique(t.begin(), t.end()), t.end());
	t.shrink_to_fit();

	return t;
}

template<typename boolean_allind>
void AABBtreenode::pushind(vector<nodedata*>& t, const nodedata& d, const vector<nodedata*>& nd, const boolean_allind& allind)
{
	vector<nodedata*>::const_iterator it;
	if (allind)
	{
		it = nd.begin();
	}
	else
	{
		it = upper_bound(nd.begin(), nd.end(), &d);
	}

#if !defined (__useAVX__)
	for (; it != nd.end(); ++it)
	{
		if (d.iscross(**it))
		{
			t.push_back(*it);
		}
	}
#else
	__m256d dm = _mm256_loadu_pd(d.m.p);
	dm.m256d_f64[3] = NAN;
	__m256d dM = _mm256_loadu_pd(d.M.p);
	dM.m256d_f64[3] = NAN;

	for (; it != nd.end(); ++it)
	{
		if (mycmp3(dm, dM, *it))
		{
			t.push_back((*it));
		}
	}
#endif
}

AABBtreenode::flag_pair AABBtreenode::f0(const nodedata& d, const Point3D& cen)
{	
	return make_pair((d.m < cen)&(cen < d.M), d.M < cen);
}

AABBtree::~AABBtree()
{
	if (root)
	{
		delete root;
	}
}

AABBtree::AABBtree(vector<nodedata>& s) : root(nullptr)
{
	sourcendata.swap(s);

	vector<nodedata*> d(sourcendata.size());
	for (size_t i = 0; i < sourcendata.size(); i++)
	{
		d[i] = sourcendata.data() + i;
	}

	root = new AABBtreenode(d);

	Buildtree();
}

AABBtree::AABBtree(vector<nodedata>&& s) : root(nullptr), sourcendata(s)
{
	vector<nodedata*> d(sourcendata.size());
	for (size_t i = 0; i < sourcendata.size(); i++)
	{
		d[i] = sourcendata.data() + i;
	}

	root = new AABBtreenode(d);

	Buildtree();
}

void AABBtree::Buildtree()
{
#ifndef useParallel
	stack<AABBtreenode*> temp;
	temp.push(root);

	while (!temp.empty())
	{
		auto& i = *temp.top();
		temp.pop();

		i.calpoint();

		i.separate(512);

		if (i.is_leaf())
		{
			
		}
		else
		{
			for (unsigned int j = 0; j < 8; j++)
			{
				temp.push(i.child + j);
			}
		}
	}
#else
	std::function<void(AABBtreenode*)> cal_separate = [&](AABBtreenode* i)->void
	{
		i->calpoint();//114,104

		i->separate(512);//116,129

		if (i->is_leaf())
		{

		}
		else
		{
			parallel_for(0, 8, [&](int j)
			//for (int j = 0; j < 8; j++)
			{
				cal_separate(i->child + j);
			});
		}
	};
	mtime3.clear();
	mtime3.start();
	cal_separate(root);//230,235//89,93
	mtime3.end();
	mtime3.print("cal2: ");
#endif
}

void AABBtree::reBuildtree()
{
	vector<nodedata*> d;

	if (root)
	{
		d.swap(root->ndata);
		delete root;
	}
	else
	{
		for (size_t i = 0; i < sourcendata.size(); i++)
		{
			d[i] = sourcendata.data() + i;
		}
	}

	root = new AABBtreenode(d);

	Buildtree();
}

inddata transform_nodedata(const vector<nodedata*>& a, const nodedata* b)
{
// 	inddata result(a.size());
// 	transform(a.begin(), a.end(), result.begin(), [&](nodedata* t)
// 	{
// 		return t - b;
// 	});

	inddata result(a.size());
	auto it = result.begin();
	for (auto&& t : a)
	{
		*it++ = t - b;
	}

	return result;
}

inddata AABBtree::returnind(const BoundBox& bbox) const
{
	return transform_nodedata(root->returnind(bbox), sourcendata.data());
}

inddatas AABBtree::returnind() const
{
	return returnind(*this);
}

inddatas AABBtree::returnind(const AABBtree& rhs) const
{
	if (this == &rhs)
	{
		return returnind_impl(rhs, constant<bool, false>());
	}
	else
	{
		return returnind_impl(rhs, constant<bool, true>());
	}
}

inddata AABBtree::returnind(const Line3D& rhs) const
{
	auto pp = root->Intersection(rhs);

	auto nn = pp.first - pp.second;
	nn.normalize();
	auto pp2 = make_pair(pp.first + nn, pp.second - nn);
	pp = pp2;

	auto bb = BoundBox({Point3D::min(pp.first, pp.second), Point3D::max(pp.first, pp.second)});
	auto abt = AABBtree(vector<nodedata>(1, bb));
	auto inds = abt.returnind(*this);

	return inds[0];
}

template<typename boolean_allind>
inddatas AABBtree::returnind_impl(const AABBtree& rhs, const boolean_allind& allind) const
{
	auto& not_same_model = allind;
	inddatas resultd;
	size_t rns = root->ndata.size();
	resultd.resize(rns);

#ifndef useParallel
#pragma omp parallel for shared(resultd)
	for (size_t iii = 0; iii < rns; iii++)//for (auto&& i : *root->ndata)
	{
		auto ii = root->ndata[iii];
#else
	parallel_for_each(root->ndata.begin(), root->ndata.end(), [&](const nodedata* ii)
	{	
#endif
		auto& i = *ii;

		AABBtreenode* currentnode = rhs.root;

		while (!not_same_model || currentnode->iscross(i))
		{
			auto checkd = AABBtreenode::f0(i, currentnode->cen);
			
			if (currentnode->is_leaf() || checkd.first)
			{
				resultd[ii - root->ndata[0]] = transform_nodedata(currentnode->returnind(i, allind), rhs.root->ndata[0]);
				break;
			}
			currentnode = currentnode->child + checkd.second;
		}
	}
#ifdef useParallel
	);
#endif

	return resultd;
}