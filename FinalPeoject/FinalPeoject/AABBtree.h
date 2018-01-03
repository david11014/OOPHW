#pragma once

#ifndef __AABBtree_H__
#define __AABBtree_H__

#include "BoundBox.h"
#include "myModelbase.h"
#include "VBOobj.h"

namespace myModel
{
	using namespace std;
	using namespace myMath;

	typedef BoundBox nodedata;
	
	template <class T, T V>
	struct constant
	{
		operator T() const { return V; }
	};
	
	class AABBtreenode: public BoundBox
	{
	public:
		typedef pair<unsigned int, unsigned int> flag_pair;

		AABBtreenode* child;
		vector<nodedata*> ndata;
		
		Point3D cen;

		AABBtreenode();
		AABBtreenode(const AABBtreenode&) = delete;
		~AABBtreenode();
		AABBtreenode(vector<nodedata*>& d);

		void Add(nodedata* d);
		void calpoint();
		bool is_leaf() const;
		void separate(size_t lim);
		
		vector<nodedata*> returnind(const nodedata& d) const
		{
			return returnind(d, constant<bool, true>());
		}

		template<typename boolean_allind>
		vector<nodedata*> returnind(const nodedata& d, const boolean_allind& allind) const;

		template<typename boolean_allind>
		static void pushind(vector<nodedata*>& t, const nodedata& d, const vector<nodedata*>& nd, const boolean_allind& allind);

		static flag_pair f0(const nodedata& d, const Point3D& cen);
	};

	class AABBtree :public VBOobj
	{
	public:
		AABBtreenode* root;
		vector<nodedata> sourcendata;

		virtual ~AABBtree();
		AABBtree(vector<nodedata>& s);
		AABBtree(vector<nodedata>&& s);

		void Buildtree();
		void reBuildtree();

		inddata returnind(const BoundBox& bbox) const;
		inddatas returnind() const;
		inddatas returnind(const AABBtree& rhs) const;
		inddata returnind(const Line3D& rhs) const;
	private:
		template<typename boolean_allind>
		inddatas returnind_impl(const AABBtree& rhs, const boolean_allind& allind) const;
	};
}

#endif
