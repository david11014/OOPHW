#pragma once

#ifndef __OctCloud_H__
#define __OctCloud_H__

#include "BoundBox.h"
#include "myModelbase.h"
#include "VBOobj.h"

namespace myModel
{
	struct triid_pid
	{
		// 	struct
		// 	{
		// 		unsigned second : 2;
		// 		unsigned first : 30;
		// 	};
		size_t v;

		triid_pid() :v(0){}
		triid_pid(const triid_pid& t) :v(t.v){}
		triid_pid(const size_t a) :v(a){}
		triid_pid(const tri_ind_t a, const vertex_ind_t b) :v(a * 4 + b){}
		bool operator< (const triid_pid& t) const
		{
			return v / 4 < t.v / 4;
		}
		bool operator== (const triid_pid& t) const
		{
			return v / 4 == t.v / 4;
		}
		tri_ind_t tid() const
		{
			return v / 4;
		}
		vertex_ind_t pid() const
		{
			return v & 0x3;
		}
	};

	class Point3DNode
	{
	private:
		Point3DNode(){};
	public:
		Point3DNode(myMath::Point3D* Data, const triid_pid& tpid);
		myMath::Point3D pp;
		myMath::Point3D* pdata;
		triid_pid triind_;
		std::vector<Point3DNode*> adj_Point;
		std::vector<triid_pid> adj_tri;
		
		inddata find_adjtri(size_t turns) const;
		inddatas find_adjtris(size_t turns) const;
	};

	class OctCloudNode :public BoundBox
	{
	public:
		OctCloudNode* Nodeptr;
		std::vector<Point3DNode*> p;

		myMath::Point3D cen;
		
		OctCloudNode();
		OctCloudNode(const OctCloudNode&) = delete;
		~OctCloudNode();
		OctCloudNode(vector<Point3DNode*>& d);

		void Add(Point3DNode* d);
		void calpoint();
		void calpoint2();
		bool is_leaf() const;
		void separate(size_t lim);

		bool isinbox(const myMath::Point3D& P, double lim = 10) const;
		void remove_duplicate_point();
		//void remove_duplicate_point2();
		void remove_duplicate_point3();
	};

	class OctCloud :public VBOobj
	{
	public:
		bool hasduplicate;
		OctCloudNode* root;
		std::vector<Point3DNode*> pnarr;
		std::vector<OctCloudNode*> LeafNode;

		virtual ~OctCloud();
		explicit OctCloud(std::vector<Point3DNode*>& PN);

		Point3DNode* searchNode(const myMath::Point3D& P, double lim = 10);
		void segmentation();
		void resegmentation();
	};
}
#endif