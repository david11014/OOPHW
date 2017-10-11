/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <math.h>
#include <iostream>
#include <memory>  
using namespace std;
//#define DEBUG
#ifndef K_DTREE_H
#define K_DTREE_H

class Point2D {

public:
	union
	{
		double p[2];
		struct
		{
			double x;
			double y;
		};
	};
	int l;
	Point2D() {};
	Point2D(double a, double b, int c)
	{
		x = a;
		y = b;
		l = c;
	}

	double& operator[](unsigned int i)
	{
		return p[i];
	}
	Point2D operator+(Point2D P)
	{
		Point2D p;
		p.x = this->x + P.x;
		p.y = this->y + P.y;
		p.l = 0;
		return p;
	}
	Point2D operator-(Point2D P)
	{
		Point2D p;
		p.x = this->x - P.x;
		p.y = this->y - P.y;
		p.l = 0;
		return p;
	}
	Point2D operator*(double x)
	{
		Point2D p;
		p.x = (this->x) * x;
		p.y = (this->y) * x;
		p.l = 0;
		return p;
	}
	Point2D operator/(double x)
	{
		Point2D p;
		p.x = (this->x) / x;
		p.y = (this->y) / x;
		p.l = 0;
		return p;
	}

	double Distant(Point2D P)
	{
		return sqrt((p[0] - P[0])*(p[0] - P[0]) + (p[1] - P[1])*(p[1] - P[1]));
	}

	void show() {
		std::cout << p[0] << " " << p[1] << " " << l;
	}

	friend ostream& operator<<(ostream&, const Point2D&);
	
};

ostream& operator<<(ostream& os, const Point2D& p)
{
	os << p.x << "\t" << p.y << "\t" << p.l;
	return os;
}

class Node {

public:
	Node* up;
	Node* child[2];
	Point2D P;
	int layer;
	double axis;
	int visited;

	Node() 
	{
		child[0] = nullptr;
		child[1] = nullptr;
		visited = 0;
	}
	Node(Node* u)
	{
		up = u;
		child[0] = nullptr;
		child[1] = nullptr;
		visited = 0;
	}

	void show()
	{
		cout << P << "\taxis:"<< axis << " address:" << this << " layer:" << layer << " up:" << up << " child1:" << child[0] << " child2:" << child[1] << endl;		
	}

	friend ostream& operator<<(ostream&, const Node&);

	bool isLeaf()
	{
		if (child[0] == nullptr && child[1] == nullptr)
			return true;
		else 
			return false;
	}
};

ostream& operator<<(ostream& os, const Node& N)
{
	os << N.P << "\taxis:" << N.axis << " address:" << &N << " layer:" << N.layer << " up:" << N.up << " child1:" << N.child[0] << " child2:" << N.child[1] << endl;
	return os;
}

class KDTree {

public:
	Node* Root;	
	Point2D *trainP;

	KDTree(Point2D trainData[],int trSize)
	{
		trainP = new Point2D[trSize];

		memcpy(trainP, trainData, trSize * (sizeof(Point2D)));

		SortP(trainP, trSize, 0);

		int mid = (int)(trSize / 2);

		Root = new Node;		
		Root->layer = 0;
		Root->axis = (trainP[mid].p[0] + trainP[mid + 1].p[0]) / 2;
		Root->up = nullptr;
		

		MakeTree(Root, trainP, trainP + mid, 0);
		MakeTree(Root, trainP + mid, trainP + trSize, 1);
#ifdef DEBUG
		Root->show();
#endif 

		stack = new Node*[trSize * 2];
		
	}

	Node* Add(Node* pa, Point2D *P,int i)
	{
		Node* N(new Node);
		N->up = pa;
		N->P = *P;
		N->layer = pa->layer + 1;
		pa->child[i] = N;
		
		return N;
	}

	Point2D FindNear(Point2D P)
	{
		end = -1;
		ResetVisit(this->Root);
		
		double nearD = DBL_MAX;
		Point2D *nearP = new Point2D;
		FindNear(&P, this->Root, nearP, &nearD);	
		
		return *nearP;
	}
	
	void show()
	{
		show(Root);
	}

	void show(Node* N)
	{
		N->show();
		cout << endl;

		if (N->child[0] != nullptr)
			show(N->child[0]);

		if (N->child[1] != nullptr)
			show(N->child[1]);

		return;
	}

private:
	Node* *stack;
	int end = -1;

	void MakeTree(Node*parrent, Point2D * start, Point2D *end, int i) //來自child[i]
	{
		int n = (end - start);
		Node* N(new Node);
		N->up = parrent;
		N->layer = parrent->layer + 1;
		parrent->child[i] = N;

		if (n == 1)
		{				
			N->P = *start;	
		}
		else
		{	
			int d = (N->layer) % 2;
			SortP(start, n, d);

			int mid = (int)(n / 2);
			//求中位數
			if (n % 2 == 0)
			{
				N->axis = (start[mid].p[d] + start[mid - 1].p[d]) / 2;
			}
			else
			{
				N->axis = start[mid].p[d];
			}
			
			MakeTree(N, start, start + mid, 0);			
			MakeTree(N, start + mid, end, 1);
				
		}

		return;
	}

	void FindNear(Point2D *P, Node* N, Point2D *nearP, double *nearD)
	{
		if (N->isLeaf())
		{
			if (P->Distant(N->P) < *nearD)
			{
				*nearD = P->Distant(N->P);
				*nearP = N->P;
			}
		}
		else
		{		
			int d = (N->layer) % 2;

			if (P->p[d] < N->P[d]) //軸在左邊
			{

				if (P->p[d] - *nearD < N->axis && CanVisit(N->child[0])) 
					FindNear(P, N->child[0], nearP, nearD);

				if (P->p[d] + *nearD >= N->axis && CanVisit(N->child[1]))
					FindNear(P, N->child[1], nearP, nearD);

			}
			else
			{
				if (P->p[d] + *nearD >= N->axis && CanVisit(N->child[1]))
					FindNear(P, N->child[1], nearP, nearD);

				if (P->p[d] - *nearD < N->axis && CanVisit(N->child[0]))
					FindNear(P, N->child[0], nearP, nearD);
			}
		
		}
		
	}

	void SortP(Point2D* P, int n, int m){
		//x: m = 0 , y: m=1
		Point2D temp;
		for (int i = 0; i < n - 1; i++)
		{
			for (int j = 0; j < n - i - 1; j++)
			{
				if ((P[j])[m] > (P[j + 1])[m])
				{
					temp = (P[j]);
					(P[j]) = (P[j + 1]);
					(P[j + 1]) = temp;
				}
			}
		}
	}

	void Push(Node* N)
	{
#ifdef DEBUG
		std::cout << "push:" << end+1 << " ";
		N->show();
#endif 
		
		end++;
		stack[end] = N;
	}
	Node* Pop()
	{
#ifdef DEBUG
		std::cout << "pop:" << end+1 << " ";
		stack[end]->show();
#endif
		
		stack[end]->visited = 1;
		end--;
		return stack[end + 1];
	}

	bool CanVisit(Node* N)
	{
		if (N == nullptr)
			return false;
		else if (N->visited == 1)
			return false;
		else
			return true;
	}

	void ResetVisit(Node* N)
	{
		N->visited = 0;
		if (N->child[0] != nullptr)
			ResetVisit(N->child[0]);

		if (N->child[1] != nullptr)
			ResetVisit(N->child[1]);

		return;		
	}

	

};

#endif // !K_DTREE
