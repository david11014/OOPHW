/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <math.h>
#include <iostream>
#include <memory>  
using namespace std;
#define DEBUG
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

	bool isLeaf()
	{
		if (child[0] == nullptr && child[1] == nullptr)
			return true;
		else 
			return false;
	}
};

class KDTree {

public:
	Node* Root;
	
	KDTree(Point2D trainP[],int trSize)
	{
		SortP(trainP, trSize, 0);

		int mid = (int)(trSize / 2);

		Root = new Node;		
		Root->layer = 0;
		Root->axis = trainP[mid].p[0];
		Root->up = nullptr;		
		

		MakeTree(Root, trainP, trainP + mid, 0);
		MakeTree(Root, trainP + mid, trainP + trSize + 1, 1);
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

		//N->show();
		

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

	void MakeTree(Node*parrent, Point2D * start, Point2D *end, int i)
	{
		int n = (end - start);

		if (n == 1)
		{
			Add(parrent, start, i);
			return;
		}
		else if (n < 1)
		{
			return;
		}
		else
		{
			int mid = (int)(n / 2);

			SortP(start, n, i);
			Node*N = new Node;
			N->up = parrent;			
			N->layer = parrent->layer + 1;
			int d = (N->layer) % 2;
			N->axis = start[mid].p[d];
			parrent->child[i] = N;

			i = (i == 0 ? 1 : 0); //change sort x or y in next node
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

	Point2D FindNear(Point2D *P, Node* N)
	{
		Point2D nearP = N->P;
		double nearD = DBL_MAX;
		Push(N);
		do {

			//判斷是否可往下走
			if (stack[end]->isLeaf())
			{
				//更新最近點
				if (stack[end]->P.Distant(nearP) < nearD)
				{
					nearD = stack[end]->P.Distant(nearP);
					nearP = stack[end]->P;
				}
				Pop();
			}			
			else {
				int d = (stack[end]->layer) % 2;
				if (stack[end]->axis > P->p[d]) //點在軸左邊
				{
					if (P->p[d] - nearD <= stack[end]->axis && CanVisit(stack[end]->child[0]))
						Push(stack[end]->child[0]);

					if (P->p[d] + nearD > stack[end]->axis && CanVisit(stack[end]->child[1]))
						Push(stack[end]->child[1]);
					else
						Pop();


				}
				else //點在軸右邊
				{
					if (P->p[d] + nearD > stack[end]->axis && CanVisit(stack[end]->child[1]))
						Push(stack[end]->child[0]);

					if (P->p[d] - nearD <= stack[end]->axis && CanVisit(stack[end]->child[0]))
						Push(stack[end]->child[1]);
					else
						Pop();
				}				
			}

		} while (end != -1);
		return nearP;
	}

	void SortP(Point2D *P, int n, int m){
		//x: m = 0 , y: m=1
		Point2D temp;
		for (int i = n; i > 0; i--)
		{
			for (int j = 0; j < i; j++)
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
