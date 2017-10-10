/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <math.h>
#include <iostream>
#include <memory>  
using namespace std;
#ifndef K_DTREE
#define K_DTREE

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
	os << p.x << ", " << p.y << " " << p.l;
	return os;
}

class Node {

public:
	Node* up;
	Node* child[2];
	Point2D P;
	int layer;
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
		cout << " address:" << this << " layer:" << layer << " up:" << up << " child1:" << child[0] << " child2:" << child[1] << endl;
	}
};

class KDTree {

public:
	Node* Root;
	
	KDTree(Point2D trainP[],int trSize)
	{
		int mid = (int)(trSize / 2);
		Root = new Node;
		Root->P = trainP[mid];
		Root->layer = 0;
		Root->up = nullptr;
		
		MakeTree(Root, trainP, trainP + mid, 0);
		MakeTree(Root, trainP + mid + 1, trainP + trSize + 1, 1);
#ifdef DEBUG
		Root->show();
#endif 

		stack = new Node*[trSize];
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
		return FindNear(&P, this->Root);
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
			Node*N = Add(parrent, start + mid, i);//add divide point in node
			i = (i == 0 ? 1 : 0); //change sort x or y in next node
			MakeTree(N, start, start + mid, 0);
			MakeTree(N, start + mid + 1, end, 1);
		}

		return;
	}

	Point2D FindNear(Point2D *P, Node* N)
	{
		Point2D nearP = N->P;
		Push(N);
		do {

			//判斷是否可往下走，
			if (CanVisit(stack[end]->child[0]) == false && CanVisit(stack[end]->child[1]) == false)
			{
				//更新最近點
				if (stack[end]->P.Distant(nearP) < P->Distant(nearP))
				{
					nearP = stack[end]->P;
				}
				Pop();
			}
			else if (CanVisit(stack[end]->child[0]) == false && CanVisit(stack[end]->child[1]) == true)
			{
				double D = P->Distant(stack[end]->child[1]->P);

				if (P->Distant(nearP) > D)
				{
					Push(stack[end]->child[1]);
				}
				else
				{
					Pop();
				}
			}
			else if (CanVisit(stack[end]->child[0]) == true && CanVisit(stack[end]->child[1]) == false)
			{
				double D = P->Distant(stack[end]->child[0]->P);

				if (P->Distant(nearP) > D)
				{
					Push(stack[end]->child[0]);
				}
				else
				{
					Pop();
				}
			}
			else
			{
				int d = (stack[end]->layer) % 2;
				if (P->p[d] < stack[end]->P[d])
				{
					if (stack[end]->child[0]->visited == 0)
						Push(stack[end]->child[0]);
				}
				else
				{
					if (stack[end]->child[1]->visited == 0)
						Push(stack[end]->child[1]);
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
		std::cout << "push:" << end << " ";
		N->show();
#endif 
		
		end++;
		stack[end] = N;
	}
	Node* Pop()
	{
#ifdef DEBUG
		std::cout << "pop:" << end << " ";
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
