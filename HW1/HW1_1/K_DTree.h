/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 1
Write by david1104
********************************************************/
#include<math.h>
#include<iostream>
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

	double Distant(Point2D P)
	{
		return sqrt((p[0] - P[0])*(p[0] - P[0]) + (p[1] - P[1])*(p[1] - P[1]));
	}

	void show() {
		std::cout << p[0] << " " << p[1] << " " << l;
	}
};

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
};

class KDTree {

public:
	Node *Root;


	KDTree(Point2D trainP[],int trSize)
	{
		int mid = (int)(trSize / 2);
		Root = new Node;
		Root->P = trainP[mid];
		Root->layer = 0;
		Root->up = nullptr;

		MakeTree(Root, trainP, trainP + trSize + 1, 0);

		stack = new Node*[trSize];
	}

	Node* Add(Node* pa, Point2D *P,int i)
	{
		Node *N = new Node;
		N->P = *P;
		N->layer = pa->layer + 1;
		pa->child[i] = N;

		return N;
	}

	void MakeTree(Node *parrent, Point2D * start, Point2D *end,int i)
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
			Node *N = Add(parrent, start + mid, i);//add divide point in node
			i = (i == 0 ? 1 : 0); //change sort x or y in next node
			MakeTree(N, start, start + mid, i);
			MakeTree(N, start + mid + 1, end, i);

		}
		
		return;
	}

	Point2D FindNear(Point2D P)
	{
		return FindNear(&P, this->Root);
	}
	//Point2D FindNear(Point2D *P, Node* N)
	//{
	//	//leaf
	//	if (N->child[0] == nullptr&&N->child[1] == nullptr)
	//	{			
	//		return N->P;
	//	}
	//	else //node
	//	{
	//		int d = (N->layer) % 2; //0:x , 1:y
	//		Point2D nearP;
	//		if (P->p[d] > N->P[d])
	//		{
	//			nearP = FindNear(P, N->child[1]);
	//			double nd = P->Distant(nearP);
	//			if (nd > abs(P->p[d] - N->P[d]))
	//			{
	//				Point2D nearPTemp = FindNear(P, N->child[0]);
	//				if (P->Distant(nearP) > P->Distant(nearPTemp))
	//					nearP = nearPTemp;
	//			}
	//		}				
	//		else
	//		{
	//			nearP = FindNear(P, N->child[0]);
	//			double nd = P->Distant(nearP);
	//			if (nd > abs(P->p[d] - N->P[d]))
	//			{
	//				Point2D nearPTemp = FindNear(P, N->child[1]);
	//				if (P->Distant(nearP) > P->Distant(nearPTemp))
	//					nearP = nearPTemp;
	//			}
	//		}
	//		return nearP;
	//	}
	//}

	Point2D FindNear(Point2D *P, Node* N)
	{
		Point2D nearP;
		Push(N);
		do {
			
			
		} while (end == -1);
	}
private:
	Node* *stack;
	int end = -1;
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
		end++;
		stack[end] = N;
	}
	Node* Pop()
	{
		end--;
		return stack[end + 1];
	}

};

#endif // !K_DTREE

