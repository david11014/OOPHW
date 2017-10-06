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

};

class Node {

public:
	Node* up;
	Node* child[2] = {nullptr};
	Point2D P;
	int layer;
	Node() {}
	Node(Node* u)
	{
		up = u;
	}
};

class KDTree {

public:
	Node *Root;

	KDTree()
	{
		Root = new Node;
		Root->layer = 0;
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
			Node *N = Add(parrent, start + mid, i);

			MakeTree(N, start, start + mid, 0);
			MakeTree(N, start + mid + 1, end, 1);

		}
		
		return;
	}

private:
	void SortP(Point2D *P, int n, int m){
		//x: m = 0 , y: m=1
		Point2D temp;
		for (int i = n; i > 0; i--)
		{
			for (int j = 0; j < i; j++)
			{
				if ((P[j])[m] < (P[j + 1])[m])
				{
					temp = (P[j]);
					(P[j]) = (P[j + 1]);
					(P[j + 1]) = temp;
				}
			}
		}
	}

};

#endif // !K_DTREE

