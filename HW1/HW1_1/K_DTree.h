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

class node {

public:
	node* up;
	node* child[2] = {nullptr};
	Point2D* P;
	int layer;
	node() {		
		P = nullptr;
	}
	node(node* u)
	{		
		P = nullptr;
		up = u;
	}
};

class KDTree {

public:
	node Root;

	KDTree() {
		Root.layer = 0;
	}

	void Add(node* pa, Point2D *P,int i)
	{
		node *N = new node(pa);
		N->P = P;
		pa->child[i] = N;
	}

	void MakeTree(node *parrent, Point2D * start, Point2D *end)
	{

	}

private:
	void SortP(Point2D *P, int n, int m)
	{
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

