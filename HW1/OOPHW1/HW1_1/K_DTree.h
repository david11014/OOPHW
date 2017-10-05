#ifndef K_DTREE
#define K_DTREE

class KDTree {
	
	node root;

	void add(node* pa, Point2D * P)
	{
		
	}

};

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
	node* left;
	node* right;
	Point2D* P;
	node() {
		left = nullptr;
		right = nullptr;
		P = nullptr;
	}
	node(node* u)
	{		
		left = nullptr;
		right = nullptr;
		P = nullptr;

		up = u;
	}
};

#endif // !K_DTREE

