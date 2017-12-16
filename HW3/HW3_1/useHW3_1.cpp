/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 3
Write by david1104
github: https://github.com/david11014
********************************************************/
#include <iostream>
#include <random>
#include <ctime>
#include "HW3_1.h"

int main()
{
	srand(time(NULL));

	do{

		IGeometry** arr = new IGeometry*[3];
		arr[0] = new Pyramid();
		arr[1] = new Cuboid();
		arr[2] = new Cylinder();
		Pyramid *Py = dynamic_cast<Pyramid*>(arr[0]);
		Cuboid *Cu = dynamic_cast<Cuboid*>(arr[1]);
		Cylinder *Cy = dynamic_cast<Cylinder*>(arr[2]);


		cout << "Pyramid\nvertices\n";
		Py->show();
		cout << "Area:" << arr[0]->Area() << " Perimete:" << arr[0]->Perimeter() << " Volume:" << arr[0]->Volume() << endl;
		cout << "Center:" << Py->Center() << endl;
		cout << "====================================================\n";

		cout << "Cuboid\nvertices\n";
		Cu->show();
		cout << "Area:" << arr[1]->Area() << " Perimete:" << arr[1]->Perimeter() << " Volume:" << arr[1]->Volume() << endl;
		double *l = Cu->SideLength();
		double *A = Cu->SideArea();
		cout << "SideLength:" << l[0] << " " << l[1] << " " << l[2] << " " << endl;
		cout << "SideArea:" << A[0] << " " << A[1] << " " << A[2] << " " << endl;
		cout << "====================================================\n";

		cout << "Cylinder\n";
		Cy->show();
		cout << "Area:" << arr[2]->Area() << " Perimete:" << arr[2]->Perimeter() << " Volume:" << arr[2]->Volume() << endl;
		cout << "Height:" << Cy->Height() << " BottomArea:" << Cy->BottomArea() << " SideArea:" << Cy->SideArea() << endl;
		cout << "====================================================\n\n";

		delete[] l;
		delete[] A;

		cout << "input enter to generate new data\ninput 'q' to exit\n";
	} while ('q' != getchar());

		
	return 0;
}