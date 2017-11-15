#include <iostream>
#include <windows.h>

using namespace std;

#define EPSILON 0.01
#define HALF_SQRT3 0.8660254

struct Point
{
	double x, y;
};

double determinant(Point *p1, Point *p2, Point *p3)
{
	double dx2 = p2->x - p1->x;
	double dx3 = p3->x - p1->x;
	double dy2 = p2->y - p1->y;
	double dy3 = p3->y - p1->y;

	return dx2*dy3 - dx3*dy2;
}

bool isPointInTriangle(Point *p, Point *tp1, Point *tp2, Point *tp3)
{
	double det1 = determinant(tp1, tp2, p);
	if (det1 == 0)
	{
		double minX = tp1->x < tp2->x ? tp1->x : tp2->x;
		double maxX = tp1->x > tp2->x ? tp1->x : tp2->x;
		return p->x >= minX && p->x <= maxX;
	}
	double det2 = determinant(tp2, tp3, p);
	if (det2 == 0)
	{
		double minX = tp2->x < tp3->x ? tp2->x : tp3->x;
		double maxX = tp2->x > tp3->x ? tp2->x : tp3->x;
		return p->x >= minX && p->x <= maxX;
	}
	double det3 = determinant(tp3, tp1, p);
	if (det3 == 0)
	{
		double minX = tp3->x < tp1->x ? tp3->x : tp1->x;
		double maxX = tp3->x > tp1->x ? tp3->x : tp1->x;
		return p->x >= minX && p->x <= maxX;
	}

	if (det1 < 0) return det2 < 0 && det3 < 0;
	else return det2 > 0 && det3 > 0;
}

// checks if a triangle (array of 3 points) meets the requirements
bool isCool(Point* triangle)
{
	Point *p1 = triangle;
	Point *p2 = triangle + 1;
	Point *p3 = triangle + 2;
	if (abs(p1->y - p2->y) < EPSILON)
	{
		if (p3->y < p1->y) return false;
		if (abs(p1->x + p2->x - 2.0*p3->x) / 2.0 > EPSILON) return false; // is the top exactly between the other two vertices
		return (abs(abs(p1->x - p2->x)*HALF_SQRT3 + p1->y - p3->y) < EPSILON); // is the top at the proper height
	}
	else if (abs(p1->y - p3->y) < EPSILON)
	{
		if (p2->y < p1->y) return false;
		if (abs(p1->x + p3->x - 2.0*p2->x) / 2.0 > EPSILON) return false;
		return (abs(abs(p1->x - p3->x)*HALF_SQRT3 + p1->y - p2->y) < EPSILON);
	}
	else if (abs(p2->y - p3->y) < EPSILON)
	{
		if (p1->y < p2->y) return false;
		if (abs(p2->x + p3->x - 2.0*p1->x) / 2.0 > EPSILON) return false;
		return (abs(abs(p2->x - p3->x)*HALF_SQRT3 + p2->y - p1->y) < EPSILON);
	}
	return false;
}

int main()
{
	Point triangle1[3];
	Point triangle2[3];

	double minX = 99999;
	double maxX = -99999;
	double minY = 99999;
	double maxY = -99999;

	for (int i = 0; i < 3; i++)
	{
		cin >> triangle1[i].x >> triangle1[i].y;
		if (triangle1[i].x > maxX) maxX = triangle1[i].x;
		if (triangle1[i].x < minX) minX = triangle1[i].x;
		if (triangle1[i].y > maxY) maxY = triangle1[i].y;
		if (triangle1[i].y < minY) minY = triangle1[i].y;
	}

	if (!isCool(triangle1))
	{
		cout << "Input does not satisfy requirements" << endl;
		return 0;
	}

	for (int i = 0; i < 3; i++)
	{
		cin >> triangle2[i].x >> triangle2[i].y;
		if (triangle2[i].x > maxX) maxX = triangle2[i].x;
		if (triangle2[i].x < minX) minX = triangle2[i].x;
		if (triangle2[i].y > maxY) maxY = triangle2[i].y;
		if (triangle2[i].y < minY) minY = triangle2[i].y;
	}

	if (!isCool(triangle2))
	{
		cout << "Input does not satisfy requirements" << endl;
		return 0;
	}

	double eps = 0.007;
	//cout << "Enter the epsilon: ";
	//cin >> eps;

	double epsArea = eps*eps;
	double area = 0;

	Point testPoint = { minX + eps / 2.0, minY + eps / 2.0 };

	while (testPoint.x < maxX)
	{
		testPoint.y = minY + eps / 2.0;
		while (testPoint.y < maxY)
		{
			if (isPointInTriangle(&testPoint, triangle1, triangle1 + 1, triangle1 + 2) &&
				isPointInTriangle(&testPoint, triangle2, triangle2 + 1, triangle2 + 2))
				area += epsArea;
			testPoint.y += eps;
		}
		testPoint.x += eps;
	}

	cout << area << endl;
	//cout << "Area: " << area << endl;
	//system("pause");
	return 0;
}