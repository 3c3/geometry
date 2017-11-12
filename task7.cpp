#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

// tolerance
#define EPSILON 0.01

struct Segment
{
	double x1, y1;
	double x2, y2;
};

struct Point
{
	double x, y;
	bool inFirst, inSecond; // does the point belong to the first/second triangle
};

bool comparePoints(Point *a, Point *b)
{
	if (a->x != b->x) return a->x < b->x;
	return a->y < b->y;
}

// Computes the determinant of this matrix
// x1 y1  1
// x2 y2  1
// x3 y3  1
// abs(det)/2 is the area of a triangle with the points from the matrix
// the sign of the result indicates the arrangement of the points (clockwise or not)
// if the result is 0, then the three points are on one line
double determinant(Point *p1, Point *p2, Point *p3)
{
	double dx2 = p2->x - p1->x;
	double dx3 = p3->x - p1->x;
	double dy2 = p2->y - p1->y;
	double dy3 = p3->y - p1->y;

	return dx2*dy3 - dx3*dy2;
}

bool isPointFromLineOnSegment(Point *p, Point *segPoint1, Point *segPoint2)
{
	double minX, maxX, minY, maxY;
	if (segPoint1->x < segPoint2->x)
	{
		minX = segPoint1->x;
		maxX = segPoint2->x;
	}
	else
	{
		minX = segPoint2->x;
		maxX = segPoint1->x;
	}

	if (segPoint1->y < segPoint2->y)
	{
		minY = segPoint1->y;
		maxY = segPoint2->y;
	}
	else
	{
		minY = segPoint2->y;
		maxY = segPoint1->y;
	}
	// we know the point is from the same line the segment is from
	return (p->x >= minX && p->x <= maxX) && (p->y >= minY && p->y <= maxY);
}

// Returns the intersection points of two line segments
// if the result is false, the segments do not intersect
// if the result is true, the intersection point has coordinates outX and outY
// if overlap is true, there is a second point in outX2 and outY2 (the segments are on the same line)
bool intersect(Segment *seg1, Segment *seg2, double *outX, double *outY, 
				double *outX2, double *outY2, bool *overlap)
{
	double dx1 = seg1->x1 - seg1->x2;
	double dy1 = seg1->y1 - seg1->y2;
	double dx2 = seg2->x1 - seg2->x2;
	double dy2 = seg2->y1 - seg2->y2;

	if (dx1 == 0 && dx2 == 0) // two vertical lines
	{
		if (seg1->x1 != seg2->x1) return false;

		double max1 = seg1->y1 > seg1->y2 ? seg1->y1 : seg1->y2;
		double min1 = seg1->y1 < seg1->y2 ? seg1->y1 : seg1->y2;

		double max2 = seg2->y1 > seg2->y2 ? seg2->y1 : seg2->y2;
		double min2 = seg2->y1 < seg2->y2 ? seg2->y1 : seg2->y2;

		if (max1 < min2) return false;
		if (max2 < min1) return false;
		// now we're sure there is intersection

		*outX = *outX2 = seg1->x1;
		
		// so we sort the y-coords and take the middle 2
		double sorted[4] = { min1, max1, min2, max2 };
		sort(sorted, sorted + 4);

		*outY = sorted[1];
		*outY2 = sorted[2];

		// if the second and third coord happen to be the same, it means we have only 1 point
		*overlap = !(sorted[1] == sorted[2]);
		return true;
	}
	else if (dx1 == 0) // first line is vertical, second one isn't
	{
		double x = seg1->x1;

		double minX = seg2->x1 < seg2->x2 ? seg2->x1 : seg2->x2;
		double maxX = seg2->x1 > seg2->x2 ? seg2->x1 : seg2->x2;

		if (x < minX || x > maxX) return false;

		*outX = x;

		double a = dy2 / dx2;
		double b = seg2->y1 - a*(seg2->x1);

		*outY = a*x + b;

		if (seg1->y1 < seg1->y2)
		{
			if (*outY < seg1->y1 || *outY > seg1->y2) return false;
		}
		else if (*outY < seg1->y2 || *outY > seg1->y1) return false;

		*overlap = false;
	}
	else if (dx2 == 0) // second line is vertical, first one isn't
	{
		double x = seg2->x1;

		double minX = seg1->x1 < seg1->x2 ? seg1->x1 : seg1->x2;
		double maxX = seg1->x1 > seg1->x2 ? seg1->x1 : seg1->x2;

		if (x < minX || x > maxX) return false;

		*outX = x;

		double a = dy1 / dx1;
		double b = seg1->y1 - a*(seg1->x1);

		*outY = a*x + b;

		if (seg2->y1 < seg2->y2)
		{
			if (*outY < seg2->y1 || *outY > seg2->y2) return false;
		}
		else if (*outY < seg2->y2 || *outY > seg2->y1) return false;

		*overlap = false;
	}
	else // both lines are not vertical
	{
		// calculate the line equations for both segments
		double a1 = dy1 / dx1;
		double b1 = seg1->y1 - a1*seg1->x1;

		double a2 = dy2 / dx2;
		double b2 = seg2->y1 - a2*seg2->x1;

		if (a1 == a2) // parallel lines
		{
			if (b1 != b2) return false;

			double start1X, start1Y;
			double end1X, end1Y;
			if (seg1->x1 < seg1->x2)
			{
				start1X = seg1->x1;
				end1X = seg1->x2;
				start1Y = seg1->y1;
				end1Y = seg1->y2;
			}
			else
			{
				start1X = seg1->x2;
				end1X = seg1->x1;
				start1Y = seg1->y2;
				end1Y = seg1->y1;
			}

			double start2X, start2Y;
			double end2X, end2Y;
			if (seg2->x1 < seg2->x2)
			{
				start2X = seg2->x1;
				end2X = seg2->x2;
				start2Y = seg2->y1;
				end2Y = seg2->y2;
			}
			else
			{
				start2X = seg2->x2;
				end2X = seg2->x1;
				start2Y = seg2->y2;
				end2Y = seg2->y1;
			}

			if (end1X < start2X) return false;
			if (end2X < start1X) return false;

			// sort the points so we can just take the middle two
			Point p1 = { seg1->x1 , seg1->y1 };
			Point p2 = { seg1->x2 , seg1->y2 };
			Point p3 = { seg2->x1 , seg2->y1 };
			Point p4 = { seg2->x2 , seg2->y2 };
			Point* sorted[4] = { &p1, &p2, &p3, &p4 };
			sort(sorted, sorted + 4, comparePoints);

			*outX = sorted[1]->x;
			*outY = sorted[1]->y;
			*outX2 = sorted[2]->x;
			*outY2 = sorted[2]->y;

			// if the two midpoints are the same, then we shouldn't care about the second
			*overlap = !(*outX == *outX2 && *outY == *outY2);
			return true;
		}

		// find the intersection
		*outX = (b2 - b1)/(a1 - a2);
		*outY = (*outX)*a1 + b1;

		double minX = seg1->x1 < seg1->x2 ? seg1->x1 : seg1->x2;
		double maxX = seg1->x1 > seg1->x2 ? seg1->x1 : seg1->x2;
		if (*outX < minX) return false;
		if (*outX > maxX) return false;

		minX = seg2->x1 < seg2->x2 ? seg2->x1 : seg2->x2;
		maxX = seg2->x1 > seg2->x2 ? seg2->x1 : seg2->x2;
		if (*outX < minX) return false;
		if (*outX > maxX) return false;

		*overlap = false;
		return true;
	}
}

bool isPointInTriangle(Point *p, Point *tp1, Point *tp2, Point *tp3)
{
	double det1 = determinant(tp1, tp2, p);
	if (det1 == 0) return isPointFromLineOnSegment(p, tp1, tp2);
	double det2 = determinant(tp2, tp3, p);
	if (det2 == 0) return isPointFromLineOnSegment(p, tp2, tp3);
	double det3 = determinant(tp3, tp1, p);
	if (det3 == 0) return isPointFromLineOnSegment(p, tp3, tp1);

	if (det1 < 0) return det2 < 0 && det3 < 0;
	else return det2 > 0 && det3 > 0;
}

double findPolygonArea(Point* points, int nPoints)
{
	Point** sorted = new Point*[nPoints+1];
	for (int i = 0; i < nPoints; i++) sorted[i] = points + i;

	// sort the points according to their x coordinate
	sort(sorted, sorted + nPoints, comparePoints);

	Point** hull = new Point*[nPoints+1];
	int k = 0;
	int hullSize;

	// finds the convex hull of the given points
	// (graham monotonne chain)
	for (int i = 0; i < nPoints; ++i)
	{
		while (k >= 2 && (determinant(hull[k-2], hull[k-1], sorted[i]) <= 0)) k--;
		hull[k++] = sorted[i];
	}

	for (int i = nPoints - 2, t = k + 1; i >= 0; i--)
	{
		while (k >= t && (determinant(hull[k - 2], hull[k - 1], sorted[i]) <= 0)) k--;
		hull[k++] = sorted[i];
	}

	hullSize = k;
	hullSize--; // the first and last point are the same

	#ifdef INFO
	cout << "Convex hull:" << endl;
	for (int i = 0; i < hullSize; i++)
	{
		cout << hull[i]->x << "; " << hull[i]->y << endl;
	}
	#endif

	double area = 0;
	for (int i = 2; i < hullSize; i++)
	{
		double triangleArea = determinant(hull[0], hull[i - 1], hull[i]);
		if (triangleArea < 0) area -= triangleArea;
		else area += triangleArea;
	}

	area /= 2.0;

	delete[] sorted;
	delete[] hull;

	return area;
}

bool isLevel(Point *a, Point *b, Point *c)
{
	if (a->y == b->y) return true;
	if (a->y == c->y) return true;
	if (b->y == c->y) return true;
	return false;
}

double distance(Point *a, Point *b)
{
	double dx = a->x - b->x;
	double dy = a->y - b->y;
	return sqrt(dx*dx + dy*dy);
}

bool isEquiliteral(Point *a, Point *b, Point *c)
{
	double d1 = distance(a, b);
	double d2 = distance(a, c);
	double d3 = distance(b, c);

	if (abs(d1 - d2) > EPSILON) return false;
	if (abs(d1 - d3) > EPSILON) return false;
	if (abs(d2 - d3) > EPSILON) return false;

	return true;
}

int main()
{
	Point points[20];
	Point* triangle1 = points;
	Point* triangle2 = points+3;

	// Input
	for (int i = 0; i < 3; i++)
	{
		cin >> triangle1[i].x >> triangle1[i].y;
		triangle1[i].inFirst = true;
		triangle1[i].inSecond = false;
	}

	// output the error message quickly
	if (!isLevel(triangle1, triangle1 + 1, triangle1 + 2) || !isEquiliteral(triangle1, triangle1 + 1, triangle1 + 2))
	{
		cout << "Input does not satisfy conditions" << endl;
		return 0;
	}

	for (int i = 0; i < 3; i++)
	{
		cin >> triangle2[i].x >> triangle2[i].y;
		triangle2[i].inFirst = false;
		triangle2[i].inSecond = true;
	}
	// End of input

	
	if (!isLevel(triangle2, triangle2 + 1, triangle2 + 2) || !isEquiliteral(triangle2, triangle2 + 1, triangle2 + 2))
	{
		cout << "Input does not satisfy conditions" << endl;
		return 0;
	}
	

	// Checks whether points from one triangle belong to the other as well
	for (int i = 0; i < 3; i++)
	{
		Point *current = triangle1 + i;
		if (isPointInTriangle(current, &triangle2[0], &triangle2[1], &triangle2[2]))
		{
			current->inSecond = true;
		}
	}

	for (int i = 0; i < 3; i++)
	{
		Point *current = triangle2 + i;
		if (isPointInTriangle(current, &triangle1[0], &triangle1[1], &triangle1[2]))
		{
			current->inFirst = true;
		}
	}

	// Adds all the intersection points
	int pointIdx = 6;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Segment seg1;
			seg1.x1 = triangle1[i].x;
			seg1.y1 = triangle1[i].y;
			seg1.x2 = triangle1[(i+1)%3].x;
			seg1.y2 = triangle1[(i+1)%3].y;

			Segment seg2;
			seg2.x1 = triangle2[j].x;
			seg2.y1 = triangle2[j].y;
			seg2.x2 = triangle2[(j+1)%3].x;
			seg2.y2 = triangle2[(j+1)%3].y;

			double x1,x2,y1,y2;
			bool overlap;
			if (intersect(&seg1, &seg2, &x1, &y1, &x2, &y2, &overlap))
			{
				Point *current = points + pointIdx;
				current->x = x1;
				current->y = y1;
				current->inFirst = true;
				current->inSecond = true;
				pointIdx++;

				#ifdef INFO
				cout << "Found intersection: " << x1 << " " << y1 << endl;
				#endif

				if (overlap)
				{
					current = points + pointIdx;
					current->x = x2;
					current->y = y2;
					current->inFirst = true;
					current->inSecond = true;
					pointIdx++;

					#ifdef INFO
					cout << "Found another intersection: " << x2 << " " << y2 << endl;
					#endif
				}
			}
		}
	}

	Point *ptrPoints[20];
	for (int i = 0; i < 20; i++) ptrPoints[i] = points + i;

	// sort the points so we can remove duplicates
	sort(ptrPoints, ptrPoints + pointIdx, comparePoints);

	Point sortedPoints[20];
	int sortedIdx = 0;
	sortedPoints[0] = *ptrPoints[0];

	// remove duplicates
	for (int i = 1; i < pointIdx; i++)
	{
		Point *current = ptrPoints[i];
		Point *sorted = sortedPoints + sortedIdx;
		if (current->x == sorted->x && current->y == sorted->y)
		{ // if two points are the same, we add up the flags
			sorted->inFirst |= current->inFirst;
			sorted->inSecond |= current->inSecond;
		}
		else sortedPoints[++sortedIdx] = *current;
	}
	sortedIdx++;

	#ifdef INFO
	cout << sortedIdx << " points to scan" << endl;
	#endif

	// put in only the points that belong to the cross section
	Point commonPoints[10];
	int commonIdx = 0;
	for (int i = 0; i < sortedIdx; i++)
	{
		Point *current = sortedPoints + i;
		if (current->inFirst && current->inSecond)
			commonPoints[commonIdx++] = *current;
	}

	#ifdef INFO
	cout << "Found " << commonIdx << " common points" << endl;
	for (int i = 0; i < commonIdx; i++)
	{
		cout << commonPoints[i].x << '\t' << commonPoints[i].y << endl; 
	}
	#endif

	if (commonIdx < 3) // We don't have enough points for an area
	{
		#ifdef INFO
		cout << "Less than 3 common points => area is " << endl;
		#endif
		cout << 0 << endl;
	}
	else if (commonIdx == 3) // the cross section is a triangle, easy
	{
		double area = determinant(commonPoints, commonPoints + 1, commonPoints + 2) / 2.0;
		if (area < 0) area = -area;
		
		#ifdef INFO
		cout << "Common area is a triangle with area " << endl;
		#endif
		cout << area << endl;
	}
	else // the cross section is a polygon, gotta use a convex hull
	{
		double area = findPolygonArea(commonPoints, commonIdx);

		#ifdef INFO
		cout << "More than 3 common points, applying convex hull" << endl;
		#endif		
		cout << area << endl;
	}

	return 0;
}