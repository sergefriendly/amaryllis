#include <iostream>
#include "abcd.h"

class CurveTest {
public:
	std::vector<Point> points = {
			{2.0, 6.0, 5.0},
			{0.0, 1.0, 1.0},
			{3.0, 7.0, 11.0},
			{1.0, 4.0, 3.0},
	};

	std::vector<Segment> segments = {
			{3, 4, 5, 6, 0, 4},
			{6, 3, 2, 7, 0, 4},
			{7, 3, 8, 2, 0, 4},
	};

	void pointTest() {
		for (const auto& p : points) {
			std::cout << p << std::endl;
		}
	}

	void segmentTest() {
		for (const auto& s : segments) {
			std::cout << s << std::endl;
		}
	}

	void curveSortTest() {
		Curve curve(points);
		for (const auto& p : curve.ordpoints) {
			std::cout << p << std::endl;
		}
	}

};



int main()
{
    std::cout << "Hello World!\n";
	CurveTest ct;
	ct.pointTest();
	std::cout << "Good." << std::endl;
	ct.segmentTest();
	std::cout << "Good." << std::endl;
	ct.curveSortTest();
	std::cout << "Good." << std::endl;



}
