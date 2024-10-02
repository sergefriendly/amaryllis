#ifndef SPLINES_H_
#define SPLINES_H_

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


struct Point {
	double t, x, y;
	Point(double t, double x, double y) : t(t), x(x), y(y) {}

	friend std::ostream& operator<<(std::ostream& os, const Point& point) {
		os << "( " << point.t << ": " << point.x << ", " << point.y << " )";
		return os;
	}
};


struct Segment {
	double a, b, c, d, t_min, t_max;
	Segment() {}
	Segment(double a, double b, double c, double d, double t_min, double t_max)
		: a(a), b(b), c(c), d(d), t_min(t_min), t_max(t_max) {}

	friend std::ostream& operator<<(std::ostream& os, const Segment& segment) {
		os << "( a = " << segment.a << ", b = " << segment.b
			<< ", c = " << segment.c << ", d = " << segment.d
			<< ", [" << segment.t_min << ", " << segment.t_max << "] )";
		return os;
	}
};


class Curve {
public:
	std::vector<Point> points;
	std::vector<Segment> segments;

	// Метод прогонки
	std::vector<double> solveTridiagonal(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
		size_t n = b.size();
		std::vector<double> x(n), p(n), q(n);

		p[0] = A[0][1] / A[0][0];
		q[0] = b[0] / A[0][0];

		for (size_t i = 1; i < n - 1; i++) {
			p[i] = A[i][i + 1] / (A[i][i] - A[i][i -1] * p[i - 1]);
			q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] - A[i][i - 1] * p[i - 1]);
		}

		x[n - 1] = (b[n - 1] - A[n - 1][n - 2] * q[n - 2]) / (A[n - 1][n - 1] - A[n - 1][n - 2] * p[n - 2]);

		for (int i = n - 2; i >= 0; i--) {
			x[i] = q[i] - p[i] * x[i + 1];
		}

		return x;
	}

	Curve(const std::vector<Point>& _points) {
		if (_points.size() < 2) {
			throw std::invalid_argument("Необходимы как минимум две точки.");
		}

		points = _points;

		std::sort(points.begin(), points.end(), [](const Point& p1, const Point& p2) {
			return p1.t < p2.t;
		});

		size_t N_pnt = points.size();
		size_t N_seg = points.size() - 1;

		std::vector<double> h(N_seg);
		for (size_t i = 0; i < N_seg; i++) {
			h[i] = points[i + 1].t - points[i].t;
			//std::cout << "h[" << i << "] = " << h[i] << std::endl;
			if (h[i] <= 0) {
				throw std::invalid_argument("После сортировки точки должны быть в возрастании параметра t.");
			}
		}

		std::vector<std::vector<double>> A(N_pnt, std::vector<double>(N_pnt, 0));
		std::vector<double> b(N_pnt, 0);

		A[0][0] = 1;
		A[N_seg][N_seg] = 1;

		for (size_t i = 1; i < N_seg; i++) {
			A[i][i - 1] = h[i - 1];
			A[i][i] = 2 * (h[i - 1] + h[i]);
			A[i][i + 1] = h[i];
			b[i] = 3 * ((points[i + 1].x - points[i].x) / h[i] - (points[i].x - points[i - 1].x) / h[i - 1]);
		}

		std::vector<double> c = solveTridiagonal(A, b);

		segments.resize(N_seg);
		for (size_t i = 0; i < N_seg; i++) {
			double a = points[i].x;
			double b = (points[i + 1].x - points[i].x) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
			double d = (c[i + 1] - c[i]) / (3 * h[i]);
			segments[i] = { a, b, c[i], d, points[i].t, points[i + 1].t };
		}
	}

	double interp_x_by_t(double t) const {
		if (t < points.front().t || t > points.back().t) {
			throw std::out_of_range("t не в интерполяционном интервале.");
		}

		// Бинарный поиск
		auto it = std::lower_bound(points.begin(), points.end(), t, [](const Point& p, double val) {return p.t < val; });

		// Индекс перед элементом it
		size_t i = std::distance(points.begin(), it) - 1;

		const Segment& seg = segments[i];
		double dt = t - seg.t_min;
		return seg.a + seg.b * dt + seg.c * dt * dt + seg.d * dt * dt * dt;
	}

	double interp_y_by_x(double x) const {
		if (x < points.front().x || x > points.back().x) {
			throw std::out_of_range("x не в интерполяционном интервале.");
		}

		auto it = std::lower_bound(points.begin(), points.end(), x, [](const Point& p, double val) {return p.x < val; });
		size_t i = std::distance(points.begin(), it) - 1;

		const Point& p1 = points[i];
		const Point& p2 = points[i + 1];
		double t = p1.t + (x - p1.x) * (p2.t - p1.t) / (p2.x - p1.x);
		return p1.y + (p2.y - p1.y) * (t - p1.t) / (p2.t - p1.t);
	}

	static std::vector<Point> findIntersections(const Curve& curve1, const Curve& curve2, double epsilon = 1e-10) {
		std::vector<Point> intersections;
		double t_min = std::max(curve1.points.front().t, curve2.points.front().t);
		double t_max = std::min(curve1.points.back().t, curve2.points.back().t);

		for (double t = t_min; t <= t_max; t += epsilon) {
			double x1 = curve1.interp_x_by_t(t);
			double y1 = curve1.interp_y_by_x(x1);
			double x2 = curve2.interp_x_by_t(t);
			double y2 = curve2.interp_y_by_x(x2);

			if (std::abs(x1 - x2) < epsilon && std::abs(y1 - y2) < epsilon) {
				intersections.push_back({t, x1, y1});
			}
		}
		return intersections;
	}

	static double findMinimalDistance(const Curve& curve1, const Curve& curve2, double epsilon = 1e-18) {
		double min_distance = std::numeric_limits<double>::max();
		double t_min = std::max(curve1.points.front().t, curve2.points.front().t);
		double t_max = std::min(curve1.points.back().t, curve2.points.back().t);

		for (double t = t_min; t <= t_max; t += epsilon) {
			double x1 = curve1.interp_x_by_t(t);
			double y1 = curve1.interp_y_by_x(x1);
			double x2 = curve2.interp_x_by_t(t);
			double y2 = curve2.interp_y_by_x(x2);

			double distance = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
			min_distance = std::min(min_distance, distance);
		}
		return min_distance;
	}
};


#endif /* SPLINES_H_ */
