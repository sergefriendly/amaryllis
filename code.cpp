#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>


struct Point {
    double t, x, y;
    Point(double t, double x, double y) : t(t), x(x), y(y) {}
};

struct Segment {
    double a, b, c, d, t_start, t_end;
};

class Curve {
private:
    std::vector<Point> sorted_points;
    std::vector<Segment> segments;

public:
    Curve(const std::vector<Point>& points) {
        if (points.size() < 2) {
            throw std::invalid_argument("At least two points are required.");
        }

        sorted_points = points;
        std::sort(sorted_points.begin(), sorted_points.end(), [](const Point& p1, const Point& p2) {
            return p1.t < p2.t;
        });

        size_t n = sorted_points.size() - 1;
        std::vector<double> h(n);
        for (size_t i = 0; i < n; ++i) {
            h[i] = sorted_points[i + 1].t - sorted_points[i].t;
            if (h[i] <= 0) {
                throw std::invalid_argument("Points must have strictly increasing t values after sorting.");
            }
        }

        std::vector<std::vector<double>> A(n + 1, std::vector<double>(n + 1, 0));
        std::vector<double> b(n + 1, 0);

        A[0][0] = 1;
        A[n][n] = 1;

        for (size_t i = 1; i < n; ++i) {
            A[i][i - 1] = h[i - 1];
            A[i][i] = 2 * (h[i - 1] + h[i]);
            A[i][i + 1] = h[i];
            b[i] = 3 * ((sorted_points[i + 1].x - sorted_points[i].x) / h[i] -
                        (sorted_points[i].x - sorted_points[i - 1].x) / h[i - 1]);
        }

        std::vector<double> c = solveTridiagonal(A, b);

        segments.resize(n);
        for (size_t i = 0; i < n; ++i) {
            double a = sorted_points[i].x;
            double b = (sorted_points[i + 1].x - sorted_points[i].x) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
            double d = (c[i + 1] - c[i]) / (3 * h[i]);
            segments[i] = {a, b, c[i], d, sorted_points[i].t, sorted_points[i + 1].t};
        }
    }

    double interpolateX(double t) const {
        if (t < sorted_points.front().t || t > sorted_points.back().t) {
            throw std::out_of_range("t is out of the interpolation range");
        }

        auto it = std::lower_bound(sorted_points.begin(), sorted_points.end(), t,
                                   [](const Point& p, double val) { return p.t < val; });
        size_t i = std::distance(sorted_points.begin(), it) - 1;
        const Segment& seg = segments[i];
        double dt = t - seg.t_start;
        return seg.a + seg.b * dt + seg.c * dt * dt + seg.d * dt * dt * dt;
    }

    double interpolateY(double x) const {
        if (x < sorted_points.front().x || x > sorted_points.back().x) {
            throw std::out_of_range("x is out of the interpolation range");
        }

        // Find the segment containing x
        auto it = std::lower_bound(sorted_points.begin(), sorted_points.end(), x,
                                   [](const Point& p, double val) { return p.x < val; });
        size_t i = std::distance(sorted_points.begin(), it) - 1;

        // Linear interpolation for y
        const Point& p1 = sorted_points[i];
        const Point& p2 = sorted_points[i + 1];
        double t = p1.t + (x - p1.x) * (p2.t - p1.t) / (p2.x - p1.x);
        return p1.y + (p2.y - p1.y) * (t - p1.t) / (p2.t - p1.t);
    }

    static std::vector<Point> findIntersections(const Curve& curve1, const Curve& curve2, double epsilon = 1e-6) {
        std::vector<Point> intersections;
        double t_min = std::max(curve1.sorted_points.front().t, curve2.sorted_points.front().t);
        double t_max = std::min(curve1.sorted_points.back().t, curve2.sorted_points.back().t);

        for (double t = t_min; t <= t_max; t += epsilon) {
            double x1 = curve1.interpolateX(t);
            double y1 = curve1.interpolateY(x1);
            double x2 = curve2.interpolateX(t);
            double y2 = curve2.interpolateY(x2);

            if (std::abs(x1 - x2) < epsilon && std::abs(y1 - y2) < epsilon) {
                intersections.push_back({t, x1, y1});
            }
        }
        return intersections;
    }

    static double findMinimalDistance(const Curve& curve1, const Curve& curve2, double epsilon = 1e-6) {
        double min_distance = std::numeric_limits<double>::max();
        double t_min = std::max(curve1.sorted_points.front().t, curve2.sorted_points.front().t);
        double t_max = std::min(curve1.sorted_points.back().t, curve2.sorted_points.back().t);

        for (double t = t_min; t <= t_max; t += epsilon) {
            double x1 = curve1.interpolateX(t);
            double y1 = curve1.interpolateY(x1);
            double x2 = curve2.interpolateX(t);
            double y2 = curve2.interpolateY(x2);

            double distance = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
            min_distance = std::min(min_distance, distance);
        }
        return min_distance;
    }

private:
    std::vector<double> solveTridiagonal(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
        size_t n = b.size();
        std::vector<double> x(n);
        std::vector<double> p(n);
        std::vector<double> q(n);

        p[0] = A[0][1] / A[0][0];
        q[0] = b[0] / A[0][0];

        for (size_t i = 1; i < n - 1; ++i) {
            p[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * p[i - 1]);
            q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] - A[i][i - 1] * p[i - 1]);
        }

        x[n - 1] = (b[n - 1] - A[n - 1][n - 2] * q[n - 2]) / (A[n - 1][n - 1] - A[n - 1][n - 2] * p[n - 2]);

        for (int i = n - 2; i >= 0; --i) {
            x[i] = q[i] - p[i] * x[i + 1];
        }

        return x;
    }

    static void parsePoints(const std::string& line, std::vector<Point>& points) {
        std::istringstream iss(line);
        std::string token;
        double t = 0.0; // Initialize t to 0 and increment for each point

        while (std::getline(iss, token, ',')) {
            std::istringstream tokenStream(token);
            double x, y;
            char dummy; // To consume the space between numbers

            tokenStream >> x >> dummy >> y;
            points.emplace_back(t, x, y);
            t++;
        }
    }


public:
    static void yahoo() {
    	std::string filename = "/home/me/_/homes/eclipse-cpp-2024-06-22-workspace/ASCONA/src/points.txt";
		  std::ifstream file(filename);
		  std::string line;

  		if (file.is_open()) {
  			std::vector<Point> points1, points2;
  			while (std::getline(file, line)) {
  				std::size_t semicolonPos = line.find(';');
  				if (semicolonPos != std::string::npos) {
  					std::string points1Str = line.substr(0, semicolonPos);
  					std::string points2Str = line.substr(semicolonPos + 1);
  
  					parsePoints(points1Str, points1);
  					parsePoints(points2Str, points2);
  				}
  			}
  
  			// For demonstration: print the points from both vectors
  			for (const auto& point : points1) {
  				std::cout << "Point1: " << point.t << " " << point.x << " " << point.y << std::endl;
  			}
  			for (const auto& point : points2) {
  				std::cout << "Point2: " << point.t << " " << point.x << " " << point.y << std::endl;
  			}
  
  			file.close();
  		} else {
  			std::cerr << "Unable to open file" << std::endl;
  		}
    }
};

int test() {
  
  return 0;
}
int __main__() {
    std::vector<Point> points1 = {{0, 0, 0}, {1, 1, 2}, {2, 4, 3}, {3, 9, 5}};
    std::vector<Point> points2 = {{0, 0, 5}, {1, 1, 3}, {2, 4, 2}, {3, 9, 0}};

    Curve curve1(points1);
    Curve curve2(points2);

    std::cout << "Intersections:\n";
    auto intersections = Curve::findIntersections(curve1, curve2);
    for (const auto& p : intersections) {
        std::cout << "t: " << p.t << ", x: " << p.x << ", y: " << p.y << "\n";
    }

    if (intersections.empty()) {
        double min_distance = Curve::findMinimalDistance(curve1, curve2);
        std::cout << "Minimal distance between curves: " << min_distance << "\n";
    }

    return 0;
}
