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

namespace spline {

/*
 * Кривую можно разбить на сегметы,
 * пара последовательно идущих точек
 * определяет этот сегмент, который
 * аппроксимируется кубическим спланом.
 * В данном случае определяется модель точки:
 * две координаты x и y и параметр t.
 * Данное определение более функциональное, чем
 * определение без параметра. Параметрически
 * заданная функция или правильнее сказать
 * вектор-функция r(t) = {x(t), y(t)}
 * с точки зрения функции y(x)
 * может быть не всегда однозначней, а это
 * в свою очередь ограничивает пределы применимости
 * такого рода функций и поэтому реализация
 * алгоритма для вектор-функции будет перспективнее,
 * ибо алгоритм в таком случае будет годен для любой
 * последовательности точек, в том числе, если
 * собственно функция y(x) будет неоднозначна
 * как таковая.
 */
struct Point {
	double t, x, y; // Параметр и две координаты для плоской кривой
	Point(double t, double x, double y) : t(t), x(x), y(y) {} // Конструктор

	friend std::ostream& operator<<(std::ostream& os, const Point& point) { // Принтер
		os << "( " << point.t << ": " << point.x << ", " << point.y << " )";
		return os;
	}
};

/*
 * Сегмент кривой будет просчитываться кубическим сплайном.
 * Для определения сегмента понадобится 5 параметров.
 * Собственно параметры a, b, c, d для определения кубического
 * многочлена a + b*t + c*t^2 + d*t^3 и интервал множества
 * первых элементов [t_min, t_max] на котором будет строиться
 * многочлен собственно который будет являться сплайном.
 */
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
	/*
	 *
	 */
	std::vector<Point> points; // Точки, если их N + 1,
	std::vector<Segment> segments; // то сегментов будет меньше на один, т.е. N

	Curve() {} // Тривиальный конструктор решено оставить для теста

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

	/*
	 *
	 */
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

	/*
	 * Метод, считывающий строку в текстовом файле
	 */
	static std::vector<Point> parsePoints(const std::string& line) {

		std::vector<Point> points;
		/*
		 * Оставлен как буферный,
		 * но при необходимости может быть исключён
		 * посредством применения ссылок и указателей,
		 * дабы такая структура уже обрабатывается
		 * в конструкторе Curve.
		 */

	    std::istringstream stream(line); // Необходим для считавания строки текста из файла
	    double x, y;
	    char separator;
	    double t = 0.0;
	    /*
	     * Начальное значение параметра t выбранно 0.0,
	     * но в принципе может быть любым: вся важность
	     * заключена в параметрах (координатах) (x, y) которые определяют
	     * своего рода "скелет" кривой будучи достроенной
	     * кубическими сплайнами. Параметр же t определяет
	     * порядок следования координат (x, y), поэтому решено
	     * нумеровать пары (x,y) с t = 0.0 и далее инкрементируя
	     * на единицу.
	     */

	    while (stream >> x >> y) {
	    	/*
	    	 * Считывается два последующих значения, которые,
	    	 * как предполагется, идут всегда один за другим
	    	 * и затем символ запятой.
	    	 */

	        points.push_back(Point{t, x, y}); // Затем считанное добавляется в вектор
	        t += 1.0;
	        stream >> separator; // Запятая или точка с запятой проскакиваются
	    }

	    return points;
	}

	/*
	 * Метод выбран как статический, в силу того, что задача ограничевается только лишь
	 * соим решением и её дальшейшая реализация не предполагается в данном
	 * конкретном случае. И всё же при необходимости метод может быть реорганизован, как
	 * обычный метод для экземпляра класса. Также он может быть организован в специальном
	 * классе, например, CurveOperations при необходимиости. В данном случае
	 * методов для операций с кривыми не там много, поэтому остановимся на статических
	 * методах, поскольку в такой организации у метода есть симметрия, которая более гармонична
	 * синтаксическом представлении.
	 */
	static void fillVectorsFromTextFile(std::string file_name,
			std::vector<Point>*const points1, std::vector<Point>*const points2,
			bool putsInCommandLine=false) {

		std::ifstream file(file_name);
		if (!file.is_open()) {
			std::cerr << "Ошибка отрытия файла" << std::endl;
		}

		std::string line;
		bool isFirstSet = true;

		/*
		 * В этом блоке кода реализуется простейший конечынй автомат
		 * который меняет своё состояние по символу `;`
		 */
		while (getline(file, line, ';')) {
			/*
			 * Из файла file считывается текст до символа `;`
			 * в строку line.
			 */

			if (isFirstSet) {
				*points1 = parsePoints(line);
				isFirstSet = false;
			} else {
				*points2 = parsePoints(line);
			}
		}

		file.close();

		/* Следующий блок кода печатает в командной строке вектор
		 * точек при необходимости в зависимости от параметра putsInCommandLine:
		 */
		if (putsInCommandLine) {
			std::cout << "Последовательность точек 1:" << std::endl;
			for (const auto& point : *points1) {
				std::cout << point << std::endl;
			}

			std::cout << "Последовательность точек 2:" << std::endl;
			for (const auto& point : *points2) {
				std::cout << point << std::endl;
			}
		}
	}
};

}

/*
 * std::vector<spline::Point> points1, points2;
	spline::Curve::fillVectorsFromTextFile("/home/me/_/homes/eclipse-cpp-2024-06-22-workspace/ASCONA/src/lala.txt",
			&points1, &points2, true);
 */

#endif /* SPLINES_H_ */
