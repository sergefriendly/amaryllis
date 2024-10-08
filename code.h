/**
 * Структура проекта по решению тестового задания:
 * 1. Структура Point (Точка).
 * 2. Структура Segment (Сегмент).
 * 3. Класс Curve (Кривая) — абстрактный класс, содержит виртуальные методы,
 *	содержит реализованными общие методы.
 * 4. Класс CubicCurve (Кубическая кривая) — класс, реализованный на базе класса Curve,
 *	содержит реализованными конкретные методы.
 * 5. Класс HermitCurve (Кривая Эрмита) — класс, присутствующий формально,
 *	по своей сути реализует иной, отличный от кибического, алгоритм аппроксимации.
 * 6. Класс CurveFactory (Фабрика кривых) — класс, реализующий шаблон "Фабричный метод"
 *	для создания специализированных кривых, т.е. кривых, в которых реализованы
 *	специализированные алгоритмы аппроксимации.
 * 7. Функция __main__ — наличествует для реализации демонстрационного кода — кода
 *	на базе сущностей которого решается тестовое задание.
 *
 *
 * Векторы (std::vectors<...>) от структур Point и Segment находятся с классом Curve в отношении композиции.
 * Классы CubicCurve и HermitCurve находятся с Curve в отношении наследования.
 * Класс CurveFactory находится с классами CubicCurve и HermitCurve в отношении ассоциации.
 *
 * Проект реализован на шаблоне "Фабричный метод".
 * В проекте используется C++17 диалект.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>


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
 * последовательности точек, в том числе, если,
 * собственно, функция y(x) будет многозначна.
 */
struct Point {
	double t, x, y; // Параметр и две координаты для плоской кривой
	Point(double t, double x, double y) :
			t(t), x(x), y(y) {
	} // Конструктор

	friend std::ostream& operator<<(std::ostream &os, const Point &point) { // Принтер
		os << "( " << point.t << ": " << point.x << ", " << point.y << " )";
		return os;
	}

	/*
	 * Метод, считывающий строку в текстовом файле
	 */
	static std::vector<Point> parsePoints(const std::string &line) {

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

			points.push_back(Point { t, x, y }); // Затем считанное добавляется в вектор
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
	 * в синтаксическом представлении.
	 */
	static void fillVectorsFromTextFile(std::string file_name,
			std::vector<Point> &points1, std::vector<Point> &points2,
			bool putsInCommandLine = false) {

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
				points1 = parsePoints(line);
				isFirstSet = false;
			} else {
				points2 = parsePoints(line);
			}
		}

		file.close();

		/* Следующий блок кода печатает в командной строке вектор
		 * точек при необходимости в зависимости от параметра putsInCommandLine:
		 */
		if (putsInCommandLine) {
			std::cout << "Последовательность точек 1:" << std::endl;
			for (const auto &point : points1) {
				std::cout << point << std::endl;
			}

			std::cout << "Последовательность точек 2:" << std::endl;
			for (const auto &point : points2) {
				std::cout << point << std::endl;
			}
		}
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
	double a = 0, b = 0, c = 0, d = 0, t_min = 0, t_max = 1; // Параметры
	Segment() {
	} // Тривиальный конструктор решено оставить для теста
	Segment(double a, double b, double c, double d, double t_min, double t_max) :
			a(a), b(b), c(c), d(d), t_min(t_min), t_max(t_max) {} // Конструктор

	friend std::ostream& operator<<(std::ostream &os, const Segment &segment) { // Принтер
		os << "( a = " << segment.a << ", b = " << segment.b << ", c = "
				<< segment.c << ", d = " << segment.d << ", [" << segment.t_min
				<< ", " << segment.t_max << "] )";
		return os;
	}
};

/*
 * Класс кривая — абстрактный класс, в котором методы
 * buildOn, interpXbyT, interpYbyT, interpYbyX обозначаются как
 * виртуальные для того, чтобы реализовать их в конкретном дочернем классе,
 * т.к. их содержимое может меняться в зависимости от метода интерполяции сплайна:
 * кубический ли, Эрмита ли и т.д.
 */
class Curve {
public:
	/*
 	 * Статическая переменная, служит для обозначения точности,
   	 * хотя может быть определена конкретно к экземпляру —
     	 * зависит от дальнейших потребностей применения задачи.
   	 */
	static double accuracy_tolerance;

	std::vector<Point> points;
	std::vector<Segment> segments;

	virtual ~Curve() = default;

	virtual void buildOn(const std::vector<Point> &_points) = 0;
	virtual double interpXbyT(double t) const = 0;
	virtual double interpYbyT(double t) const = 0;
	virtual double interpYbyX(double x) const = 0;

	size_t findSegment(double t) const {
		/*
		 * Бинарный поиск. Цель которого найти такую точку, параметр t
		 * которой был бы ближе к значению value
		 */
		auto it = std::lower_bound(points.begin(), points.end(), t,
				[](const Point &p, double value) {
					return p.t < value;
				});

		// Находим и возвращаем индекс точки it на единицу мешьший
		return std::distance(points.begin(), it) - 1;
	}

	static std::vector<Point> findIntersections(const Curve * const curve1,
			const Curve * const curve2, double epsilon = accuracy_tolerance) {
		/*
  		 * В этот вектор будем собирать точки пересечения с заданной точностью epsilon
     		 */
		std::vector<Point> intersections;

		/*
  		 * Вычисляем область для последующего поиска общих точек пересечения
     		 */
		double t_min = std::max(curve1->points.front().t, curve2->points.front().t);
		double t_max = std::min(curve1->points.back().t,  curve2->points.back().t);

		/*
  		 * Интерполируем необходимые занчения, пробегая повсем точкам области, находя при этом
     		 * абсолютные занчения, по которым в всою очередь и определяется является ли пара точкой пересечения или нет.
     		 */
		for (double t = t_min; t <= t_max; t += epsilon) {
			double x1 = curve1->interpXbyT(t);
			double y1 = curve1->interpYbyX(x1);
			double x2 = curve2->interpXbyT(t);
			double y2 = curve2->interpYbyX(x2);

			if (std::abs(x1 - x2) < epsilon && std::abs(y1 - y2) < epsilon) {
				intersections.push_back( { t, x1, y1 });
			}
		}
		return intersections;
	}

	static double findMinimalDistance(const Curve * const curve1, const Curve * const curve2,
			double epsilon = accuracy_tolerance) {
		double min_distance = std::numeric_limits<double>::max();
		double t_min = std::max(curve1->points.front().t,
				curve2->points.front().t);
		double t_max = std::min(curve1->points.back().t, curve2->points.back().t);

		for (double t = t_min; t <= t_max; t += epsilon) {
			double x1 = curve1->interpXbyT(t);
			double y1 = curve1->interpYbyX(x1);
			double x2 = curve2->interpXbyT(t);
			double y2 = curve2->interpYbyX(x2);

			double distance = std::sqrt(
					std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
			min_distance = std::min(min_distance, distance);
		}
		return min_distance;
	}
};

class CubicCurve: public Curve {
public:
	/*
	 * Точки и сегменты. Если точек N, то сегментов N - 1,
	 * поскольку один сегмент кривой, восстановленный посредством
	 * кубического сплайна, смежен с двумя точками - таким образом
	 * есть набор точек, а между ними восстанавливаются дуги коих
	 * и будет всегда меньше на один. Поскольку введёное количество
	 * точек пользователем может быть условно говоря любым, то
	 * в даммон случае лучше выбрать для хранения такую структуру данных
	 * как std::vector и следовательно для хранения сегментов тоже.
	 */

	void buildOn(const std::vector<Point> &_points) override {
		if (_points.size() < 2) {
			throw std::invalid_argument("Необходимы как минимум две точки.");
		}

		points = _points;

		// Сортировка точек в порядке возрастания параметра t.
		std::sort(points.begin(), points.end(),
				[](const Point &p1, const Point &p2) {
					return p1.t < p2.t;
				});

		size_t N_pnt = points.size();
		size_t N_seg = points.size() - 1;

		/*
		 * Вектор для храния длин покрывающих отрезок [t_min, t_max] интервалов,
		 * получившихся в результате ввода координат (x, y) пользователем
		 */
		std::vector<double> h(N_seg);

		// Заполнение данного вектора
		for (size_t i = 0; i < N_seg; i++) {
			h[i] = points[i + 1].t - points[i].t;
			if (h[i] <= 0) {
				throw std::invalid_argument(
						"После сортировки точки должны быть в возрастании параметра t.");
			}
		}

		/*
		 * Далее моделируется квадратная матрица A и вектор-столбец b
		 * которые в последствии будут определяющими компонентами
		 * для матричного уравнения Ax = b = 0, которе в последствии предстоит
		 * решить для получения коэфициентов a, b, c, d, определяющих
		 * собственно говоря кубический сплайн.
		 *
		 * Примечание: следует иметь ввиду, что вектор-столбец b и
		 * коэффициент сплайна b разные сущности.
		 */
		std::vector<std::vector<double>> A(N_pnt,
				std::vector<double>(N_pnt, 0));
		std::vector<double> b(N_pnt, 0);

		A[0][0] = 1;
		A[N_seg][N_seg] = 1;

		// Заполнение матрицы в соответствии с методом
		for (size_t i = 1; i < N_seg; i++) {
			A[i][i - 1] = h[i - 1];
			A[i][i] = 2 * (h[i - 1] + h[i]);
			A[i][i + 1] = h[i];
			b[i] = 3 * ((points[i + 1].x - points[i].x) / h[i] - (points[i].x - points[i - 1].x) / h[i - 1]);
		}

		// Решаем матричное уравнение Аx = b = 0 методом прогонки
		std::vector<double> c = solveTridiagonal(A, b);

		segments.resize(N_seg);

		/*
		 * Заполняем высчитанными коэффициентами параметрами,
		 * определяющими каждый сегмент.
		 */
		for (size_t i = 0; i < N_seg; i++) {
			double a = points[i].x;
			double b = (points[i + 1].x - points[i].x) / h[i]
					- h[i] * (c[i + 1] + 2 * c[i]) / 3;
			double d = (c[i + 1] - c[i]) / (3 * h[i]);
			segments[i] = { a, b, c[i], d, points[i].t, points[i + 1].t };
		}
	}

	// Метод прогонки
	std::vector<double> solveTridiagonal(
			const std::vector<std::vector<double>> &A,
			const std::vector<double> &b) {
		size_t n = b.size();
		std::vector<double> x(n), p(n), q(n);

		p[0] = A[0][1] / A[0][0];
		q[0] = b[0] / A[0][0];

		for (size_t i = 1; i < n - 1; i++) {
			p[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * p[i - 1]);
			q[i] = (b[i] - A[i][i - 1] * q[i - 1])
					/ (A[i][i] - A[i][i - 1] * p[i - 1]);
		}

		x[n - 1] = (b[n - 1] - A[n - 1][n - 2] * q[n - 2])
				/ (A[n - 1][n - 1] - A[n - 1][n - 2] * p[n - 2]);

		for (int i = n - 2; i >= 0; i--) {
			x[i] = q[i] - p[i] * x[i + 1];
		}

		return x;
	}

	double interpXbyT(double t) const override {
		/*
		 * Если параметр берётся не из зоны покртия, которую определяют
		 * точки введеные пользователем и для которой определился алгоритм,
		 * то в этом случае естественно нужно выкинуть ошибку.
		 */
		if (t < points.front().t || t > points.back().t) {
			throw std::out_of_range("t не в интерполяционном интервале.");
		}

		/*
		 * Бинарный поиск. Цель которого найти такую точку, параметр t
		 * которой был бы ближе к значению val
		 */
		auto it = std::lower_bound(points.begin(), points.end(), t,
				[](const Point &p, double val) {
					return p.t < val;
				});

		// Находим индекс точки it на единицу мешьший
		size_t i = std::distance(points.begin(), it) - 1;

		const Segment &seg = segments[i]; // Берём этот сегмент
		double dt = t - seg.t_min; // Определяем расстояние от t до t_min

		/*
		 * В выбранном сегменте возвращаем значение из аппроксимированной
		 * кубическим сплайном дуги
		 */
		return seg.a + seg.b * dt + seg.c * dt * dt + seg.d * dt * dt * dt;
	}

	double interpYbyT(double t) const override {
		if (t < points.front().t || t > points.back().t) {
			throw std::out_of_range("Interpolation value t is out of range.");
		}

		size_t i = findSegment(t);
		double y0 = points[i].y;
		double y1 = points[i + 1].y;
		double t0 = points[i].t;
		double t1 = points[i + 1].t;

		// Линейная интерполяция для Y
		return y0 + (y1 - y0) * (t - t0) / (t1 - t0);
	}

	double interpYbyX(double x) const override {
		if (x < points.front().x || x > points.back().x) {
			throw std::out_of_range("x не в интерполяционном интервале.");
		}

		auto it = std::lower_bound(points.begin(), points.end(), x,
				[](const Point &p, double val) {
					return p.x < val;
				});
		size_t i = std::distance(points.begin(), it) - 1;

		const Point &p1 = points[i];
		const Point &p2 = points[i + 1];
		double t = p1.t + (x - p1.x) * (p2.t - p1.t) / (p2.x - p1.x);
		return p1.y + (p2.y - p1.y) * (t - p1.t) / (p2.t - p1.t);
	}

};

/*
 * Присутствует формально
 */
class HermiteCurve: public Curve {
	void buildOn(const std::vector<Point> &_points) override {
		// a code
	}

	double interpXbyT(double t) const override {
		return 0.0;
	}

	double interpYbyT(double t) const override {
		return 0.0;
	}

	double interpYbyX(double x) const override {
		return 0.0;
	}
};

class CurveFactory {
public:
	enum CurveType {
		CUBIC_SPLINE, HERMITE_SPLINE
	};

	static Curve* createCurve(CurveType type) {
		switch (type) {
		case CUBIC_SPLINE:
			return new CubicCurve();
		case HERMITE_SPLINE:
			return new HermiteCurve();
		default:
			return nullptr;
		}
	}
};

double Curve::accuracy_tolerance = 1e-5;

int __main__(int argc, char *argv[]) {


	std::cout << std::fixed;
	std::cout << std::setprecision(9); // необходим для указания количества показываемых знаков после запятой для числа

	std::vector<Point> points1, points2; // определяем два набора векторов точек, т.к. производим операции над двумя кривыми
	Point::fillVectorsFromTextFile(argv[1], points1, points2, true); // Имя файла передаётся как параметр командной строки через argv[1]

	/*
	 * Считываем второй параметр командной строки, в который передаётся точность
	 */
	char* end;
	Curve::accuracy_tolerance = std::strtod(argv[2], &end);
	if (*end != '\0') {
		std::cerr << "Недопустимое число: " << argv[2] << std::endl;
		return 1;
	}

	/*
	 * Создаём две кривые на базе считанных координат точек, используя фабричный метод,
  	 * который создаёт объект динамически.
	 */
	auto curve1 = CurveFactory::createCurve(CurveFactory::CUBIC_SPLINE);
	auto curve2 = CurveFactory::createCurve(CurveFactory::CUBIC_SPLINE);
	curve1->buildOn(points1);
	curve2->buildOn(points2);

	/*
	 * Для примера проводится интерполирование для первой и второй кривой значений, взятых из их областей определния
	 */
	std::cout << "val |-> interp" << std::endl;
	double val1 = 1.5;
	double interp1 = curve1->interpYbyX(val1);
	std::cout << val1 << " |-> " << interp1 << std::endl;

	double val2 = 2.5;
	double interp2 = curve2->interpYbyX(val2);
	std::cout << val2 << " |-> " << interp2 << std::endl;

	// Затем, вычислям прересечения
	std::vector<Point> interpPoints = Curve::findIntersections(curve1, curve2);
	std::cout << "Пересечения: " << std::endl;
	for (const auto& p : interpPoints) {
		std::cout << p << std::endl;
	}
	if (interpPoints.empty()) {
		std::cout << "Пересечений нет." << std::endl;
	}

	// И затем, находим минимальную дистанцию между кривыми.
	double min_dist = Curve::findMinimalDistance(curve1, curve2);
	std::cout << "Минимальная дистанция: " << min_dist << std::endl;
	return 0;
}
