#pragma once
#undef NDEBUG
#define _CRT_SECURE_NO_WARNINGS
#include <cassert>
#include <math.h>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>

/*
double half_division(double left_edge, double right_edge, double epsilon=1e-3);

double half_division_test(double left_edge, double right_edge, double epsilon = 1e-3);

double simple_iteration(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c=0.011930646237607);

double simple_iteration_t(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c = 0.0018327550019109917);

double simple_iteration_test(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c = 0.0018327550019109917);

void line_sum(std::vector<std::vector<double>>& mat, size_t src_index, size_t dst_index, double k, bool transpose=false);

std::vector<std::vector<double>>* get_submat(const std::vector<std::vector<double>>& mat, size_t s_min, size_t c_min, size_t s_max, size_t c_max);
*/

// Реализация разложения Холецкого
// Метод принимает симметрическую, и положительно определённую матрицу
// Возвращает указатель на нижнетреугольную матрицу, полученную при разложении
std::vector<std::vector<double>>* cholesky(std::vector<std::vector<double>>& mat);

// Правило суммирования в разложении для элементов на диагонали
// Метод принимает исходную матрицу t_mat, известный элемент a, индекс линии line
// Возвращает значение диагонального элемента в разложении
double ii_element(std::vector<std::vector<double>>& t_mat, double a, size_t line); 

// Правило суммирования в разложении для элементов ниже диагонали
// Метод принимает исходную матрицу t_mat, известный элемент a, индекс линии line, индекс колонки col
// Возвращает значение поддиагонального элемента в разложении
double ij_element(std::vector<std::vector<double>>& t_mat, double a, size_t line, size_t col);

// Основной метод, компанующий все остальные
void lab_2();

// Метод для перебора числа обусловленности и построения норм ошибок и невязки, пишет результаты в файл out.csv
// Принимает начало и конец отрезка [min, max] для линейного распределения значений на главной диагонали
// n - размерность матрицы, stride - шаг числа обусловленности
void dependency(double min, double max, size_t n, double stride = 1);

// Метод генерации симметричной положительно определённой матрицы
// n - размерность матрицы, mod - модуль случайных целых чисел, [min, max] - отрезок для линспейса на главной диагонали
std::vector<std::vector<double>>* generate_matrix(
	size_t n,
	size_t mod = 10,
	double min = 0.01,
	double max = 1,
	Eigen::MatrixXd* diag = nullptr
);

// Печать матрицы и столбца, если есть
// mat - матрица для вывода, b - вектор для вывода(опционален)
void print_mat(const std::vector<std::vector<double>>& mat, const std::vector<double>& b = std::vector<double>());

// Выбор случайного целого числа в зависимости от делителя числа, и сдвига в сторону положительных чисел
// mod - модуль случайных целых чисел, positive - сдвиг в отрицательную сторону
int int_random(int mod, bool positive = false);

// Критерий Сильвестра для матрицы mat
bool sylvester_criterion(const std::vector<std::vector<double>>& mat);

// Вычисление детерминанта
double det(const std::vector<std::vector<double>>& mat);

// Метод вернёт транспонированную матрицу
// mat - транспонируемая матрица
std::vector<std::vector<double>>* transpose(const std::vector<std::vector<double>>& mat);

// Генерация случайного вектора из целых чисел, метод вернёт указатель
// len - длина вектора, mod и positive прокинуты для int_random()
std::vector<double>* generate_vector(size_t len, int mod, bool positive=false);

// Решение треугольной системы методом Гаусса SS'x=b, mat = S, v = b
// Возвращает указатель на вектор - корень. mat - решаемая матрица системы, v - вектор системы.
std::vector<double>* solve_system(const std::vector<std::vector<double>>& mat, const std::vector<double>& v);

// Методы перехода для векторов и матриц между std::vector и Eigen
Eigen::VectorXd vector_to_eigen(const std::vector<double>& v);

std::vector<double>* eigen_to_vector(Eigen::VectorXd v);

Eigen::MatrixXd matrix_to_eigen(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>>* eigen_to_matrix(Eigen::MatrixXd mat);

// Точка входа для 3 лабораторной
void lab_3();

// Метод для перебора ε и логирования бесконечных ||x - x*||, ||Ax - b||, числа итераций.
// Пишет результаты в файл out.csv
// Принимает начало и конец диапазона [min_eps, max_eps] для итераций по нему с шагом min_eps
// n - размерность матрицы, beta - параметр релаксации
void dependency_3(double min_eps, double max_eps, size_t n, double beta = 1.5);

// Метод релаксации
// mat -  квадратная, симметричная, положительно определённая матрица
// b - вектор свободных значений в СЛАУ
// beta - коэффициент релаксации из отрезка (0, 2)
// eps - точность решения
// Возвращает пару элементов: вектор и число итераций
std::pair<std::vector<double>*, size_t>* SOR(
	const std::vector<std::vector<double>>& mat,
	const std::vector<double>& b,
	double beta,
	double eps
);

// Вычитание двух векторов
std::vector<double> v_sub(const std::vector<double>& x0, const std::vector<double>& x1);

// Вычисление бесконечной нормы в матрице
double matrix_inf_norm(const std::vector<std::vector<double>>& mat);

// Вычисление бесконечной нормы вектора
double matrix_inf_norm(const std::vector<double>& v);

// Вернёт максимальный элемент вектора
double max(const std::vector<double>& v);

// Точка входа для 4 лабораторной
void lab_4();

// Метод для перебора ε и логирования ||x - x*||, ||Ax - λx||, числа итераций.
// Пишет результаты в файл out.csv
// Принимает начало и конец диапазона [min_eps, max_eps] для итераций по нему с шагом min_eps
// n - размерность матрицы
void dependency_4(double min_eps, double max_eps, size_t n);

// Вычисление произведения матриц
std::vector<std::vector<double>>* matrix_mul(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b);

// Метод реализующий обнуление оптимального элемента a_ij
// Возвращает указатель на матрицу поворота
// Изменяет mat в ходе обнуления a_ij
std::vector<std::vector<double>>* Jacobi_rotation(std::vector<std::vector<double>>& mat);

// Метод вращений Якоби
// mat -  квадратная, симметричная матрица
// eps - требуемая точность решения
// Возвращает пару элементов: диагональную матрицу собственных значений и вектор из матриц вращений
std::pair<std::vector<std::vector<double>>*, std::vector<std::vector<std::vector<double>>*>*>* Jacobi(
	const std::vector<std::vector<double>> mat,
	double eps
);

// Вычисление суммы квадратов элементов вне диагонали
double ndss(const std::vector<std::vector<double>>& mat);

// Выбор оптимального элемента для обнуления
// mat - матрица для выбора
// Метод вернёт пару индексов (i строки, j столбца)
std::pair<size_t, size_t> select_opt(const std::vector<std::vector<double>>& mat);
