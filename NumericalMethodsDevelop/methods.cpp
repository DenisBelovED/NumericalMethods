#include "methods.h"
#include <iostream>
#include <Eigen/Dense>

/*
double half_division(double left_edge, double right_edge, double epsilon)
{
	double (*f)(double) = [](double x) { return pow(x, 5) + pow(x, 2) - 5; };
	assert (left_edge < right_edge);
	assert (f(left_edge) * f(right_edge) < 0);
	double c;
	while (right_edge - left_edge > epsilon)
	{
		c = (right_edge + left_edge) / 2;
		if (f(right_edge) * f(c) < 0)
			left_edge = c;
		else
			right_edge = c;
	}
	return c;
}

double half_division_test(double left_edge, double right_edge, double epsilon)
{
	double (*f)(double) = [](double x) { return pow(x, 2) - 1; };
	assert(left_edge < right_edge);
	assert(f(left_edge) * f(right_edge) < 0);
	double c;
	while (right_edge - left_edge > epsilon)
	{
		c = (right_edge + left_edge) / 2;
		if (f(right_edge) * f(c) <= 0)
			left_edge = c;
		else
			right_edge = c;
	}
	return c;
}

double simple_iteration(double left_edge, double right_edge, double x_0, double epsilon, double c)
{
	double (*f)(double) = [](double x) { return pow(x, 5) + pow(x, 2) - 5; };
	double x_i = x_0, x_j = x_0 - c * f(x_0);
	while (fabs(x_j - x_i) > epsilon)
	{
		x_i = x_j;
		x_j = x_i - c * f(x_i);
	}
	return x_j;
}

double simple_iteration_t(double left_edge, double right_edge, double x_0, double epsilon, double c)
{
	double (*f)(double) = [](double x) { return pow(tan(x), 3) - x + 1; };
	double x_i = x_0, x_j = x_0 - c * f(x_0);
	while (fabs(x_j - x_i) > epsilon)
	{
		x_i = x_j;
		x_j = x_i - c * f(x_i);
		std::cout << x_j << std::endl;
	}
	return x_j;
}

double simple_iteration_test(double left_edge, double right_edge, double x_0, double epsilon, double c)
{
	double (*f)(double) = [](double x) { return pow(x, 2) - 1; };
	double x_i = x_0, x_j = x_0 - c * f(x_0);
	while (fabs(x_j - x_i) > epsilon)
	{
		x_i = x_j;
		x_j = x_i - c * f(x_i);
		std::cout << x_j << std::endl;
	}
	return x_j;
}
*/

std::vector<std::vector<double>>* cholesky(std::vector<std::vector<double>>& mat)
{
	auto l_mat = new std::vector<std::vector<double>>(mat.size(), std::vector<double>(mat.size(), 0));
	for (size_t i = 0; i < mat.size(); i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			if (i == j)
				(*l_mat)[i][j] = ii_element(*l_mat, mat[i][j], i);
			else if (j < i)
				(*l_mat)[i][j] = ij_element(*l_mat, mat[i][j], i, j);
		}
	}
	return l_mat;
}

double ii_element(std::vector<std::vector<double>>& t_mat, double a, size_t line)
{
	double sum = 0;
	for (int i = 0; i < line; i++)
		sum += t_mat[line][i] * t_mat[line][i];
	return std::sqrt(a - sum);
}

double ij_element(std::vector<std::vector<double>>& t_mat, double a, size_t line, size_t col)
{
	double sum = 0;
	for (int i = 0; i <= col; i++)
		sum += t_mat[line][i] * t_mat[col][i];
	return (a - sum) / t_mat[col][col];
}

void lab_2()
{
	size_t n;
	double min, max;
	std::cout << "Input n-dim:" << std::endl;
	std::cin >> n;
	std::cout << "Input real D_min and D_max for generate condition number:" << std::endl;
	std::cin >> min >> max;
	std::cout << "Generate matrix A: " << n << "x" << n << std::endl;
	auto mat = generate_matrix(n, 100, min, max);
	print_mat(*mat);
	auto S = cholesky(*mat);
	auto S_t = transpose(*S);
	std::cout << "L" << std::endl;
	print_mat(*S);
	std::cout << "L'" << std::endl;
	print_mat(*S_t);
	delete mat;
	delete S;
	delete S_t;
}

double det(const std::vector<std::vector<double>>& mat)
{
	if (mat.size() == 1)
		return mat[0][0];
	if (mat.size() == 2)
		return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
	double d = 0;
	for (size_t j = 0; j < mat.size(); j++)
	{
		std::vector<std::vector<double>> sub_mat(mat.size() - 1, std::vector<double>());
		for (size_t i = 1; i < mat.size(); i++)
			for (size_t k = 0; k < mat.size(); k++)
				if (k != j)
					sub_mat[i-1].push_back(mat[i][k]);
		d += pow(-1, j) * mat[0][j] * det(sub_mat);
	}
	return d;
}

std::vector<std::vector<double>>* generate_matrix(size_t n, size_t mod, double min, double max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(-1.0, 1.0);

	Eigen::MatrixXd matx = Eigen::MatrixXd::NullaryExpr(n, n, [&]() {return dis(gen); });
	Eigen::MatrixXd Q;
	Eigen::MatrixXd D = Eigen::VectorXd::LinSpaced(n, min, max).asDiagonal();
	matx *= mod;
	while (matx.determinant() == 0)
	{
		matx.setRandom();
		matx *= mod;
	}
	Eigen::HouseholderQR<Eigen::MatrixXd> qr(matx);
	Q = qr.householderQ();
	matx = Q * D * Q.transpose();

	std::cout << "Condition number = " << matx.norm() * matx.reverse().norm() << std::endl;

	auto mat = new std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < n; j++)
			(*mat)[i][j] = matx(i, j);
	return mat;
}

void line_sum(std::vector<std::vector<double>>& mat, size_t src_index, size_t dst_index, double k, bool transpose)
{
	if (transpose)
		for (size_t i = 0; i < mat.size(); i++)
			mat[i][dst_index] += k * mat[i][src_index];
	else
		for (size_t i = 0; i < mat.size(); i++)
			mat[dst_index][i] += k * mat[src_index][i];
}

void print_mat(const std::vector<std::vector<double>>& mat)
{
	for (size_t i = 0; i < mat.size(); i++)
	{
		for (size_t j = 0; j < mat.size(); j++)
			printf("%.8f ", mat[i][j]);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

bool sylvester_criterion(const std::vector<std::vector<double>>& mat)
{
	std::cout << "Sylvester criterion:" << std::endl;
	double d = mat[0][0];
	std::cout << "Minor 1x1 = " << d << std::endl;
	if (d <= 0)
		return false;
	for (size_t i = 1; i < mat.size(); i++)
	{
		std::vector<std::vector<double>> sub_mat(i+1, std::vector<double>(i+1, 0));
		for (size_t g = 0; g < i+1; g++)
			for (size_t k = 0; k < i+1; k++)
				sub_mat[g][k] = mat[g][k];
		d = det(sub_mat);
		std::cout << "Minor " << i+1 << "x" << i+1 << " = " << d << std::endl;
		if (d <= 0)
			return false;
	}
	std::cout << "Succsess!" << std::endl;
	return true;
}

std::vector<std::vector<double>>* get_submat(const std::vector<std::vector<double>>& mat, size_t s_min, size_t c_min, size_t s_max, size_t c_max)
{
	auto sub_mat = new std::vector<std::vector<double>>(s_max - s_min + 1, std::vector<double>(c_max - c_min + 1, 0));
	for (size_t s = s_min; s <= s_max; s++)
		for (size_t c = c_min; c <= c_max; c++)
			(*sub_mat)[s][c] = mat[s][c];
	return sub_mat;
}

int int_random(int mod, bool positive)
{
	std::random_device rd;
	return (positive ? (rd() % mod) + 1 : (rd() % mod) - (int)(mod / 2));
}

std::vector<std::vector<double>>* transpose(const std::vector<std::vector<double>>& mat)
{
	auto t = new std::vector<std::vector<double>>(mat[0].size(), std::vector<double>(mat.size(), 0));
	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = 0; j < mat[0].size(); j++)
			(*t)[j][i] = mat[i][j];
	return t;
}
