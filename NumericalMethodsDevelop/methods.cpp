#include "methods.h"

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

void line_sum(std::vector<std::vector<double>>& mat, size_t src_index, size_t dst_index, double k, bool transpose)
{
	if (transpose)
		for (size_t i = 0; i < mat.size(); i++)
			mat[i][dst_index] += k * mat[i][src_index];
	else
		for (size_t i = 0; i < mat.size(); i++)
			mat[dst_index][i] += k * mat[src_index][i];
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
*/

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
					sub_mat[i - 1].push_back(mat[i][k]);
		d += pow(-1, j) * mat[0][j] * det(sub_mat);
	}
	return d;
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
		std::vector<std::vector<double>> sub_mat(i + 1, std::vector<double>(i + 1, 0));
		for (size_t g = 0; g < i + 1; g++)
			for (size_t k = 0; k < i + 1; k++)
				sub_mat[g][k] = mat[g][k];
		d = det(sub_mat);
		std::cout << "Minor " << i + 1 << "x" << i + 1 << " = " << d << std::endl;
		if (d <= 0)
			return false;
	}
	std::cout << "Succsess!" << std::endl;
	return true;
}


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

void dependency(double min, double max, size_t n, double stride)
{
	FILE* out = fopen("out.csv", "w");
	fprintf(out, "cond,||x - x*||,||Ax - b||\n");

	for (double c = min; c < max; c += stride)
	{
		auto A = generate_matrix(n, 100, min, c);
		auto b = generate_vector(n, 10);
		auto S = cholesky(*A);
		auto S_t = transpose(*S);
		auto x = solve_system(*S, *b);
		auto e_x = vector_to_eigen(*x);
		auto A_e = matrix_to_eigen(*A);
		auto b_e = vector_to_eigen(*b);
		auto x_c = A_e.inverse() * b_e;
		double cond = A_e.norm() * A_e.norm();

		fprintf(out, "%.16f,%.16f,%.16f\n", cond, (x_c - e_x).norm(), (A_e * e_x - b_e).norm());

		delete A;
		delete S;
		delete S_t;
		delete b;
		delete x;
	}

	fclose(out);
}

void lab_2()
{
	size_t n;
	double min, max;
	std::cout << "Input n-dim:" << std::endl;
	std::cin >> n;
	std::cout << "Input real D_min and D_max for generate condition number:" << std::endl;
	std::cin >> min >> max;
	std::cout << "Generate linear system Ax = b: " << n << "x" << n << std::endl;
	auto A = generate_matrix(n, 100, min, max);
	if (!sylvester_criterion(*A))
	{
		std::cout << "Matrix not positive defined" << std::endl;
		return;
	}
	auto b = generate_vector(n, 10);
	auto A_e = matrix_to_eigen(*A);
	auto b_e = vector_to_eigen(*b);
	auto x_c = A_e.inverse() * b_e;
	double A_cond = A_e.norm() * A_e.inverse().norm();
	std::cout << "Condition number = " << A_cond << std::endl;
	print_mat(*A, *b);

	auto S = cholesky(*A);
	auto S_t = transpose(*S);
	std::cout << "L" << std::endl;
	print_mat(*S);
	std::cout << "L'" << std::endl;
	print_mat(*S_t);

	auto x = solve_system(*S, *b);
	auto e_x = vector_to_eigen(*x);
	std::cout << "Solve x:" << std::endl << e_x << std::endl;
	std::cout << "||x - x*|| = " << (x_c - e_x).norm() << std::endl;
	std::cout << "||Ax - b|| = " << (A_e * e_x - b_e).norm() << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(-1.0, 1.0);
	auto b_e_delta = Eigen::VectorXd::NullaryExpr(n, [&]() {return dis(gen) / 100; });
	auto x_e_delta = A_e.inverse() * b_e_delta;
	std::cout << "||delta_x|| = " << x_e_delta.norm() << std::endl;
	std::cout << "||x|| = " << e_x.norm() << std::endl;
	std::cout << "cond = " << A_cond << std::endl;
	std::cout << "||delta_b|| = " << b_e_delta.norm() << std::endl;
	std::cout << "||b|| = " << b_e.norm() << std::endl;
	double _left = x_e_delta.norm() / e_x.norm();
	double _right = b_e_delta.norm() / b_e.norm();
	std::cout << _left << " <= " << A_cond << "*" << _right << " - " << (_left <= A_cond * _right ? "true" : "false") << std::endl;

	delete A;
	delete S;
	delete S_t;
	delete b;
	delete x;
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
	return eigen_to_matrix(matx);
}

void print_mat(const std::vector<std::vector<double>>& mat, const std::vector<double>& b)
{
	if (b.empty())
		for (size_t i = 0; i < mat.size(); i++)
		{
			for (size_t j = 0; j < mat.size(); j++)
				printf("%.8f ", mat[i][j]);
			std::cout << std::endl;
		}
	else
		for (size_t i = 0; i < mat.size(); i++)
		{
			for (size_t j = 0; j < mat.size(); j++)
				printf("%.8f ", mat[i][j]);
			printf("| %.8f\n", b[i]);
		}
	std::cout << std::endl;
}

std::vector<double>* generate_vector(size_t len, int mod, bool positive)
{
	auto v = new std::vector<double>();
	for (size_t s = 0; s < len; s++)
		(*v).push_back(int_random(mod, positive));
	return v;
}

std::vector<double>* solve_system(const std::vector<std::vector<double>>& mat, const std::vector<double>& v)
{
	auto S = matrix_to_eigen(mat);
	auto b = vector_to_eigen(v);
	Eigen::VectorXd y = S.triangularView<Eigen::Lower>().solve(b);
	Eigen::VectorXd x = S.transpose().triangularView<Eigen::Upper>().solve(y);
	return eigen_to_vector(x);
}

Eigen::VectorXd vector_to_eigen(const std::vector<double>& v)
{
	Eigen::VectorXd e_v(v.size());
	for (size_t i = 0; i < v.size(); i++)
		e_v(i, 0) = v[i];
	return e_v;
}

std::vector<double>* eigen_to_vector(Eigen::VectorXd v)
{
	auto v_v = new std::vector<double>(v.rows(), 0);
	for (size_t i = 0; i < v.rows(); i++)
		(*v_v)[i] = v(i, 0);
	return v_v;
}

Eigen::MatrixXd matrix_to_eigen(const std::vector<std::vector<double>>& mat)
{
	Eigen::MatrixXd e_mat(mat.size(), mat[0].size());
	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = 0; j < mat[0].size(); j++)
			e_mat(i, j) = mat[i][j];
	return e_mat;
}

std::vector<std::vector<double>>* eigen_to_matrix(Eigen::MatrixXd mat)
{
	auto v_mat = new std::vector<std::vector<double>>(mat.rows(), std::vector<double>(mat.cols(), 0));
	for (size_t i = 0; i < mat.rows(); i++)
		for (size_t j = 0; j < mat.cols(); j++)
			(*v_mat)[i][j] = mat(i, j);
	return v_mat;
}

void lab_3()
{
	size_t n;
	double min, max, beta;
	std::cout << "Input n-dim:" << std::endl;
	std::cin >> n;
	std::cout << "Input real D_11 and D_nn for generate determinant:" << std::endl;
	std::cin >> min >> max;
	std::cout << "Input real beta between (0, 2):" << std::endl;
	std::cin >> beta;
	if (max - min <= 0)
	{
		std::cout << "Determinant must be > 0" << std::endl;
		return;
	}
	std::cout << "Generate linear system Ax = b: " << n << "x" << n << std::endl;
	auto A = generate_matrix(n, 100, min, max);
	if (!sylvester_criterion(*A))
	{
		std::cout << "Matrix not positive defined" << std::endl;
		return;
	}
	auto b = generate_vector(n, 10);
	auto A_e = matrix_to_eigen(*A);
	auto b_e = vector_to_eigen(*b);
	auto x_c = A_e.inverse() * b_e;
	std::cout << "det(A) = " << A_e.determinant() << std::endl;
	print_mat(*A, *b);

	std::cout << "SOR" << std::endl << std::endl;
	auto pair = SOR(*A, *b, beta, 1e-3);
	auto x = pair->first;
	auto iters = pair->second;
	auto e_x = vector_to_eigen(*x);
	std::cout << "Solve x:" << std::endl << e_x << std::endl;
	auto dx = eigen_to_vector(x_c);
	auto nx = eigen_to_vector(A_e * e_x - b_e);
	std::cout << "||x - x*||inf = " << matrix_inf_norm(v_sub(*x, *dx)) << std::endl;
	std::cout << "||Ax - b||inf = " << matrix_inf_norm(*nx) << std::endl;
	std::cout << "||x - x*|| = " << (x_c - e_x).norm() << std::endl;
	std::cout << "||Ax - b|| = " << (A_e * e_x - b_e).norm() << std::endl;
	std::cout << "Iterations = " << iters << std::endl;

	delete A;
	delete b;
	delete x;
	delete pair;
	delete dx;
	delete nx;
}

void dependency_3(double min, double max, size_t n, double stride)
{
}

std::pair<std::vector<double>*, double>* SOR(
	const std::vector<std::vector<double>>& mat,
	const std::vector<double>& b,
	double beta,
	double eps
)
{
	if ((0 < beta) && (beta < 1))
		std::cout << beta << " - lower relaxation" << std::endl;
	else if ((1 < beta) && (beta < 2))
		std::cout << beta << " - upper relaxation" << std::endl;
	else if (beta == 1)
		std::cout << beta << " - Seidel method" << std::endl;
	else
	{
		std::cout << beta << " - must be in (0, 2)" << std::endl;
		return nullptr;
	}
	if (eps <= 0)
	{
		std::cout << eps << " - must be > 0" << std::endl;
		return nullptr;
	}

	auto x0 = new std::vector<double>(mat.size(), 0);
	auto x1 = new std::vector<double>(mat.size(), 1);
	double iters = 0;
	double delta = matrix_inf_norm(v_sub(*x0, *x1));
	while (delta >= eps)
	{
		iters += 1;
		std::cout << delta << " = delta " << iters << " = iters" << std::endl;
		for (size_t i = 0; i < (*x0).size(); i++)
			(*x0)[i] = (*x1)[i];
		for (size_t i = 0; i < mat.size(); i++)
		{
			(*x1)[i] = 0;
			for (size_t j = 0; j < mat.size(); j++)
			{
				if (i == j)
					(*x1)[i] += (1 - beta) * (*x0)[j];
				if (i < j)
					(*x1)[i] -= beta * (*x0)[j] * mat[i][j] / mat[i][i];
				if (i > j)
					(*x1)[i] -= beta * (*x1)[j] * mat[i][j] / mat[i][i];
			}
			(*x1)[i] += beta * b[i] / mat[i][i];
		}
		delta = matrix_inf_norm(v_sub(*x0, *x1));
	}
	auto res = new std::pair<std::vector<double>*, double>(x1, iters);
	return res;
}

std::vector<double> v_sub(const std::vector<double>& x0, const std::vector<double>& x1)
{
	auto new_v = std::vector<double>();
	for (size_t i = 0; i < x0.size(); i++)
		new_v.push_back(x0[i] - x1[i]);
	return new_v;
}

double matrix_inf_norm(const std::vector<std::vector<double>>& mat)
{
	auto v = std::vector<double>();
	for (size_t i = 0; i < mat.size(); i++)
	{
		v.push_back(0);
		for (size_t j = 0; j < mat[i].size(); j++)
			v[i] += std::abs(mat[i][j]);
	}
	return max(v);
}

double matrix_inf_norm(const std::vector<double>& v)
{
	auto new_v = std::vector<double>();
	for (size_t i = 0; i < v.size(); i++)
		new_v.push_back(std::abs(v[i]));
	return max(new_v);
}

double max(const std::vector<double>& v)
{
	double m = std::numeric_limits<double>::min();
	for (auto e : v)
		if (m < e)
			m = e;
	return m;
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
