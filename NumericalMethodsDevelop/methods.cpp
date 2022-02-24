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
			{
				(*l_mat)[i][j] = ii_element(*l_mat, mat[i][j], i);
				if ((*l_mat)[i][j] == -1)
				{
					std::cout << "Matrix contain negative eigenvalues";
					return nullptr;
				}
			}
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
	if (a - sum <= 0)
		return -1;
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

std::vector<std::vector<double>>* generate_matrix(
	size_t n,
	size_t mod,
	double min,
	double max,
	Eigen::MatrixXd* diag,
	Eigen::MatrixXd* eigen_vectors
)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(-1.0, 1.0);

	Eigen::MatrixXd matx = Eigen::MatrixXd::NullaryExpr(n, n, [&]() {return dis(gen); });
	Eigen::MatrixXd Q, D;
	if (diag == nullptr)
		D = Eigen::VectorXd::LinSpaced(n, min, max).asDiagonal();
	else
	{
		*diag = Eigen::VectorXd::LinSpaced(n, min, max).asDiagonal();
		D = *diag;
	}	
	while (matx.determinant() == 0)
		matx = Eigen::MatrixXd::NullaryExpr(n, n, [&]() {return dis(gen); });
	matx *= mod;
	Eigen::HouseholderQR<Eigen::MatrixXd> qr(matx);
	if (eigen_vectors == nullptr)
		Q = qr.householderQ();
	else
	{
		*eigen_vectors = qr.householderQ();
		Q = *eigen_vectors;
	}
	matx = Q * D * Q.transpose();
	return eigen_to_matrix(matx);
}

void print_mat(const std::vector<std::vector<double>>& mat, const std::vector<double>& b)
{
	if (b.empty())
		for (size_t i = 0; i < mat.size(); i++)
		{
			for (size_t j = 0; j < mat.size(); j++)
				printf("%.1f, ", mat[i][j]);
			std::cout << std::endl;
		}
	else
		for (size_t i = 0; i < mat.size(); i++)
		{
			for (size_t j = 0; j < mat.size(); j++)
				printf("%.1f, ", mat[i][j]);
			printf("| %.1f\n", b[i]);
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
	auto A = generate_matrix(n, 10, min, max);
	auto x_c = generate_vector(n, 10);
	auto A_e = matrix_to_eigen(*A);
	auto x_e = vector_to_eigen(*x_c);
	auto b_e = A_e * x_e;
	auto b = eigen_to_vector(b_e);
	print_mat(*A, *x_c);

	/*auto pair = SOR(*A, *b, beta, 1e-1);
	auto x = pair->first;
	auto iters = pair->second;
	auto e_x = vector_to_eigen(*x);
	std::cout << "Solve x:" << std::endl << e_x << std::endl;
	auto nx = eigen_to_vector(A_e * e_x - b_e);
	std::cout << "||x - x*|| = " << matrix_inf_norm(v_sub(*x, *x_c)) << std::endl;
	std::cout << "||Ax - b|| = " << matrix_inf_norm(*nx) << std::endl;
	std::cout << "Iterations = " << iters << std::endl;

	delete A;
	delete b;
	delete x;
	delete pair;
	delete x_c;
	delete nx;*/
}

void dependency_3(double min_eps, double max_eps, size_t n, double beta)
{
	FILE* out = fopen("C:/Users/Denis/YandexDisk/Polytech/NumMethods/Lab3/out3.csv", "w");
	fprintf(out, "eps,||x - x*||,||Ax - b||,iters,%.2f\n", beta);

	Eigen::MatrixXd A_e(10, 10);
	Eigen::VectorXd x_e(10);
	/*
	A_e <<
		55.8, 3.1, -6.7, -16.5, -7.7, 1.6, 4.2, -6.9, -4.4, 7.8,
		3.1, 39.5, -14.4, 16.7, -11.5, 2.5, -13.7, -2.6, 15.2, -22.0,
		-6.7, -14.4, 67.3, -6.4, 11.3, 13.2, -4.8, 3.7, 4.5, 4.6,
		-16.5, 16.7, -6.4, 57.8, 5.4, -0.3, -7.9, -5.5, 0.5, -0.1,
		-7.7, -11.5, 11.3, 5.4, 56.6, -2.0, -0.4, -16.3, -11.2, 13.1,
		1.6, 2.5, 13.2, -0.3, -2.0, 33.5, -32.7, 4.1, 2.8, -0.7,
		4.2, -13.7, -4.8, -7.9, -0.4, -32.7, 47.4, -7.4, -0.7, -8.1,
		-6.9, -2.6, 3.7, -5.5, -16.3, 4.1, -7.4, 42.6, 0.9, 3.8,
		-4.4, 15.2, 4.5, 0.5, -11.2, 2.8, -0.7, 0.9, 40.4, -2.0,
		7.8, -22.0, 4.6, -0.1, 13.1, -0.7, -8.1, 3.8, -2.0, 64.1;
	x_e <<
		41.0, 38.0, -36.0, -45.0, 40.0, 41.0, 17.0, -8.0, -7.0, -16.0;
	*/

	A_e <<
		67.4, 2.3, -5.0, -0.1, -0.7, -0.7, -1.8, 2.1, 4.1, -9.5, 
		2.3, 67.1, 7.3, 6.4, 6.6, 2.1, 9.9, 0.9, -3.6, -2.5, 
		-5.0, 7.3, 71.6, 1.8, -0.0, 6.1, 1.9, -0.1, 5.7, -1.2, 
		-0.1, 6.4, 1.8, 78.4, 0.5, 4.0, -3.5, -7.5, 2.6, 1.7, 
		-0.7, 6.6, -0.0, 0.5, 88.3, -8.0, -9.5, 6.2, 2.9, -9.1,
		- 0.7, 2.1, 6.1, 4.0, -8.0, 74.0, -3.8, 4.1, 1.2, -4.7,
		- 1.8, 9.9, 1.9, -3.5, -9.5, -3.8, 82.2, 2.3, -7.6, -5.6,
		2.1, 0.9, -0.1, -7.5, 6.2, 4.1, 2.3, 65.8, -1.0, -0.3,
		4.1, -3.6, 5.7, 2.6, 2.9, 1.2, -7.6, -1.0, 72.9, -1.8,
		- 9.5, -2.5, -1.2, 1.7, -9.1, -4.7, -5.6, -0.3, -1.8, 82.3;
	x_e <<
		-1.0, -5.0, -2.0, 3.0, 1, 2, -1, -2, 2, 2;

	Eigen::VectorXd b_e = A_e * x_e;

	for (double eps = min_eps; max_eps + min_eps > eps; eps *= 10)
	{
		auto pair = SOR(A_e, b_e, beta, eps);
		fprintf(out, "%.16f,%.16f,%.16f,%d\n", eps, (pair.first - x_e).norm(), (A_e * pair.first - b_e).norm(), pair.second);
	}
	fclose(out);
}

std::pair<Eigen::VectorXd, size_t> SOR(
	const Eigen::MatrixXd& mat,
	const Eigen::VectorXd& b,
	double beta,
	double eps
)
{
	Eigen::VectorXd v0 = Eigen::VectorXd::Zero(mat.cols());
	Eigen::VectorXd v1 = Eigen::VectorXd::Ones(mat.cols());
	size_t iters = 0;
	//(mat * v1 - b)
	while ((v1 - v0).norm() >= eps)
	{
		iters += 1;
		v0 = v1;
		for (size_t i = 0; i < mat.cols(); i++)
		{
			v1(i) = 0;
			for (size_t j = 0; j < mat.cols(); j++)
			{
				double z = 0;
				if (i == j)
					z = ((1 - beta) * v0(j));
				if (i < j)
					z = ((-beta) * v0(j) * (mat(i, j) / mat(i, i)));
				if (i > j)
					z = ((-beta) * v1(j) * (mat(i, j) / mat(i, i)));
				v1(i) += z;
			}
			v1(i) += ((beta * b(i)) / mat(i, i));
		}
	}
	return std::pair<Eigen::VectorXd, size_t>(v1, iters);
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

void lab_4()
{
	size_t n;
	double min, max, eps_j = 1e-1, eps_inv = 1e-16;
	std::cout << "Input n-dim:" << std::endl;
	std::cin >> n;
	std::cout << "Input E_1 and E_n for generate [E_1, E_n] linspace eigen numbers:" << std::endl;
	std::cin >> min >> max;
	if (min > max)
	{
		max = min - max;
		min = min - max;
		max = min + max;
	}
	std::cout << "Epsilon Jacobi = " << eps_j << std::endl << "Epsilon Inverse iterations = " << eps_inv << std::endl;
	auto c_eigen = new Eigen::MatrixXd();
	auto c_eigen_mat = new Eigen::MatrixXd();
	auto A = generate_matrix(n, 100, min, max, c_eigen, c_eigen_mat);
	std::cout << "Generated matrix A:\n";
	print_mat(*A);
	auto A_e = matrix_to_eigen(*A);
	auto res_pair = Jacobi(*A, eps_j);
	auto eigen_e = matrix_to_eigen(*res_pair->first);
	auto rotation_vectors = res_pair->second;

	std::cout << "Strict eigen vectors:\n" << *c_eigen_mat << "\n\n";

	Eigen::MatrixXd T = Eigen::MatrixXd::Identity(n, n);
	for (size_t i = 0; i < (*rotation_vectors).size(); i++)
		T *= matrix_to_eigen(*(*rotation_vectors)[i]);
	std::cout << "Rotation matrix G:\n" << T << "\n\n";

	auto d1 = eigen_e.diagonal();
	auto rotation_matrix_columns_map = std::map<double, size_t>();
	for (size_t i = 0; i < d1.size(); i++)
		rotation_matrix_columns_map[d1(i)] = i;
	auto d2 = (*c_eigen).diagonal();

	std::sort(d1.begin(), d1.end(), std::less<double>());
	std::sort(d2.begin(), d2.end(), std::less<double>());

	std::cout << "Eigen value : eigen vector\n";
	for (size_t i = 0; i < d1.size(); i++)
	{
		std::cout << d1(i) << " : (" << T.col(rotation_matrix_columns_map[d1(i)]).transpose() << ")'\n";
		std::cout << "||A*x - " << d1(i) << "*x|| = " << (A_e * T.col(rotation_matrix_columns_map[d1(i)]) - d1(i) * T.col(rotation_matrix_columns_map[d1(i)])).norm() << "\n";
	}

	std::cout << "Iterations = " << rotation_vectors->size() << "\n\n";
	std::cout << "Strict eigen values: (" << d2.transpose() << ")'\n";
	std::cout << "Approximate eigen values: (" << d1.transpose() << ")'\n";
	std::cout << "||eigen - eigen*|| = " << (d1 - d2).norm() << "\n\n";

	std::cout << "Inverse iteration corrected eigen values and eigen vectors:\n";
	size_t iters = 0;
	/*for (size_t i = 0; i < d1.size(); i++)
	{
		auto t = inv_iters(*A, d1(i), eps_inv);
		double exact_eigen_value = std::get<0>(*t);
		auto exact_eigen_vector = std::get<1>(*t);
		iters += std::get<2>(*t);
		std::cout << exact_eigen_value << " : (" << exact_eigen_vector.transpose() << ")'\n";
		double nrm = get_signed_max_abs((*c_eigen_mat).col(i));
		std::cout << "||x - x*|| = " << (((*c_eigen_mat).col(i) / nrm).cwiseAbs() - exact_eigen_vector.cwiseAbs()).norm() << "\n";
		std::cout << "||A*x - " << d2(i) << "*x|| = " << (A_e * exact_eigen_vector - exact_eigen_value * exact_eigen_vector).norm() << "\n";
	}
	std::cout << "\nEigen value and corresponding ||x - x*||:\n";
	for (size_t i = 0; i < d1.size(); i++)
		std::cout << d2(i) << " : " << ((*c_eigen_mat).col(i).cwiseAbs() - T.col(rotation_matrix_columns_map[d1(i)]).cwiseAbs()).norm() << "\n";*/
}

void dependency_4(double min_eps, double max_eps, size_t n)
{
	FILE* out = fopen("C:/Users/Denis/YandexDisk/Polytech/NumMethods/Lab4/out4.csv", "w");
	fprintf(
		out,
		"eps,||x - x*||,||eigen - eigen*||,||A*x - eigen*x||,iters,||x - x*||e,||eigen - eigen*||e,||A*x - eigen*x||e,inv_iters\n"
	);
	size_t itc = 0;
	auto c_eigen = new Eigen::MatrixXd();
	auto c_eigen_mat = new Eigen::MatrixXd();
	auto A = generate_matrix(n, 20, 1, 10, c_eigen, c_eigen_mat);
	auto A_e = matrix_to_eigen(*A);
	for (double eps = min_eps; max_eps + min_eps > eps; eps *= 10)
	{
		if (itc % 1000 == 0)
			printf("%.2f\n", (itc * 100 / (max_eps / min_eps)));
		itc++;
		auto res_pair = Jacobi(*A, eps);
		auto eigen_e = matrix_to_eigen(*res_pair->first);
		auto rotation_vectors = res_pair->second;

		Eigen::MatrixXd T = Eigen::MatrixXd::Identity(n, n);
		for (size_t i = 0; i < (*rotation_vectors).size(); i++)
			T *= matrix_to_eigen(*(*rotation_vectors)[i]);

		auto d1 = eigen_e.diagonal();
		auto rotation_matrix_columns_map = std::map<double, size_t>();
		for (size_t i = 0; i < d1.size(); i++)
			rotation_matrix_columns_map[d1(i)] = i;
		auto d2 = (*c_eigen).diagonal();

		std::sort(d1.begin(), d1.end(), std::less<double>());
		std::sort(d2.begin(), d2.end(), std::less<double>());

		double euality_error_sum = 0, vectors_error_sum = 0, exact_vectors_error_sum = 0, exact_euality_error_sum = 0;
		for (size_t i = 0; i < d1.size(); i++)
		{
			euality_error_sum += (A_e * T.col(rotation_matrix_columns_map[d1(i)]) - d1(i) * T.col(rotation_matrix_columns_map[d1(i)])).norm();
			vectors_error_sum += ((*c_eigen_mat).col(i).cwiseAbs() - T.col(rotation_matrix_columns_map[d1(i)]).cwiseAbs()).norm();
		}

		Eigen::VectorXd d3(d1.size());
		size_t iters = 0;
		for (size_t i = 0; i < d1.size(); i++)
		{
			auto t = inv_iters(*A, T.col(rotation_matrix_columns_map[d1(i)]), d1(i), eps);
			double exact_eigen_value = std::get<0>(*t);
			auto exact_eigen_vector = std::get<1>(*t);
			iters += std::get<2>(*t);
			d3(i) = exact_eigen_value;
			exact_euality_error_sum += (A_e * exact_eigen_vector - exact_eigen_value * exact_eigen_vector).norm();
			double nrm = get_signed_max_abs((*c_eigen_mat).col(i));
			exact_vectors_error_sum += (((*c_eigen_mat).col(i) / nrm) - exact_eigen_vector).norm();
			delete t;
		}

		fprintf(
			out,
			"%.16f,%.16f,%.16f,%.16f,%d,%.16f,%.16f,%.16f,%d\n",
			eps, vectors_error_sum, (d1 - d2).norm(), euality_error_sum, rotation_vectors->size(),
			exact_vectors_error_sum, (d3 - d2).norm(), exact_euality_error_sum, iters / n
		);

		delete res_pair->first;
		for (size_t i = 0; i < (*(*res_pair).second).size(); i++)
			delete (*(*res_pair).second)[i];
		delete (*res_pair).second;
		delete res_pair;
	}
	delete A;
	delete c_eigen;

	fclose(out);
}

std::vector<std::vector<double>>* matrix_mul(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b)
{
	auto m = new std::vector<std::vector<double>>(a.size(), std::vector<double>(b[0].size(), 0));
	for (size_t i = 0; i < a.size(); i++)
		for (size_t j = 0; j < b[0].size(); j++)
			for (size_t k = 0; k < a[0].size(); k++)
				(*m)[i][j] += a[i][k] * b[k][j];
	return m;
}

std::vector<std::vector<double>>* Jacobi_rotation(
	std::vector<std::vector<double>>& mat,
	double sin_phi,
	double cos_phi,
	size_t opt_i,
	size_t opt_j
)
{
	auto rot_matrix = new std::vector<std::vector<double>>(mat.size(), std::vector<double>(mat.size(), 0));
	for (size_t i = 0; i < mat.size(); i++)
		(*rot_matrix)[i][i] = 1;
	(*rot_matrix)[opt_i][opt_i] = cos_phi;
	(*rot_matrix)[opt_j][opt_j] = cos_phi;
	(*rot_matrix)[opt_i][opt_j] = -sin_phi;
	(*rot_matrix)[opt_j][opt_i] = sin_phi;
	for (size_t k = 0; k < mat.size(); k++)
	{
		 double x = cos_phi * mat[opt_i][k] + sin_phi * mat[opt_j][k];
		 double y = -sin_phi * mat[opt_i][k] + cos_phi * mat[opt_j][k];
		 mat[opt_i][k] = x;
		 mat[opt_j][k] = y;
	}
	for (size_t k = 0; k < mat.size(); k++)
	{
		double x = cos_phi * mat[k][opt_i] + sin_phi * mat[k][opt_j];
		double y = -sin_phi * mat[k][opt_i] + cos_phi * mat[k][opt_j];
		mat[k][opt_i] = x;
		mat[k][opt_j] = y;
	}
	return rot_matrix;
}

std::pair<std::vector<std::vector<double>>*, std::vector<std::vector<std::vector<double>>*>*>* Jacobi(
	const std::vector<std::vector<double>>& mat,
	double eps
)
{
	auto A = new std::vector<std::vector<double>>(mat.size(), std::vector<double>(mat.size(), 0));
	auto rotation_matrixes = new std::vector<std::vector<std::vector<double>>*>();
	size_t iters = 0;

	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = 0; j < mat[0].size(); j++)
			(*A)[i][j] = mat[i][j];

	while (ndss(*A, true) >= eps)
	{
		auto opt = select_opt(*A);
		size_t opt_i = opt.first, opt_j = opt.second;
		double sin_phi = 1, cos_phi = 0;
		if ((*A)[opt_i][opt_i] != (*A)[opt_j][opt_j])
		{
			double tan_2phi = 2 * (*A)[opt_i][opt_j] / ((*A)[opt_i][opt_i] - (*A)[opt_j][opt_j]);
			sin_phi = ((tan_2phi > 0) - (tan_2phi < 0)) * sqrt(0.5 * (1 - (1 / sqrt(1 + tan_2phi * tan_2phi))));
			cos_phi = sqrt(0.5 * (1 + (1 / sqrt(1 + tan_2phi * tan_2phi))));
		}
		rotation_matrixes->push_back(Jacobi_rotation(*A, sin_phi, cos_phi, opt_i, opt_j));
	}
	auto res = new std::pair<std::vector<std::vector<double>>*, std::vector<std::vector<std::vector<double>>*>*>(A, rotation_matrixes);
	return res;
}

std::tuple<double, Eigen::VectorXd, size_t>* inv_iters(const std::vector<std::vector<double>>& mat, Eigen::VectorXd v_input, double lambda, double eps)
{
	Eigen::MatrixXd A, E;
	E = lambda * Eigen::VectorXd::Ones(mat.size()).asDiagonal();
	A = matrix_to_eigen(mat);
	Eigen::FullPivLU<Eigen::MatrixXd> lu(A - E);
	Eigen::VectorXd y0 = v_input;
	Eigen::VectorXd y1 = Eigen::VectorXd::Zero(mat.size());

	double u = get_signed_max_abs(y0);
	size_t iters = 0;
	for (size_t stop_iter = 0; stop_iter < 1000; stop_iter++)
	{
		iters++;
		y1 = lu.solve(y0 / u);
		u = get_signed_max_abs(y1);
		if ((A * y1 - (lambda + (1 / u)) * y1).norm() < 1e-6)
			break;
		y0 = y1;
	}
	return new std::tuple<double, Eigen::VectorXd, size_t>(lambda + (1 / u), y1 / u, iters);
}

double get_signed_max_abs(const Eigen::VectorXd& v)
{
	double m = std::numeric_limits<double>::min();
	for (size_t i = 0; i < v.size(); i++)
		if (std::abs(v(i)) > std::abs(m))
			m = v(i);
	return m;
}

double ndss(const std::vector<std::vector<double>>& mat, bool t)
{
	double sum = 0, m = std::numeric_limits<double>::min();
	if (t)
	{
		for (size_t i = 0; i < mat.size(); i++)
			for (size_t j = i + 1; j < mat[0].size(); j++)
				if (mat[i][j] * mat[i][j] > m)
					m = mat[i][j] * mat[i][j];
		return m;
	}
	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = i + 1; j < mat[0].size(); j++)
			sum += mat[i][j] * mat[i][j];
	return sum;
}

std::pair<size_t, size_t> select_opt(const std::vector<std::vector<double>>& mat)
{
	size_t opt_i = 0, opt_j = 0;
	double r_i, max_r = 0;
	for (size_t i = 0; i < mat.size(); i++)
	{
		r_i = 0;
		for (size_t j = 0; j < mat[0].size(); j++)
			if (i != j)
				r_i += mat[i][j] * mat[i][j];
		if (max_r < r_i)
		{
			max_r = r_i;
			opt_i = i;
		}
	}
	max_r = 0;
	for (size_t j = 0; j < mat[0].size(); j++)
		if ((j != opt_i) && (max_r < std::abs(mat[opt_i][j])))
		{
			opt_j = j;
			max_r = std::abs(mat[opt_i][j]);
		}
	return std::pair<size_t, size_t>(opt_i, opt_j);
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
