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

std::vector<std::vector<double>>* generate_matrix(size_t n, size_t mod, double min, double max, Eigen::MatrixXd* diag)
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

	auto pair = SOR(*A, *b, beta, 1e-3);
	auto x = pair->first;
	auto iters = pair->second;
	auto e_x = vector_to_eigen(*x);
	std::cout << "Solve x:" << std::endl << e_x << std::endl;
	auto dx = eigen_to_vector(x_c);
	auto nx = eigen_to_vector(A_e * e_x - b_e);
	std::cout << "||x - x*|| = " << matrix_inf_norm(v_sub(*x, *dx)) << std::endl;
	std::cout << "||Ax - b|| = " << matrix_inf_norm(*nx) << std::endl;
	std::cout << "Iterations = " << iters << std::endl;

	delete A;
	delete b;
	delete x;
	delete pair;
	delete dx;
	delete nx;
}

void dependency_3(double min_eps, double max_eps, size_t n, double beta)
{
	FILE* out = fopen("out3.csv", "w");
	fprintf(out, "eps,||x - x*||,||Ax - b||,iters\n");

	for (double eps = min_eps; eps < max_eps; eps += min_eps)
	{
		auto A = generate_matrix(n, 100, 1, 100);
		auto b = generate_vector(n, 10);
		auto A_e = matrix_to_eigen(*A);
		auto b_e = vector_to_eigen(*b);
		auto x_c = A_e.inverse() * b_e;
		auto pair = SOR(*A, *b, beta, eps);
		auto x = pair->first;
		auto iters = pair->second;
		auto e_x = vector_to_eigen(*x);
		auto dx = eigen_to_vector(x_c);
		auto discrepancy = eigen_to_vector(A_e * e_x - b_e);

		fprintf(out, "%.16f,%.16f,%.16f,%d\n", eps, matrix_inf_norm(v_sub(*x, *dx)), matrix_inf_norm(*discrepancy), iters);

		delete A;
		delete b;
		delete x;
		delete pair;
		delete dx;
		delete discrepancy;
	}

	fclose(out);
}

std::pair<std::vector<double>*, size_t>* SOR(
	const std::vector<std::vector<double>>& mat,
	const std::vector<double>& b,
	double beta,
	double eps
)
{
	if ((beta <= 0) || (2 <= beta))
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
	size_t iters = 0;
	while (matrix_inf_norm(v_sub(*x0, *x1)) >= eps)
	{
		iters += 1;
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
	}
	auto res = new std::pair<std::vector<double>*, size_t>(x1, iters);
	delete x0;
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

void lab_4()
{
	size_t n;
	double min, max, eps_j = 1e-1, eps_inv = 1e-16;
	std::cout << "Input n-dim:" << std::endl;
	std::cin >> n;
	std::cout << "Input E_1 and E_n for generate [E_1, E_n] linspace eigen numbers:" << std::endl;
	std::cin >> min >> max;
	std::cout << "Epsilon Jacobi = " << eps_j << std::endl << "Epsilon Inverse iterations = " << eps_inv << std::endl;
	auto c_eigen = new Eigen::MatrixXd();
	auto A = generate_matrix(n, 100, min, max, c_eigen);
	std::cout << "Generated matrix A:\n";
	print_mat(*A);
	auto A_e = matrix_to_eigen(*A);
	auto res_pair = Jacobi(*A, eps_j);
	auto eigen_e = matrix_to_eigen(*res_pair->first);
	auto rotation_vectors = res_pair->second;

	auto T = matrix_to_eigen(*(*rotation_vectors)[0]);
	for (size_t i = 1; i < (*rotation_vectors).size(); i++)
		T *= matrix_to_eigen(*(*rotation_vectors)[i]);
	std::cout << "Rotation matrix G:\n" << T << "\n\n";

	auto d1 = eigen_e.diagonal();
	auto rotation_matrix_columns_map = std::map<double, size_t>();
	for (size_t i = 0; i < d1.size(); i++)
		rotation_matrix_columns_map[d1(i)] = i;
	auto d2 = (*c_eigen).diagonal();

	for (size_t i = 0; i < d1.size(); i++)
	{
		std::cout << "Eigen value : eigen vector\n";
		std::cout << d1(i) << " : (" << T.col(i).transpose() << ")'\n";
		std::cout << "|| A*x - " << d1(i) << "*x || = " << (A_e * T.col(i) - d1(i) * T.col(i)).norm() << "\n";
	}

	std::sort(d1.begin(), d1.end(), std::less<double>());
	std::sort(d2.begin(), d2.end(), std::less<double>());

	std::cout << "Iterations = " << rotation_vectors->size() << "\n\n";
	std::cout << "Exact eigen values: (" << d2.transpose() << ")'\n";
	std::cout << "Not exact eigen values: (" << d1.transpose() << ")'\n";
	std::cout << "||eigen - eigen*|| = " << (d1 - d2).norm() << "\n\n";

	std::cout << "Inverse iteration corrected eigen values and eigen vectors:\n";
	std::vector<Eigen::VectorXd> e_v;
	for (size_t i = 0; i < d1.size(); i++)
	{
		auto p = inv_iters(*A, d1(i), eps_inv);
		e_v.push_back(p->second);
		std::cout << p->first << " : (" << p->second.transpose() << ")'\n";
		std::cout << "|| A*x - " << d2(i) << "*x || = " << (A_e * p->second - d2(i) * p->second).norm() << "\n";
	}
	std::cout << "\nEigen value and corresponding || x - *x ||:\n";
	for (size_t i = 0; i < d1.size(); i++)
	{
		std::cout << d2(i) << " : " << (e_v[i] - T.col(rotation_matrix_columns_map[d1(i)])).norm() << "\n" << e_v[i].transpose() << "\n" << T.col(rotation_matrix_columns_map[d1(i)]).transpose() << "\n";
	}
}

void dependency_4(double min_eps, double max_eps, size_t n)
{

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

	while (ndss(*A) >= eps)
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

std::pair<double, Eigen::VectorXd>* inv_iters(const std::vector<std::vector<double>>& mat, double lambda, double eps)
{
	Eigen::MatrixXd A, B, E;
	E = lambda * Eigen::VectorXd::Ones(mat.size()).asDiagonal();
	A = matrix_to_eigen(mat);
	B = (A - E).inverse();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0, 1);
	Eigen::VectorXd y0 = Eigen::VectorXd::NullaryExpr(mat.size(), [&]() {return dis(gen); });
	Eigen::VectorXd y1 = Eigen::VectorXd::Zero(mat.size());
	double u = get_signed_max_abs(y0);
	while (true)
	{
		y1 = B * (y0 / u);
		u = get_signed_max_abs(y1);
		if ((y1 - y0).norm() < eps)
			break;
		y0 = y1;
	}
	return new std::pair<double, Eigen::VectorXd>(lambda + (1 / u), y1 / u);
}

double get_signed_max_abs(const Eigen::VectorXd& v)
{
	double m = std::numeric_limits<double>::min();
	for (size_t i = 0; i < v.size(); i++)
		if (std::abs(v(i)) > std::abs(m))
			m = v(i);
	return m;
}

double ndss(const std::vector<std::vector<double>>& mat)
{
	double sum = 0;
	for (size_t i = 0; i < mat.size(); i++)
		for (size_t j = 0; j < mat[0].size(); j++)
			if (i != j)
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
	//if(opt_i > opt_j)
	//	return std::pair<size_t, size_t>(opt_j, opt_i);
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
