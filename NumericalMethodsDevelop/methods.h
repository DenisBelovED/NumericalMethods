#pragma once
#undef NDEBUG
#include <cassert>
#include <math.h>
#include <vector>
#include <random>

/*
double half_division(double left_edge, double right_edge, double epsilon=1e-3);

double half_division_test(double left_edge, double right_edge, double epsilon = 1e-3);

double simple_iteration(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c=0.011930646237607);

double simple_iteration_t(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c = 0.0018327550019109917);

double simple_iteration_test(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c = 0.0018327550019109917);
*/

std::vector<std::vector<double>>* cholesky(std::vector<std::vector<double>>& mat);

double ii_element(std::vector<std::vector<double>>& t_mat, double a, size_t line);

double ij_element(std::vector<std::vector<double>>& t_mat, double a, size_t line, size_t col);

void lab_2();

double det(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>>* generate_matrix(size_t n, size_t mod = 10, double min = 0.01, double max = 1);

void line_sum(std::vector<std::vector<double>>& mat, size_t src_index, size_t dst_index, double k, bool transpose=false);

void print_mat(const std::vector<std::vector<double>>& mat);

bool sylvester_criterion(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>>* get_submat(const std::vector<std::vector<double>>& mat, size_t s_min, size_t c_min, size_t s_max, size_t c_max);

int int_random(int mod, bool positive = false);

std::vector<std::vector<double>>* transpose(const std::vector<std::vector<double>>& mat);