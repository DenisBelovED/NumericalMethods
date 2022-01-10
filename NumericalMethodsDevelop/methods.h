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

double det(const std::vector<std::vector<double>>& mat);

void line_sum(std::vector<std::vector<double>>& mat, size_t src_index, size_t dst_index, double k, bool transpose=false);

bool sylvester_criterion(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>>* get_submat(const std::vector<std::vector<double>>& mat, size_t s_min, size_t c_min, size_t s_max, size_t c_max);
*/

std::vector<std::vector<double>>* cholesky(std::vector<std::vector<double>>& mat);

double ii_element(std::vector<std::vector<double>>& t_mat, double a, size_t line);

double ij_element(std::vector<std::vector<double>>& t_mat, double a, size_t line, size_t col);

void lab_2();

void dependency(double min, double max, size_t n, double stride = 1);

std::vector<std::vector<double>>* generate_matrix(size_t n, size_t mod = 10, double min = 0.01, double max = 1);

void print_mat(const std::vector<std::vector<double>>& mat, const std::vector<double>& b = std::vector<double>());

int int_random(int mod, bool positive = false);

std::vector<std::vector<double>>* transpose(const std::vector<std::vector<double>>& mat);

std::vector<double>* generate_vector(size_t len, int mod, bool positive=false);

std::vector<double>* solve_system(const std::vector<std::vector<double>>& mat, const std::vector<double>& v);

Eigen::VectorXd vector_to_eigen(const std::vector<double>& v);

std::vector<double>* eigen_to_vector(Eigen::VectorXd v);

Eigen::MatrixXd matrix_to_eigen(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>>* eigen_to_matrix(Eigen::MatrixXd mat);
