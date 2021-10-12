#pragma once
#undef NDEBUG
#include <cassert>
#include <math.h>

double half_division(double left_edge, double right_edge, double epsilon=1e-3);

double simple_iteration(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c=0.011930646237607);

double simple_iteration_t(double left_edge, double right_edge, double x_0, double epsilon = 1e-3, double c = 0.0018327550019109917);

