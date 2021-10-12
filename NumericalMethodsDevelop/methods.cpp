#include "methods.h"
#include <iostream>
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