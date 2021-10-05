#include <iostream>
#include "methods.h"

int main()
{
	double l, r;
	std::cin >> l >> r;
	std::cout << simple_iteration(l, r, r, 1e-5);
	return 0;
}