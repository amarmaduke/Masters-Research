#include <iostream>

#include "jacobi.h"
#include "../json/json.h"
#include "main.h"

/*
int jac_main()
{
	json::Object obj = json::parse(std::cin,10);
	vector_type v(2);
	v(0) = 1; v(1) = 1;
	vector_type u(2);
	matrix_type J(4,3);

	parameter p(obj);
	jacobi_functor f(p);

	f(v,J,0,u);

	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 3; ++j)
		{
			std::cout << J(i,j) << " ";
		}
		std::cout << std::endl;
	}

	return 0;
}*/

int main()
{
	//jac_main();

	json::Object jin = json::parse(std::cin,10);
	ode_test(jin);
	return 0;
}

