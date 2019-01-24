#include "Matrix.hpp"
#include <random>

#include <time.h>

int main(int argc, char** argv)
{
	std::default_random_engine engine(time(NULL));
	std::uniform_int_distribution<int> range(0, 4);

	// Hardcode a matrix
	Matrix myMatrix(3, 3);
	myMatrix.SetNumber(0, 0, 1);
	myMatrix.SetNumber(0, 1, 1);
	myMatrix.SetNumber(0, 2, 1);
	myMatrix.SetNumber(1, 0, 0);
	myMatrix.SetNumber(1, 1, 1);
	myMatrix.SetNumber(1, 2, 1);
	myMatrix.SetNumber(2, 0, 0);
	myMatrix.SetNumber(2, 1, 0);
	myMatrix.SetNumber(2, 2, 1);

	Matrix inverse = myMatrix;
	inverse.Invert();
	

	std::cout << myMatrix << std::endl;

	std::cout << inverse << std::endl;

	std::cout << mul(myMatrix, inverse) << std::endl;

	getchar();

	return 0;
}