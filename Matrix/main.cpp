#include "Matrix.hpp"
#include <random>

#include <time.h>

int main(int argc, char** argv)
{
	std::default_random_engine engine(time(NULL));
	std::uniform_int_distribution<int> range(0, 4);

	Matrix myMatrix(4, 4);
	for (int row = 0; row < 3; row++)
		for (int col = 0; col < 4; col++)
			myMatrix.SetNumber(row, col, range(engine));

	Matrix otherMatrix(4, 4);
	Identity(otherMatrix);

	Matrix addMatrix = myMatrix * otherMatrix;

	std::cout << myMatrix;
	std::cout << std::endl << "-" << std::endl << std::endl;
	std::cout << otherMatrix;
	std::cout << std::endl << "=" << std::endl << std::endl;
	std::cout << addMatrix << std::endl;
	std::cout << det(addMatrix) << std::endl;
	std::cout << addMatrix.IsInvertable() << std::endl;
	std::cout << std::endl;

	getchar();

	return 0;
}