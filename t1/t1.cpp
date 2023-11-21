#include <iostream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <cmath>

using namespace std;

// Функція для обчислення середнього арифметичного
vector<double> calculateArithmeticMean(const vector<vector<double>>&pairwiseComparisonMatrix) {
	int n = pairwiseComparisonMatrix.size();
	vector<double> priorityVector(n, 0.0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			priorityVector[i] += pairwiseComparisonMatrix[i][j];
		}
		priorityVector[i] /= n;
	}
	return priorityVector;
}

// Функція для обчислення середнього геометричного
vector<double> calculateGeometricMean(const vector<vector<double>>&pairwiseComparisonMatrix) {
	int n = pairwiseComparisonMatrix.size();
	vector<double> priorityVector(n, 1.0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			priorityVector[i] *= pairwiseComparisonMatrix[i][j];
		}
		priorityVector[i] = pow(priorityVector[i], 1.0 / n);
	}
	return priorityVector;
}

void random(vector<vector<double>>& random_vector, int size) {
	for (int i = 0; i < size; i++) {
		vector<double> temp;
		for (int j = 0; j < size; j++) {
			double random_value = static_cast<double>(rand()) / RAND_MAX;
			temp.push_back(random_value);
		} random_vector.push_back(temp);
	}
}

// Функція для виведення матриці на екран
void printMatrix(const vector<vector<double>>& matrix) {
	int n = matrix.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << matrix[i][j] << " ";
		} cout << endl;
	}
}

vector<double> calculateCombinationMethod(const vector<vector<double>>&
	pairwiseComparisonMatrix, int maxIterations) {
	int n = pairwiseComparisonMatrix.size();
	vector<double> priorityVector(n, 1.0);
	for (int iteration = 0; iteration < maxIterations; iteration++) {
		vector<double> newPriorityVector(n, 0.0);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				newPriorityVector[i] += pairwiseComparisonMatrix[i][j] * priorityVector[j];
			}
		}
		double sum = 0.0;
		for (int i = 0; i < n; i++) {
			sum += newPriorityVector[i];
		}
		for (int i = 0; i < n; i++) {
			newPriorityVector[i] /= sum;
		}
		priorityVector = newPriorityVector;
	}

	return priorityVector;
}

int main() {
	srand(time(NULL));
	int size = 10;
	vector<vector<double>> pairwiseComparisonMatrix;
	random(pairwiseComparisonMatrix, size);
	int maxIterations = 10;
	cout << "Pairwise Comparison Matrix:" << endl;
	printMatrix(pairwiseComparisonMatrix);
	vector<double> arithmeticPriorityVector = calculateArithmeticMean(pairwiseComparisonMatrix);
	vector<double> geometricPriorityVector = calculateGeometricMean(pairwiseComparisonMatrix);
	vector<double> combinationPriorityVector =
		calculateCombinationMethod(pairwiseComparisonMatrix, maxIterations);
	cout << "Arithmetic Priority Vector:" << endl;
	for (double value : arithmeticPriorityVector) {
		cout << value << " ";
	} cout << endl;
	cout << "Geometric Priority Vector:" << endl;
	for (double value : geometricPriorityVector) {
		cout << value << " ";
	} cout << endl;
	cout << "Combination Priority Vector:" << endl;
	for (double value : combinationPriorityVector) {
		cout << value << " ";
	} cout << endl;
	// Перший тест
	vector<vector<double>> diagnosticPairwiseComparisonMatrix1 = { {0.25132, 0.13834},
	{0.427381, 0.932188} };
	vector<double> expectedPriorityVector1 = { 0.154204, 0.845796 }; // Очікуваний результат
	vector<double> calculatedPriorityVector1 =
		calculateCombinationMethod(diagnosticPairwiseComparisonMatrix1, maxIterations);
	bool test1Passed = true;
	for (int i = 0; i < diagnosticPairwiseComparisonMatrix1.size(); i++) {
		if (abs(calculatedPriorityVector1[i] - expectedPriorityVector1[i]) > 0.001) {
			cout << "Test 1 failed!" << endl;
			test1Passed = false;
			break;
		}
	}
	if (test1Passed) {
		cout << "Test 1 passed!" << endl;
	}
	// Другий тест
	vector<vector<double>> diagnosticPairwiseComparisonMatrix2 = { {0.252, 0.134}, {0.481, 0.938}
	};
	vector<double> expectedPriorityVector2 = { 0.154, 0.796 }; // Очікуваний результат
	vector<double> calculatedPriorityVector2 =
		calculateCombinationMethod(diagnosticPairwiseComparisonMatrix2, maxIterations);
	bool test2Passed = true;
	for (int i = 0; i < diagnosticPairwiseComparisonMatrix2.size(); i++) {
		if (abs(calculatedPriorityVector2[i] - expectedPriorityVector2[i]) > 0.001) {
			cout << "Test 2 failed!" << endl;
			test2Passed = false;
			break;
		}
	}
	if (test2Passed) {
		cout << "Test 2 passed!" << endl;
	}

	return 0;
}