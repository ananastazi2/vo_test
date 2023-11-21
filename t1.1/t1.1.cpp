#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>

using namespace std;

const int V = 7;

bool elementallyAveragingVectorsMethod; //0 - arithmetic mean, 1 - geometric mean

vector<pair<int, int>> graphEdges = {
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6},
        {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6},
        {2, 3}, {2, 4}, {2, 5}, {2, 6},
        {3, 4}, {3, 5}, {3, 6},
        {4, 5}, {4, 6},
        {5, 6}
};

vector <vector <double>>  edgeWeights = {
    {1,    2,     3,    4,    2,   6,    7},
    {0.5,  1,     4,    3,    2,   4,    8},
    {0.33, 0.25,  1,    0.33, 2,   3,    4},
    {0.25, 0.33,  3,    1,    5,   7,    4},
    {0.5,  0.5,   0.5,  0.2,  1,   2,    5},
    {0.17, 0.25,  0.33, 0.14, 0.5, 1,    3},
    {0.14, 0.125, 0.25, 0.25, 0.2, 0.33, 1}
};

vector<vector<pair<int, int>>> spanningTrees;

vector<vector <double>> spanningTreeMatrix(V, vector<double>(V));

vector<vector <double>> priorityVectors;

vector <double> generalPriorityVector(V);

void makeSpanningTreeMatrixReflexive() {
    for (int i = 0; i < spanningTreeMatrix.size(); ++i) {
        for (int j = 0; j < spanningTreeMatrix[0].size(); ++j) {
            if (i == j)
                spanningTreeMatrix[i][j] = 1;
            else
                spanningTreeMatrix[i][j] = 0;
        }
    }
}

bool isSpanningTree(const vector<pair<int, int>>& edges) {
    vector<int> parent(V);
    iota(parent.begin(), parent.end(), 0); 

    // Find root of the set that vertex i is in
    auto find = [&parent](int i) {
        while (i != parent[i]) {
            i = parent[i];
        }
        return i;
    };

    // Union of two sets. Returns false if x and y are already in the same set
    auto unionSets = [&parent, &find](int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX == rootY) {
            return false;
        }
        parent[rootY] = rootX;
        return true;
    };

    // Check if all edges can be added without forming a cycle
    for (const auto& edge : edges) {
        if (!unionSets(edge.first, edge.second)) {
            return false; // Cycle detected
        }
    }

    // Check if all vertices are connected
    int root = find(0);
    for (int i = 1; i < V; ++i) {
        if (find(i) != root) {
            return false; // Not all vertices are connected
        }
    }

    return true;
}

void generateMatrixOfSpanningTree(vector<pair<int, int>>& edgeSubset) {
    makeSpanningTreeMatrixReflexive();

    // Fill the matrix according to the rule of symmetry
    for (pair<int, int> verticesPair : edgeSubset)
    {
        spanningTreeMatrix[verticesPair.first][verticesPair.second] =
            edgeWeights[verticesPair.first][verticesPair.second];

        spanningTreeMatrix[verticesPair.second][verticesPair.first] =
            1.0 / spanningTreeMatrix[verticesPair.first][verticesPair.second];
    }

    //Fill missing elements with a transitive rule
    int numberOfMissingElements = V * (V - 1) - (edgeSubset.size() * 2);
    int maxEdges, numOfEdges, k;

    while (numberOfMissingElements > 0) {
        maxEdges = 0;
        k = 0;

        for (int i = 0; i < spanningTreeMatrix.size(); ++i) {
            numOfEdges = 0;

            for (int j = 0; j < spanningTreeMatrix[0].size(); ++j) {
                if (spanningTreeMatrix[i][j])
                    numOfEdges++;
            }

            if (numOfEdges > maxEdges) {
                maxEdges = numOfEdges;
                k = i;
            }
        }

        for (int i = 0; i < spanningTreeMatrix.size(); ++i) {
            for (int j = 0; j < spanningTreeMatrix[0].size(); ++j) {
                if (i == k || k == j || i == j)
                    continue;
                if (!spanningTreeMatrix[i][j] && 
                    spanningTreeMatrix[i][k] && 
                    spanningTreeMatrix[k][j]) 
                {
                    spanningTreeMatrix[i][j] = spanningTreeMatrix[i][k] * spanningTreeMatrix[k][j];
                    numberOfMissingElements--;
                }
            }
        }
    }
}

vector <double> calculatePriorityVectorForSpanningTree() {
    vector <double> priorityVector(V, 0);

    double sumOfRowelements = 0;

    for (int j = 0; j < spanningTreeMatrix[0].size(); ++j) {
        sumOfRowelements += spanningTreeMatrix[0][j];
    }

    for (int j = 0; j < spanningTreeMatrix[0].size(); ++j) {
        priorityVector[j] = spanningTreeMatrix[0][j] / sumOfRowelements;
    }

    return priorityVector;
}

void ArithmeticMean() {
    double sum;

    for (int j = 0; j < priorityVectors[0].size(); ++j) {
        sum = 0;
        for (int i = 0; i < priorityVectors.size(); ++i) {
            sum += priorityVectors[i][j];
        }

        generalPriorityVector[j] = sum / spanningTrees.size();
    }
}

void GeometricMean() {
    double sum;

    for (int j = 0; j < priorityVectors[0].size(); ++j) {
        sum = 0;
        for (int i = 0; i < priorityVectors.size(); ++i) {
            sum += log(priorityVectors[i][j]);
        }

        generalPriorityVector[j] = exp(sum / priorityVectors.size());
    }
}

void calculateGeneralPriorityVector() {
    if (!elementallyAveragingVectorsMethod)
        ArithmeticMean();
    else GeometricMean();
}

void generateAllSpanningTrees() {
    // Generate all possible subsets of edges of size V-1 (the number of edges in a spanning tree)
    vector<bool> bitmask(graphEdges.size());
    fill(bitmask.end() - (V - 1), bitmask.end(), true);

    do {
        vector<pair<int, int>> edgeSubset;
        for (int i = 0; i < graphEdges.size(); ++i) {
            if (bitmask[i]) {
                edgeSubset.push_back(graphEdges[i]);
            }
        }

        // checking subset for being a spanning tree
        if (isSpanningTree(edgeSubset)) {
            spanningTrees.push_back(edgeSubset);   
        }

    } while (next_permutation(bitmask.begin(), bitmask.end()));

    cout << "Number of trees: " << spanningTrees.size();
}

void printMatrix() {
    for (int i = 0; i < spanningTreeMatrix.size(); ++i) {
        for (int j = 0; j < spanningTreeMatrix[0].size(); ++j) {
           cout << spanningTreeMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

void Algorithm() {
    double startTime = omp_get_wtime();

    generateAllSpanningTrees();

    double endTime = omp_get_wtime();

    double elapsedTimeMilliseconds = (endTime - startTime) * 1000.0;

    cout << endl << "Time spent for generating all spanning trees: " << elapsedTimeMilliseconds;

    startTime = omp_get_wtime();

    for (int i = 0; i < spanningTrees.size(); ++i) {
        generateMatrixOfSpanningTree(spanningTrees[i]);

        //cout << "SpanningTree matrix" << i << ":" << endl;
        //printMatrix();

        priorityVectors.push_back(calculatePriorityVectorForSpanningTree());
        //cout << endl << " PriorityVectorForSpanningTree " << i << " pushed" << endl;
    }

    endTime = omp_get_wtime();

    elapsedTimeMilliseconds = (endTime - startTime) * 1000.0;

    cout << endl << "Time spent for generating all matricies and calculating all priority vectors: " << elapsedTimeMilliseconds;

    calculateGeneralPriorityVector();
}

int main() {

    cout << "Choose your method for averaging vectors elementally (0 - arithmetic mean, 1 - geometric mean): "; cin >> elementallyAveragingVectorsMethod;

    double startTime = omp_get_wtime();

    Algorithm();

    double endTime = omp_get_wtime();

    double elapsedTimeMilliseconds = (endTime - startTime) * 1000.0;
    
    cout << endl << "General priority vector: ";
    for (int i = 0; i < generalPriorityVector.size(); ++i) {
        cout << generalPriorityVector[i] << " ";
    }

    cout << endl << "Time spent total: " << elapsedTimeMilliseconds;

    return 0;
}