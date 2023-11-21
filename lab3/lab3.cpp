#include <iostream>
#include <vector>
#include <set>

std::vector<std::vector<double>> graph; // Граф у вигляді списку суміжності

bool isCoveringTree(const std::set<double>& nodes, int n) {
    // Перевірка, чи утворюють вершини nodes покривне дерево для графа розмірності n
    std::vector<bool> visited(n, false);

    for (double node : nodes) {
        int index = static_cast<int>(node);
        visited[index] = true;
        for (int i = 0; i < n; ++i) {
            if (graph[index][i] > 0) {
                visited[i] = true;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            return false;
        }
    }

    return true;
}

void generateAllCoveringTrees(int n) {
    // Генерація всіх можливих підмножин вершин графа
    for (int mask = 1; mask < (1 << n); ++mask) {
        std::set<double> nodes;
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {
                nodes.insert(i);
            }
        }

        // Перевірка, чи утворюють вершини nodes покривне дерево
        if (isCoveringTree(nodes, n)) {
            // nodes є покривним деревом, можна використовувати його як поточний результат
            std::cout << "Covering Tree: ";
            for (double node : nodes) {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
    }
}

int main() {
    // Приклад графу (задайте свій граф у вигляді списку суміжності)
    graph = {
        {1, 2, 3, 4, 2, 6, 7},
        {0.5, 1, 4, 3, 2, 4, 8},
        {0.33, 0.25, 1, 0, 2, 3, 4},
        {0.25, 0, 3, 1, 5, 7, 4},
        {0.5, 0.5, 0, 0.2, 1, 2, 5},
        {0.17, 0.25, 0.33, 0.14, 0.5, 1, 0},
        {0.14, 0.125, 0.25, 0.25, 0.2, 0.33, 1}
    };
    int n = graph.size(); // Розмірність графа

    std::cout << "All Covering Trees:" << std::endl;
    generateAllCoveringTrees(n);

    return 0;
}
