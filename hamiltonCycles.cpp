#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

// Класс для удобного взаимодействия с графом
// Доступен ввод и вывод (файл, поток ввода/вывода)
class Graph {
private:
  vector<vector<int>> adjMatrix;
  int numVertices;

public:
  Graph() {}

  Graph(int numVertices) {
    this->numVertices = numVertices;
    adjMatrix.resize(numVertices, vector<int>(numVertices, 0));
  }

  Graph(vector<vector<int>> matrix) {
    numVertices = matrix.size();
    adjMatrix = matrix;
  }

  // Доступ к матрице смежности
  auto getAdjMatrix() { return &adjMatrix; }

  // Доступ к количеству вершин в графе
  int getNumVertices() { return numVertices; }

  // Установить новое значение веса грани между существующими вершинами
  void changeEdge(int i, int j, int weight) {
    if (i < numVertices && j < numVertices)
      adjMatrix[i][j] = weight;
  }

  // Добавить вершину к уже существующему графу, без ребер
  void addNode() {
    for (int i = 0; i < numVertices; i += 1)
      adjMatrix[i].push_back(0);
    numVertices += 1;
    adjMatrix.push_back(vector<int>(numVertices, 0));
  }

  // Убрать ребро(вес = 0)
  void removeEdge(int i, int j) { adjMatrix[i][j] = adjMatrix[j][i] = false; }

  // Доступ к весу любого ребра между вершинами
  int getEdge(int i, int j) { return adjMatrix[i][j]; }

  Graph &operator=(const Graph &g) {
    numVertices = g.numVertices;
    adjMatrix = g.adjMatrix;
    return *this;
  }

  friend istream &operator>>(istream &ustream, Graph &g);

  friend ostream &operator<<(ostream &ustream, const Graph &g);
};

istream &operator>>(istream &ustream, Graph &g) {
  ustream >> g.numVertices;
  g.adjMatrix.resize(g.numVertices, vector<int>(g.numVertices, 0));
  for (auto &i : g.adjMatrix)
    for (auto &j : i)
      ustream >> j;
  return ustream;
}

ostream &operator<<(ostream &ustream, const Graph &g) {
  if (typeid(ustream).name() == typeid(ofstream).name()) {
    ustream << g.numVertices << '\n';
  } else {
    ustream << "\nVertexCount: " << g.numVertices << "\n";
  }

  int spacesMax = 1;
  for (auto i : g.adjMatrix)
    for (auto j : i)
      spacesMax = max(spacesMax, (int)to_string(j).size() + 1);

  for (auto i : g.adjMatrix) {
    for (auto j : i) {
      ustream << j;
      for (int spaces = spacesMax - to_string(j).size(); spaces > 0; --spaces)
        ustream << ' ';
    }
    ustream << '\n';
  }
  return ustream;
}

// Общий класс алгоритма
// Доступен вывод в поток вывода
class HamiltonAlgorithm {
protected:
  int numVertices;
  vector<vector<int>> adjMatrix;
  vector<int> path;

public:
  HamiltonAlgorithm(Graph &g) {
    numVertices = g.getNumVertices();
    adjMatrix = *g.getAdjMatrix();
  }

  friend ostream &operator<<(ostream &, const HamiltonAlgorithm &HA);
};

ostream &operator<<(ostream &, const HamiltonAlgorithm &HA) {
  for (int p = 0; p < HA.numVertices; p += 1)
    cout << HA.path[p] << " ";
  cout << HA.path[0];
  cout << endl;
  return cout;
}

// Алгоритм перебора Робертса Флореса
// Доступен вывод в поток вывода
class RobertsFloresClass : public HamiltonAlgorithm {
public:
  vector<vector<int>> allHamiltonCycles;

  bool isValid(int v, int length) {
    if (!adjMatrix[path[length - 1]][v])
      return false;
    for (auto _v : path)
      if (_v == v)
        return false;
    return true;
  }

  // Вспомогательный рекурсивый метод нахождения Гамилтоновых циклов
  void findHamiltonCycleRecursive(int length = 1, bool findOne = false) {
    if (length == numVertices && adjMatrix[path[path.size() - 1]][path[0]]) {
      path.push_back(path[0]);
      allHamiltonCycles.push_back(path);
      if (findOne)
        return;
      path.pop_back();
    }
    for (int v = 0; v < numVertices; v++) {
      if (isValid(v, length)) {
        path.push_back(v);
        findHamiltonCycleRecursive(length + 1, findOne);
        path.pop_back();
      }
    }
  }

public:
  RobertsFloresClass(Graph &g) : HamiltonAlgorithm(g) {}

  // Доступ ко всем найденным Гамильтоновым циклам
  vector<vector<int>> getHamiltonCycles() { return allHamiltonCycles; }

  // Основной метод поиска Гамильтоновых циклов
  void findHamiltonCycles(int startVertex = 0, bool findOne = false) {
    path.push_back(startVertex);
    findHamiltonCycleRecursive(1, findOne);
  }

  void algTest() { findHamiltonCycles(0, true); }

  friend ostream &operator<<(ostream &, const RobertsFloresClass &RFC);
};

ostream &operator<<(ostream &, const RobertsFloresClass &RFC) {
  cout << '\n';
  for (auto cycle : RFC.allHamiltonCycles) {
    for (auto v : cycle)
      cout << v << ' ';
    cout << '\n';
  }
  return cout;
}

// Решение задачи коммивояжера динамическим методом
class TSPDynamic : public HamiltonAlgorithm {
private:
  int MAX = numeric_limits<int>::max();
  vector<vector<int>> dp;
  int VISITED_ALL = pow(2, numVertices) - 1;

public:
  TSPDynamic(Graph &g) : HamiltonAlgorithm(g) {
    dp = vector(pow(2, numVertices), vector(numVertices, -1));
  }

  // Рекурсивный метод решения задачи Коммивояжера. Mask представляет собой
  // набор городов. Например, если маска равна 01101, то города номер 1 и 4 еще
  // не были посещены.
  int getMinRoute(int mask = 1, int pos = 0) {
    if (mask == VISITED_ALL) // Все ли города посещены
      return adjMatrix[pos][0];
    if (dp[mask][pos] != -1) { // Вычислено ли уже значение
      return dp[mask][pos];
    }

    int ans = MAX;

    // Посещает все города, которые еще не посещены, и выбирает лучший путь
    for (int city = 0; city < numVertices; city++) {
      if ((mask & (1 << city)) == 0) { // Посещен ли текущий город
        int newAns =
            adjMatrix[pos][city] + getMinRoute(mask | (1 << city), city);
        ans = min(ans, newAns);
      }
    }
    return dp[mask][pos] = ans;
  }

  void algTest() { getMinRoute(); }

  friend ostream &operator<<(ostream &ustream, const TSPDynamic &tsp);
};

ostream &operator<<(ostream &ustream, const TSPDynamic &tsp) {
  cout << "Min-route: " << tsp.dp[1][0] << '\n';
  return ustream;
}

// Вспомогательная функция для подсчета времени работы алгоритма
template <class T> double getTestTime(T &alg) {
  clock_t time;
  time = clock();
  alg.algTest();
  cout << alg << "Time executed: ";
  return ((double)(clock() - time)) / CLOCKS_PER_SEC;
}

int main() {
  for (int num = 2; num <= 12; num++) {
    string path = "../tests/test" + to_string(num) + ".txt";
    cout << path;
    ifstream curTest(path);

    if (curTest.is_open()) {
      Graph testGraph;
      curTest >> testGraph;
      cout << testGraph;

      RobertsFloresClass rfc(testGraph);
      cout << getTestTime(rfc) << "\n\n";

      TSPDynamic tsp(testGraph);
      cout << getTestTime(tsp) << "\n\n";

    } else {
      cout << "The file \"../tests/test" + to_string(num) +
                  ".txt\" not found!\n";
    }
  }
  return 0;
}