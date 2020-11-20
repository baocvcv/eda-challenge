#include "pq.hpp"
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

#define N 20

int main() {
    priority_queue<int> pq;
    for (int i = 0; i < N; i++) {
        int c = rand() % 10000;
        cout << i << ": " << c << ";  ";
        pq.insert(i, c);
    }
    cout << endl << "pq: " << endl;
    for (int i = 0; i < N; i++) {
        cout << pq.vec[i] << ": " << pq.cost[pq.vec[i]] << ";  ";
    }
    cout << endl;

    pq.update(11, 9200);
    pq.update(4, 0);
    cout << endl << "pq: " << endl;
    for (int i = 0; i < N; i++) {
        cout << pq.vec[i] << ": " << pq.cost[pq.vec[i]] << ";  ";
    }
    cout << endl;

    for (int i = 0; i < 10; i++) pq.pop();
    cout << endl << "pq: " << endl;
    for (int i = 0; i < pq.size(); i++) {
        cout << pq.vec[i] << ": " << pq.cost[pq.vec[i]] << ";  ";
    }
    cout << endl;

    return 0;
}