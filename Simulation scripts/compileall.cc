#include <iostream>
#include <vector>
#include <string>
using namespace std;

int main() {
    int runs = 16;

    for (int i = 1; i <= runs; i++) {
        string command = string("g++ -std=c++11 run") + to_string(i) + string(".cc -o run") + to_string(i);
        system(command.c_str());
    }

    return 0;
}