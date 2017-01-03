// This files has been copied from https://github.com/google/benchmark/cmake/
// Distributed under the Apache License Version 2.0.
//   https://github.com/google/benchmark/LICENSE
// Modified on Dec. 22th 2015 by Jean-Gabriel Young <jean.gabriel.young@gmail.com>
#include <chrono>

int main() {
    typedef std::chrono::steady_clock Clock;
    Clock::time_point tp = Clock::now();
    ((void)tp);
}