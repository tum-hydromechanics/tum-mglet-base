#include "catalyst_adaptor.h"


constexpr unsigned ITERATIONS = 1;
constexpr float DT = 0.1;

int main(int argc, char* argv[]) {
    catalyst_adaptor::initialize(argc, argv);
    float time = 0;
    for (unsigned i = 0; i < ITERATIONS; ++i) {
        catalyst_adaptor::execute(i, time);
        time += DT;
    }
    catalyst_adaptor::finalize();
}