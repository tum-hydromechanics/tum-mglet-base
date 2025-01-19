#include <iostream>

#include "CatalystAdaptor.h"
#include "grid.h"


constexpr unsigned ITERATIONS = 1;
constexpr float DT = 0.1;

int main(int argc, char* argv[]) {
    float time = 0;

    Grid grid_5 {
        std::array{ 3u, 3u, 3u },
        std::array{ 1.0f, 1.0f, 1.0f },
        std::array{ 2.0f, 0.0f, 0.0f },
        0, // Level
        4, // Id
        -1, // Parent Id
        std::vector<unsigned>{}, // Child ids
        0.75f // Val
    };

    // The parent of all grids (this grid wrapps everything)
    Grid grid_1 {
        std::array{ 3u, 3u, 3u },
        std::array{ 1.0f, 1.0f, 1.0f },
        std::array{ 0.0f, 0.0f, 0.0f },
        0, // Level
        0, // Id
        -1, // Parent Id
        std::vector<unsigned>{1, 2}, // Child ids
        0.5f // Val
    };

    Grid grid_2 {
        std::array{ 3u, 3u, 3u },
        std::array{ 0.5f, 0.5f, 0.5f },
        std::array{ 0.0f, 0.0f, 0.0f },
        1, // Level
        1, // Id
        0, // Parent Id
        std::vector<unsigned>{}, // Child ids
        1.0f // Val
    };

    Grid grid_3 {
        std::array{ 3u, 3u, 3u },
        std::array{ 0.5f, 0.5f, 0.5f },
        std::array{ 1.0f, 1.0f, 1.0f },
        1, // Level
        2, // Id
        0, // Parent Id
        std::vector<unsigned>{3}, // Child ids
        1.0f // Val
    };

    Grid grid_4 {
        std::array{ 3u, 3u, 3u },
        std::array{ 0.25f, 0.25f, 0.25f },
        std::array{ 1.5f, 1.5f, 1.5f },
        2, // Level
        3, // Id
        2, // Parent Id
        std::vector<unsigned>{}, // Child ids
        2.0f // Val
    };


    std::vector grids {
        grid_1, grid_2, grid_3, grid_4, grid_5
    };

    catalyst_adaptor::initialize(argc, argv);
    for (unsigned i = 0; i < ITERATIONS; ++i) {
        catalyst_adaptor::execute(i, time, grids);
        time += DT;
    }
    catalyst_adaptor::finalize();
}