#include <iostream>

#include "CatalystAdaptor.h"
#include "grid.h"


constexpr unsigned ITERATIONS = 1;
constexpr float DT = 0.1;

int main(int argc, char* argv[]) {
    // Coarse grid 1
    Grid grid_0 {
        std::array{ 2u, 2u, 2u },
        std::array{ 1.0f, 1.0f, 1.0f },
        std::array{ 2.0f, 0.0f, 0.0f },
        0,
        4,
        -1,
        std::vector<unsigned>{},
        0.75f
    };

    // Coarse grid 2
    Grid grid_1 {
        std::array{ 2u, 2u, 2u },        // Cells per dimension
        std::array{ 1.0f, 1.0f, 1.0f },  // Spacing
        std::array{ 0.0f, 0.0f, 0.0f },  // Origin (global coordinate system)
        0,                               // Level
        0,                               // Grid id
        -1,                              // Parent grid id
        std::vector<unsigned>{1, 2},     // Child grid ids
        0.5f                             // Default value for cells
    };

    // Refined grid 1 in coarse grid 2
    Grid grid_2 {
        std::array{ 2u, 2u, 2u },
        std::array{ 0.5f, 0.5f, 0.5f },
        std::array{ 0.0f, 0.0f, 0.0f },
        1,
        1, 
        0,
        std::vector<unsigned>{},
        1.0f
    };

    // Refined grid 2 in coarse grid 2
    Grid grid_3 {
        std::array{ 2u, 2u, 2u },
        std::array{ 0.5f, 0.5f, 0.5f },
        std::array{ 1.0f, 1.0f, 1.0f },
        1,
        2,
        0,
        std::vector<unsigned>{3},
        1.0f
    };

    // More refined grid 1 in refined grid 2
    Grid grid_4 {
        std::array{ 2u, 2u, 2u },
        std::array{ 0.25f, 0.25f, 0.25f },
        std::array{ 1.5f, 1.5f, 1.5f },
        2,
        3,
        2,
        std::vector<unsigned>{},
        2.0f
    };

    std::vector grids {
        grid_0, grid_1, grid_2, grid_3, grid_4
    };

    catalyst_adaptor::initialize(argc, argv);
    float time = 0;
    for (unsigned i = 0; i < ITERATIONS; ++i) {
        catalyst_adaptor::execute(i, time, grids);
        time += DT;
    }
    catalyst_adaptor::finalize();
}