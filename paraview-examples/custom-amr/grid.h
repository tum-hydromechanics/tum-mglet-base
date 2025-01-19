#pragma once

#include <array>
#include <vector>

constexpr unsigned WORLD_DIMS = 3;

class Grid {
 public:
  Grid(std::array<unsigned, WORLD_DIMS> dims,
       std::array<float, WORLD_DIMS> spacing,
       std::array<float, WORLD_DIMS> origin,
       unsigned level,
       unsigned id,
       int parent_id,
       std::vector<unsigned> child_ids,
       float default_cell_value = 1.0f);
  auto get_num_values() const -> unsigned;


  std::array<unsigned, WORLD_DIMS> dims;
  std::array<float, WORLD_DIMS> spacing;
  std::array<float, WORLD_DIMS> origin;
  std::vector<float> values;
  unsigned level;
  unsigned id;
  int parent_id;
  std::vector<unsigned> child_ids;
};