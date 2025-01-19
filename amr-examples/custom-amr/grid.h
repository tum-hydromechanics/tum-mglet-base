#pragma once

#include <array>
#include <vector>

constexpr unsigned N_D = 3;

class Grid {
 public:
  Grid(std::array<unsigned, N_D> cells,
       std::array<float, N_D> spacing,
       std::array<float, N_D> origin,
       unsigned level,
       unsigned id,
       int parent_id,
       std::vector<unsigned> child_ids,
       float default_cell_value = 1.0f);
  auto get_num_values() const -> unsigned;
  auto get_dims() const -> std::array<unsigned, N_D>;
  auto get_spacing() const -> std::array<float, N_D>;
  auto get_global_origin() const -> std::array<float, N_D>;
  auto get_values() const -> std::vector<float>;
  auto get_level() const -> unsigned;
  auto get_id() const -> unsigned;
  auto get_parent_id() const -> int;
  auto get_child_ids() const -> std::vector<unsigned>;

private:
  std::array<unsigned, N_D> cells;
  std::array<float, N_D> spacing;
  std::array<float, N_D> origin;
  std::vector<float> values;
  unsigned level;
  unsigned id;
  int parent_id;
  std::vector<unsigned> child_ids;
};