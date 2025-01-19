#include "grid.h"

Grid::Grid(std::array<unsigned, N_D> cells, std::array<float, N_D> spacing,
           std::array<float, N_D> origin, unsigned level, unsigned id,
           int parent_id, std::vector<unsigned> child_ids,
           float default_cell_value)
    : cells(cells),
      spacing(spacing),
      origin(origin),
      level(level),
      id(id),
      child_ids(child_ids),
      parent_id(parent_id) {
  values = std::vector<float>(get_num_values(), default_cell_value);
}

auto Grid::get_num_values() const -> unsigned {
  return cells[0] * cells[1] * cells[2];
}

auto Grid::get_dims() const -> std::array<unsigned, N_D> {
  return std::array<unsigned, N_D>{cells[0] + 1, cells[1] + 1, cells[2] + 1};
}

auto Grid::get_spacing() const -> std::array<float, N_D> { return spacing; }

auto Grid::get_global_origin() const -> std::array<float, N_D> {
  return origin;
}

auto Grid::get_values() const -> std::vector<float> { return values; }

auto Grid::get_level() const -> unsigned { return level; }

auto Grid::get_id() const -> unsigned { return id; }

auto Grid::get_parent_id() const -> int { return parent_id; }

auto Grid::get_child_ids() const -> std::vector<unsigned> { return child_ids; }
