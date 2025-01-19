#include "grid.h"

Grid::Grid(std::array<unsigned, WORLD_DIMS> dims,
           std::array<float, WORLD_DIMS> spacing,
           std::array<float, WORLD_DIMS> origin, unsigned level, unsigned id,
           int parent_id, std::vector<unsigned> child_ids, float default_cell_value)
    : dims(dims),
      spacing(spacing),
      origin(origin),
      level(level),
      id(id),
      child_ids(child_ids),
      parent_id(parent_id) {
  values = std::vector<float>(get_num_values(), default_cell_value);
}

auto Grid::get_num_values() const -> unsigned {
  return (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
}
