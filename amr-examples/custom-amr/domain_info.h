#pragma once

#include <vector>
#include "grid.h"

struct LevelInfo {
    
};

class DomainInfo {
public:
  DomainInfo(const std::vector<Grid>& grids);

private:
  std::vector<LevelInfo> level_info;
};