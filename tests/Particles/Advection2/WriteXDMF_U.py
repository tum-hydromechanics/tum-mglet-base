from mgtools import xdmfWriterKMT8
writer = xdmfWriterKMT8("fields.h5","particle_grids.h5")
writer.writeXdmf("velocity_u.xdmf","U")
writer.close()