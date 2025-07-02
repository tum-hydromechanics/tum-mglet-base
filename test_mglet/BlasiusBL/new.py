# import all necessary modules
from mgtools import xdmfWriterKMT8
 
# define input/output filenames
fieldfile = "fields.h5"
gridfile = "grids.h5"
 
# initialise writer instance
writer = xdmfWriterKMT8(fieldfile, gridfile)
 
# write XDMF-files for each field
writer.writeXdmf(fieldfile, "U")
writer.writeXdmf(fieldfile, "V")
writer.writeXdmf(fieldfile, "W")
writer.writeXdmf(fieldfile, "P")
 
# close writer
writer.close()