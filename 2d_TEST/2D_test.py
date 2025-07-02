#-------------------------------------------------------------------------------
# Example 2D Grid (thin slab in z) with nlevels=1
#-------------------------------------------------------------------------------

from gridgen3 import (
    Domain, 
    BoundaryCondition,
    GridsFile,
    Precision, 
    UniformGridGenerator, 
    XdmfWriter
)
import math

# ---- Step 1: Define the domain -----------------------------------------------
print("Creating 2D domain as thin slab ... \n")

# x and y size
xfr = 0.0
xba = 16.0
yri = 0.0
yle = 16.0

# z direction: tiny thickness
zbo = 0.0
zto = 1e-8

# Boundary conditions
bfr = [ BoundaryCondition("PER") ]   # x front
bba = [ BoundaryCondition("PER") ]   # x back
bri = [ BoundaryCondition("PER") ]   # y lower
ble = [ BoundaryCondition("PER") ]   # y upper
bbo = [ BoundaryCondition("NOS") ]   # z bottom
bto = [ BoundaryCondition("NOS") ]   # z top

# Create domain
domain = Domain(
    xfr, xba,
    yri, yle,
    zbo, zto,
    bfr, bba,
    bri, ble,
    bbo, bto
)

print("finished.\n")

# ---- Step 2: Generate the base grid ------------------------------------------
print("Creating 2D base grid ... \n")

gridgen = UniformGridGenerator(
    domain=domain,
    shape=(16,16,1),    # z=1
    box_size=(1,1,1),   # z=1
    nlevels=1,            # Only 1 level
    mode='xyz'
)

basegrid = gridgen.generate()

print("finished.\n")

# ---- Step 3: Assign grid ids -------------------------------------------------
ngrid = 1
for level in basegrid:
    ngrid = level.enumerate_grids(ngrid)

# ---- Step 4: Write the HDF5 file ---------------------------------------------
print("Writing output2d_thin.h5 ... ", end="", flush=True)

precisions = Precision("double")

writer2 = GridsFile("output2d_thin.h5", precisions)
writer2.write(basegrid)

print("finished.\n")

# ---- Step 5: Write the XDMF file ---------------------------------------------
print("Writing output2d_thin.xdmf ... ", end="", flush=True)

writer3 = XdmfWriter("output2d_thin.xdmf")
writer3.write(basegrid)

print("finished.\n")
