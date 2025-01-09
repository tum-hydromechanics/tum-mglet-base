# # import conduit
# from conduit import blueprint
# n = Node(None)

# # n["my"]="data"

# # print(n)
# n2 = Node(None)
# n2 = mesh.examples.basic("polyhedra",3,3,3,n2)


import conduit
from conduit.relay import io
from conduit.flow import Node

n2 = Node(None)

mesh = mesh.examples.basic("polyhedra",3,3,3,n2)
output_file = "mesh_data.hdf5"
conduit.relay.io.save(mesh, output_file)