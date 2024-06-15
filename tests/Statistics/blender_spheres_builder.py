
import bpy
import math

### USAGE ###
# use this script on HYD37 with the following command:
# blender --background --python this-script-name.py

# file to read parameteres from
file_name_input = "sphere_data.txt"

# processing parameters
n_package = 1000         # number of spheres in one *.stl file
accuracy_of_stl = 3     # subdivions of icosahedral sphere

# clearing all mesh parts
bpy.ops.object.select_all(action='DESELECT')
mesh = [m for m in bpy.context.scene.objects if m.type == 'MESH']

for obj in mesh:
    obj.select_set(True)
bpy.ops.object.delete()

# initialization of counters and indices
counter = 0; total_counter = 0; i_file = 0

# read data of spheres from the input file
with open ( file_name_input ) as read_file:

    # reading all lines from the file
    all_lines = read_file.readlines()

    # extracting the number of spheres
    first_line = all_lines[0].split()
    n_spheres = int( first_line[3] )

    # some information
    print( "Number of spheres =", n_spheres )
    print( "Number per packages =", n_package )
    print( "Number of STL files =", math.ceil( n_spheres / n_package ) )
    if ( math.ceil( n_spheres / n_package ) == 1 ):
        single_file = True
    else:
        single_file = False

    # iteration over all spheres
    for i_sphere in range( 0, n_spheres ):

        # indices of the lines
        i1_line = 2 * i_sphere + 1
        i2_line = 2 * i_sphere + 2

        # validity check for the line
        if ( "SPHERE" in all_lines[i1_line] ):

            row = all_lines[i2_line].split()

            # get position and check the validity of the sphere
            r = float(row[0])  # (think about factor 0.99 to restore porosity)
            x = float(row[1])
            y = float(row[2])
            z = float(row[3])

            # increasing counter for valid spheres and create mesh
            counter = counter + 1
            total_counter = total_counter + 1
            sphere = bpy.ops.mesh.primitive_ico_sphere_add( subdivisions=accuracy_of_stl, location=(x,y,z), radius=r )

            # giving progress information
            if single_file:
                print( "- ", total_counter, " / ", n_spheres, " generated")

            # check for write point (full package or last sphere)
            if ( (counter == n_package) or (total_counter == n_spheres) ):

                # joining all spheres
                bpy.ops.object.select_all(action='DESELECT')
                mesh = [m for m in bpy.context.scene.objects if m.type == 'MESH']
                for obj in mesh:
                    obj.select_set(True)
                bpy.ops.object.join()

                # rename to "AllSpheres"
                for obj in bpy.context.selected_objects:
                    obj.name = "AllSpheres"

                # write to STL file
                i_file = i_file + 1
                path = 'spheres.stl'
                bpy.ops.export_mesh.stl( filepath=path )

                # clearing all mesh parts
                bpy.ops.object.select_all(action='DESELECT')
                mesh = [m for m in bpy.context.scene.objects if m.type == 'MESH']

                for obj in mesh:
                    obj.select_set(True)
                bpy.ops.object.delete()

                # resetting the counter
                counter = 0
                print( "- number of spheres generated =", total_counter )

        else:

            print( "INVALID LINE found (not a sphere)" )

print( "- number of spheres generated = ", total_counter )
print( "- number of stl files created = ", i_file )
print( "- check Paraview for correctness" )
