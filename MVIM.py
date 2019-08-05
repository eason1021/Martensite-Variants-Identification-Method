from ovito.io import import_file
from ovito.io import export_file
from ovito.modifiers import *
from ovito.data import *
from ovito.data import NearestNeighborFinder
import numpy as np
import sys
import scipy.linalg as sl
import math

# input file 
node = import_file('E:/MD_PostProcess/data/10nm_tensile/average.*')
node.add_to_scene()
#input reference
node_ref = import_file('E:/MD_PostProcess/data/10nm_tensile/ref')
#create directory to store dynamic variable
my_var = {}

#calculate score

tol = 0.012
ntol = -0.012

#define variants

#X direction
M1 = np.array([0,1,1,1,0,0,0,0,1,1,1,0,1,0,0])
M2 = np.array([0,1,1,1,0,0,1,1,1,0,0,0,1,0,0])
M3 = np.array([0,1,1,1,0,0,1,0,0,0,1,1,1,0,0])
M4 = np.array([0,1,1,1,0,0,0,1,0,1,0,1,1,0,0])
O1 = np.array([0,1,1,1,0,0,0,0,1,0,0,0,1,0,0])
O2 = np.array([0,1,1,1,0,0,0,0,0,0,0,1,1,0,0])
A  = np.array([0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])

#Y direction
M5 = np.array([1,0,1,0,1,0,0,1,0,1,0,1,0,1,0])
M6 = np.array([1,0,1,0,1,0,1,1,1,0,0,0,0,1,0])
M7 = np.array([1,0,1,0,1,0,1,0,0,0,1,1,0,1,0])
M8 = np.array([1,0,1,0,1,0,0,0,1,1,1,0,0,1,0])
O3 = np.array([1,0,1,0,1,0,0,1,0,0,0,0,0,1,0])
O4 = np.array([1,0,1,0,1,0,0,0,0,0,1,0,0,1,0])

#Z direction
M9 = np.array([1,1,0,0,0,1,1,0,0,0,1,1,0,0,1])
M10= np.array([1,1,0,0,0,1,1,1,1,0,0,0,0,0,1])
M11= np.array([1,1,0,0,0,1,0,0,1,1,1,0,0,0,1])
M12= np.array([1,1,0,0,0,1,0,1,0,1,0,1,0,0,1])
O5 = np.array([1,1,0,0,0,1,1,0,0,0,0,0,0,0,1])
O6 = np.array([1,1,0,0,0,1,0,0,0,1,0,0,0,0,1])

#R phase


#DECIDE REF PBC
cell_ref = node_ref.compute().cell_
with cell_ref:
    cell_ref.pbc=(True,True,False)

#calculate ref
pbc_x_ref = float(cell_ref[0][0])
pbc_y_ref = float(cell_ref[1][1])
pbc_z_ref = float(cell_ref[2][2])
data_ref = node_ref.compute()

#Found REF nearest 6 atoms around the center atom
N_ref = 6
finder_ref = NearestNeighborFinder(N_ref, data_ref)
ptypes_ref = data_ref.particles['Particle Type']
positions_ref = data_ref.particles['Position']

# Loop over all input particles:
for index_ref in range(data_ref.particles.count):

    neighbors_ref = [ (neigh_ref.index, neigh_ref.delta) for neigh_ref in finder_ref.find(index_ref) ]

    #rearrange  
    resorted_neighbors_x_ref = sorted( neighbors_ref , key=lambda k: [k[1][0], k[1][1], k[1][2]], reverse=True )
    resorted_neighbors_y_ref = sorted( neighbors_ref , key=lambda k: [k[1][1], k[1][0], k[1][2]], reverse=True )
    resorted_neighbors_z_ref = sorted( neighbors_ref , key=lambda k: [k[1][2], k[1][0], k[1][1]], reverse=True )
    #X Y Z -X -Y -Z
    my_var["list_%s"%index_ref+"_x_ref"] = [n_ref[0] for n_ref in resorted_neighbors_x_ref]
    #Y X Z -Y -X -Z
    my_var["list_%s"%index_ref+"_y_ref"] = [n_ref[0] for n_ref in resorted_neighbors_y_ref]
    #Z X Y -Z -X -Y
    my_var["list_%s"%index_ref+"_z_ref"] = [n_ref[0] for n_ref in resorted_neighbors_z_ref]
    
for i_ref in range(data_ref.particles.count):

    #STORE REMOVE ATOM
    my_var["make_remove_px__ref%s"%i_ref] = []
    my_var["make_remove_py__ref%s"%i_ref] = []
    my_var["make_remove_pz__ref%s"%i_ref] = []
    my_var["make_remove_nx__ref%s"%i_ref] = []
    my_var["make_remove_ny__ref%s"%i_ref] = []
    my_var["make_remove_nz__ref%s"%i_ref] = []
    
    #X -X
    my_var["make_list_%s"%i_ref+"_x_ref"] = my_var["list_%s"%i_ref+"_x_ref"].copy()
    my_var["make_remove_px__ref%s"%i_ref].insert(0,my_var["make_list_%s"%i_ref+"_x_ref"].pop(0))
    my_var["make_remove_nx__ref%s"%i_ref].insert(0,my_var["make_list_%s"%i_ref+"_x_ref"].pop(4))
    
    #Y -Y    
    my_var["make_list_%s"%i_ref+"_y_ref"] = my_var["list_%s"%i_ref+"_y_ref"].copy()
    my_var["make_remove_py__ref%s"%i_ref].insert(0,my_var["make_list_%s"%i_ref+"_y_ref"].pop(0))
    my_var["make_remove_ny__ref%s"%i_ref].insert(0,my_var["make_list_%s"%i_ref+"_y_ref"].pop(4))

    #Z -Z  
    my_var["make_list_%s"%i_ref+"_z_ref"] = my_var["list_%s"%i_ref+"_z_ref"].copy()
    my_var["make_remove_pz__ref%s"%i_ref].insert(0,my_var["make_list_%s"%i_ref+"_z_ref"].pop(0))
    my_var["make_remove_nz__ref%s"%i_ref].insert(0,my_var["make_list_%s"%i_ref+"_z_ref"].pop(4))
        
    #Arrange like id x y z -x -y -z
    my_var["neigh_list__ref%s"%i_ref] = []
    my_var["neigh_list__ref%s"%i_ref] = my_var["make_remove_px__ref%s"%i_ref] + my_var["make_remove_py__ref%s"%i_ref] + my_var["make_remove_pz__ref%s"%i_ref] + my_var["make_remove_nx__ref%s"%i_ref] + my_var["make_remove_ny__ref%s"%i_ref] + my_var["make_remove_nz__ref%s"%i_ref]
    my_var["neigh_list__ref%s"%i_ref].insert(0,i_ref)

for j_ref in range(data_ref.particles.count):
    
    # +X Direction
    my_var["list_%s"%j_ref+"_px_ref"] = my_var["neigh_list__ref%s"%j_ref].copy()
    remove_px_ref = my_var["list_%s"%j_ref+"_px_ref"].pop(1)
    remove_nx_ref = my_var["list_%s"%j_ref+"_px_ref"].pop(3)
    # id y z -y -z

    #find next atom along +x direction
    my_var["list_%s"%j_ref+"_npx_ref"] = my_var["neigh_list__ref%s"%remove_px_ref].copy()
    remove_cen_px_ref = my_var["list_%s"%j_ref+"_npx_ref"].pop(0)
    remove_npx_ref = my_var["list_%s"%j_ref+"_npx_ref"].pop(0)
    remove_nnx_ref = my_var["list_%s"%j_ref+"_npx_ref"].pop(2)
    #+x direction neighbor atom y z -y -z
        
    my_var["tetragonal_id_%s"%j_ref+"_px_ref"] = my_var["list_%s"%j_ref+"_px_ref"].copy() + my_var["list_%s"%j_ref+"_npx_ref"].copy() 
    # ID +Y +Z -Y -Z +Y +Z -Y -Z

    #+Y DIRECTION
    my_var["list_%s"%j_ref+"_py_ref"] = my_var["neigh_list__ref%s"%j_ref].copy()
    remove_py_ref = my_var["list_%s"%j_ref+"_py_ref"].pop(2)
    remove_ny_ref = my_var["list_%s"%j_ref+"_py_ref"].pop(4)
    #id x z -x -z

    #find next atom along +y direction
    my_var["list_%s"%j_ref+"_npy_ref"] = my_var["neigh_list__ref%s"%remove_py_ref].copy()
    remove_cen_py_ref = my_var["list_%s"%j_ref+"_npy_ref"].pop(0)
    remove_npy_ref = my_var["list_%s"%j_ref+"_npy_ref"].pop(1)
    remove_nny_ref = my_var["list_%s"%j_ref+"_npy_ref"].pop(3)
    #+y direction neighbor atom x z -x -z
        
    my_var["tetragonal_id_%s"%j_ref+"_py_ref"] = my_var["list_%s"%j_ref+"_py_ref"].copy() + my_var["list_%s"%j_ref+"_npy_ref"].copy()
    # ID +X +Z -X -Z +X +Z -X -Z
        
    #+z direction
    my_var["list_%s"%j_ref+"_pz_ref"] = my_var["neigh_list__ref%s"%j_ref].copy()
    remove_pz_ref = my_var["list_%s"%j_ref+"_pz_ref"].pop(3)
    remove_nz_ref = my_var["list_%s"%j_ref+"_pz_ref"].pop(5)
    # id x y -x -y

    #find next atom along +z direction
    my_var["list_%s"%j_ref+"_npz_ref"] = my_var["neigh_list__ref%s"%remove_pz_ref].copy()
    remove_cen_pz_ref = my_var["list_%s"%j_ref+"_npz_ref"].pop(0)
    remove_npz_ref = my_var["list_%s"%j_ref+"_npz_ref"].pop(2)
    remove_nnz_ref = my_var["list_%s"%j_ref+"_npz_ref"].pop(4)
    # +z direction x y -x -y

    my_var["tetragonal_id_%s"%j_ref+"_pz_ref"] = my_var["list_%s"%j_ref+"_pz_ref"].copy() + my_var["list_%s"%j_ref+"_npz_ref"].copy()
    # ID +X +Y -X -Y +X +Y -X -Y

    #find R phase
    #my_var["list_%s"%j_ref+"_r_ref"] = my_var["neigh_list__ref%s"%j_ref].copy()
    #remove_r1_ref = my_var["list_%s"%j_ref+"_r_ref"].pop(4)
    #remove_r2_ref = my_var["list_%s"%j_ref+"_r_ref"].pop(4)
    #remove_r3_ref = my_var["list_%s"%j_ref+"_r_ref"].pop(4)
    #my_var["tetragonal_id_%s"%j_ref+"_r_ref"] = my_var["list_%s"%j_ref+"_r_ref"].copy()



for k_ref in range(data_ref.particles.count):
        
    #+x direciton three vector
    #take +Y -Y +Z NEXT +Z
    my_var["tetragonal_id_%s_copy"%k_ref+"_px_ref"] = my_var["tetragonal_id_%s"%k_ref+"_px_ref"].copy()
        
    #catch atom +y -y +z next+z
    my_var["px_y__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_px_ref"][1]
    my_var["px_ny__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_px_ref"][3]
    my_var["px_z__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_px_ref"][2]
    my_var["px_next_z__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_px_ref"][6]
        
    #take atom position
    my_var["px_y_pos__ref%s"%k_ref] = positions_ref[my_var["px_y__ref%s"%k_ref]]
    my_var["px_ny_pos__ref%s"%k_ref] = positions_ref[my_var["px_ny__ref%s"%k_ref]]
    my_var["px_z_pos__ref%s"%k_ref] = positions_ref[my_var["px_z__ref%s"%k_ref]]
    my_var["px_next_z_pos__ref%s"%k_ref] = positions_ref[my_var["px_next_z__ref%s"%k_ref]]#STANDARD POINT
        
    #take three vector
    my_var["px_y_vector__ref%s"%k_ref] = my_var["px_y_pos__ref%s"%k_ref] - my_var["px_z_pos__ref%s"%k_ref]

    if(my_var["px_y_vector__ref%s"%k_ref][0] - pbc_x_ref*0.9 >= 0 ):
        my_var["px_y_vector__ref%s"%k_ref][0] = my_var["px_y_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["px_y_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["px_y_vector__ref%s"%k_ref][0] = my_var["px_y_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["px_y_vector__ref%s"%k_ref][1] - pbc_y_ref*0.9 >= 0 ):
        my_var["px_y_vector__ref%s"%k_ref][1] = my_var["px_y_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["px_y_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["px_y_vector__ref%s"%k_ref][1] = my_var["px_y_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["px_y_vector__ref%s"%k_ref][2] - pbc_z_ref*0.9 >= 0 ):
        my_var["px_y_vector__ref%s"%k_ref][2] = my_var["px_y_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["px_y_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["px_y_vector__ref%s"%k_ref][2] = my_var["px_y_vector__ref%s"%k_ref][2] + pbc_z_ref

    my_var["px_ny_vector__ref%s"%k_ref] = my_var["px_ny_pos__ref%s"%k_ref] - my_var["px_z_pos__ref%s"%k_ref]

    if(my_var["px_ny_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["px_ny_vector__ref%s"%k_ref][0] = my_var["px_ny_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["px_ny_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["px_ny_vector__ref%s"%k_ref][0] = my_var["px_ny_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["px_ny_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["px_ny_vector__ref%s"%k_ref][1] = my_var["px_ny_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["px_ny_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["px_ny_vector__ref%s"%k_ref][1] = my_var["px_ny_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["px_ny_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["px_ny_vector__ref%s"%k_ref][2] = my_var["px_ny_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["px_ny_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["px_ny_vector__ref%s"%k_ref][2] = my_var["px_ny_vector__ref%s"%k_ref][2] + pbc_z_ref

    my_var["px_z_vector__ref%s"%k_ref] = my_var["px_next_z_pos__ref%s"%k_ref] - my_var["px_z_pos__ref%s"%k_ref]

    if(my_var["px_z_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["px_z_vector__ref%s"%k_ref][0] = my_var["px_z_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["px_z_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["px_z_vector__ref%s"%k_ref][0] = my_var["px_z_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["px_z_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["px_z_vector__ref%s"%k_ref][1] = my_var["px_z_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["px_z_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["px_z_vector__ref%s"%k_ref][1] = my_var["px_z_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["px_z_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["px_z_vector__ref%s"%k_ref][2] = my_var["px_z_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["px_z_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["px_z_vector__ref%s"%k_ref][2] = my_var["px_z_vector__ref%s"%k_ref][2] + pbc_z_ref

    #transfer vector to matrix form
    my_var["px_y_matrix__ref%s"%k_ref] = np.array(my_var["px_y_vector__ref%s"%k_ref])
    my_var["px_ny_matrix__ref%s"%k_ref] = np.array(my_var["px_ny_vector__ref%s"%k_ref])
    my_var["px_z_matrix__ref%s"%k_ref] = np.array(my_var["px_z_vector__ref%s"%k_ref])
    
    #ref matrix
    my_var["px_matrix__ref%s"%k_ref] = np.column_stack([my_var["px_ny_matrix__ref%s"%k_ref],my_var["px_y_matrix__ref%s"%k_ref],my_var["px_z_matrix__ref%s"%k_ref]])
    
    #ref matrix inverse 
    my_var["px_matrix_inv_ref%s"%k_ref] = np.linalg.pinv(my_var["px_matrix__ref%s"%k_ref])

    #+y direciton three vector
    #take +X -X +Z NEXT +Z    
    my_var["tetragonal_id_%s_copy"%k_ref+"_py_ref"] = my_var["tetragonal_id_%s"%k_ref+"_py_ref"].copy()
        
    #take id +X -X +Z NEXT +Z 
    my_var["py_x__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_py_ref"][1]
    my_var["py_nx__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_py_ref"][3]
    my_var["py_z__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_py_ref"][2]
    my_var["py_next_z__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_py_ref"][6]
        
    #take position
    my_var["py_x_pos__ref%s"%k_ref] = positions_ref[my_var["py_x__ref%s"%k_ref]]
    my_var["py_nx_pos__ref%s"%k_ref] = positions_ref[my_var["py_nx__ref%s"%k_ref]]
    my_var["py_z_pos__ref%s"%k_ref] = positions_ref[my_var["py_z__ref%s"%k_ref]]
    my_var["py_next_z_pos__ref%s"%k_ref] = positions_ref[my_var["py_next_z__ref%s"%k_ref]]#STANDARD POINT
        
    #take three vector
    my_var["py_x_vector__ref%s"%k_ref] = my_var["py_x_pos__ref%s"%k_ref] - my_var["py_z_pos__ref%s"%k_ref]

    if(my_var["py_x_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["py_x_vector__ref%s"%k_ref][0] = my_var["py_x_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["py_x_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["py_x_vector__ref%s"%k_ref][0] = my_var["py_x_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["py_x_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["py_x_vector__ref%s"%k_ref][1] = my_var["py_x_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["py_x_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["py_x_vector__ref%s"%k_ref][1] = my_var["py_x_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["py_x_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["py_x_vector__ref%s"%k_ref][2] = my_var["py_x_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["py_x_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["py_x_vector__ref%s"%k_ref][2] = my_var["py_x_vector__ref%s"%k_ref][2] + pbc_z_ref

    my_var["py_nx_vector__ref%s"%k_ref] = my_var["py_nx_pos__ref%s"%k_ref] - my_var["py_z_pos__ref%s"%k_ref]

    if(my_var["py_nx_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["py_nx_vector__ref%s"%k_ref][0] = my_var["py_nx_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["py_nx_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["py_nx_vector__ref%s"%k_ref][0] = my_var["py_nx_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["py_nx_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["py_nx_vector__ref%s"%k_ref][1] = my_var["py_nx_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["py_nx_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["py_nx_vector__ref%s"%k_ref][1] = my_var["py_nx_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["py_nx_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["py_nx_vector__ref%s"%k_ref][2] = my_var["py_nx_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["py_nx_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["py_nx_vector__ref%s"%k_ref][2] = my_var["py_nx_vector__ref%s"%k_ref][2] + pbc_z_ref

    my_var["py_z_vector__ref%s"%k_ref] = my_var["py_next_z_pos__ref%s"%k_ref] - my_var["py_z_pos__ref%s"%k_ref]

    if(my_var["py_z_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["py_z_vector__ref%s"%k_ref][0] = my_var["py_z_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["py_z_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["py_z_vector__ref%s"%k_ref][0] = my_var["py_z_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["py_z_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["py_z_vector__ref%s"%k_ref][1] = my_var["py_z_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["py_z_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["py_z_vector__ref%s"%k_ref][1] = my_var["py_z_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["py_z_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["py_z_vector__ref%s"%k_ref][2] = my_var["py_z_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["py_z_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["py_z_vector__ref%s"%k_ref][2] = my_var["py_z_vector__ref%s"%k_ref][2] + pbc_z_ref

    #transfer vector to matrix form        
    my_var["py_x_matrix__ref%s"%k_ref] = np.array(my_var["py_x_vector__ref%s"%k_ref])
    my_var["py_nx_matrix__ref%s"%k_ref] = np.array(my_var["py_nx_vector__ref%s"%k_ref])
    my_var["py_z_matrix__ref%s"%k_ref] = np.array(my_var["py_z_vector__ref%s"%k_ref])
    
    #ref matrix
    my_var["py_matrix__ref%s"%k_ref] = np.column_stack([my_var["py_x_matrix__ref%s"%k_ref],my_var["py_nx_matrix__ref%s"%k_ref],my_var["py_z_matrix__ref%s"%k_ref]])
    #inverse ref matrix 
    my_var["py_matrix_inv_ref%s"%k_ref] = np.linalg.pinv(my_var["py_matrix__ref%s"%k_ref])

    #+z direciton three vector
    #take +X -X +Y NEXT +Y    
    my_var["tetragonal_id_%s_copy"%k_ref+"_pz_ref"] = my_var["tetragonal_id_%s"%k_ref+"_pz_ref"].copy()
        
    #take id +X -X +Y NEXT +Y 
    my_var["pz_x__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_pz_ref"][1]
    my_var["pz_nx__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_pz_ref"][3]
    my_var["pz_y__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_pz_ref"][2]
    my_var["pz_next_y__ref%s"%k_ref] = my_var["tetragonal_id_%s_copy"%k_ref+"_pz_ref"][6]
        
    #take position
    my_var["pz_x_pos__ref%s"%k_ref] = positions_ref[my_var["pz_x__ref%s"%k_ref]]
    my_var["pz_nx_pos__ref%s"%k_ref] = positions_ref[my_var["pz_nx__ref%s"%k_ref]]
    my_var["pz_y_pos__ref%s"%k_ref] = positions_ref[my_var["pz_y__ref%s"%k_ref]]
    my_var["pz_next_y_pos__ref%s"%k_ref] = positions_ref[my_var["pz_next_y__ref%s"%k_ref]]#STANDARD POINT
        
    #take three vector
    my_var["pz_x_vector__ref%s"%k_ref] = my_var["pz_x_pos__ref%s"%k_ref] - my_var["pz_y_pos__ref%s"%k_ref]

    if(my_var["pz_x_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["pz_x_vector__ref%s"%k_ref][0] = my_var["pz_x_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["pz_x_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["pz_x_vector__ref%s"%k_ref][0] = my_var["pz_x_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["pz_x_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["pz_x_vector__ref%s"%k_ref][1] = my_var["pz_x_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["pz_x_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["pz_x_vector__ref%s"%k_ref][1] = my_var["pz_x_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["pz_x_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["pz_x_vector__ref%s"%k_ref][2] = my_var["pz_x_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["pz_x_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["pz_x_vector__ref%s"%k_ref][2] = my_var["pz_x_vector__ref%s"%k_ref][2] + pbc_z_ref

    my_var["pz_nx_vector__ref%s"%k_ref] = my_var["pz_nx_pos__ref%s"%k_ref] - my_var["pz_y_pos__ref%s"%k_ref]

    if(my_var["pz_nx_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["pz_nx_vector__ref%s"%k_ref][0] = my_var["pz_nx_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["pz_nx_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["pz_nx_vector__ref%s"%k_ref][0] = my_var["pz_nx_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["pz_nx_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["pz_nx_vector__ref%s"%k_ref][1] = my_var["pz_nx_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["pz_nx_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["pz_nx_vector__ref%s"%k_ref][1] = my_var["pz_nx_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["pz_nx_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["pz_nx_vector__ref%s"%k_ref][2] = my_var["pz_nx_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["pz_nx_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["pz_nx_vector__ref%s"%k_ref][2] = my_var["pz_nx_vector__ref%s"%k_ref][2] + pbc_z_ref

    my_var["pz_y_vector__ref%s"%k_ref] = my_var["pz_next_y_pos__ref%s"%k_ref] - my_var["pz_y_pos__ref%s"%k_ref]

    if(my_var["pz_y_vector__ref%s"%k_ref][0] >= pbc_x_ref*0.9):
        my_var["pz_y_vector__ref%s"%k_ref][0] = my_var["pz_y_vector__ref%s"%k_ref][0] - pbc_x_ref
    elif(my_var["pz_y_vector__ref%s"%k_ref][0] <= pbc_x_ref*0.9*-1):
        my_var["pz_y_vector__ref%s"%k_ref][0] = my_var["pz_y_vector__ref%s"%k_ref][0] + pbc_x_ref
    if(my_var["pz_y_vector__ref%s"%k_ref][1] >= pbc_y_ref*0.9):
        my_var["pz_y_vector__ref%s"%k_ref][1] = my_var["pz_y_vector__ref%s"%k_ref][1] - pbc_y_ref
    elif(my_var["pz_y_vector__ref%s"%k_ref][1] <= pbc_y_ref*0.9*-1):
        my_var["pz_y_vector__ref%s"%k_ref][1] = my_var["pz_y_vector__ref%s"%k_ref][1] + pbc_y_ref
    if(my_var["pz_y_vector__ref%s"%k_ref][2] >= pbc_z_ref*0.9):
        my_var["pz_y_vector__ref%s"%k_ref][2] = my_var["pz_y_vector__ref%s"%k_ref][2] - pbc_z_ref
    elif(my_var["pz_y_vector__ref%s"%k_ref][2] <= pbc_z_ref*0.9*-1):
        my_var["pz_y_vector__ref%s"%k_ref][2] = my_var["pz_y_vector__ref%s"%k_ref][2] + pbc_z_ref

    #transfer vector to matrix form            
    my_var["pz_x_matrix__ref%s"%k_ref] = np.array(my_var["pz_x_vector__ref%s"%k_ref])
    my_var["pz_nx_matrix__ref%s"%k_ref] = np.array(my_var["pz_nx_vector__ref%s"%k_ref])
    my_var["pz_y_matrix__ref%s"%k_ref] = np.array(my_var["pz_y_vector__ref%s"%k_ref])

    #ref matrix
    my_var["pz_matrix__ref%s"%k_ref] = np.column_stack([my_var["pz_x_matrix__ref%s"%k_ref],my_var["pz_nx_matrix__ref%s"%k_ref],my_var["pz_y_matrix__ref%s"%k_ref]])
    
    #inverse ref matrix 
    my_var["pz_matrix_inv_ref%s"%k_ref] = np.linalg.pinv(my_var["pz_matrix__ref%s"%k_ref])


def compute_myproperty(frame, data):
    property_list = []
    data.cell_.pbc = (True, True, False)

#calculate ref
    pbc_x = float(data.cell_[0][0])
    pbc_y = float(data.cell_[1][1])
    pbc_z = float(data.cell_[2][2])

    # Initialize neighbor finder object.
    # Visit the 6 nearest neighbors of each particle.
    N = 6
    finder = NearestNeighborFinder(N, data)

    ptypes = data.particles['Particle Type']
    positions = data.particles['Position']

    # Loop over all input particles:
    for index in range(data.particles.count):
        #print(neigh.index, neigh.distance, neigh.delta)
        neighbors = [ (neigh.index, neigh.delta) for neigh in finder.find(index) ]
        #print neighbor indices  
        resorted_neighbors_x = sorted( neighbors , key=lambda k: [k[1][0], k[1][1], k[1][2]], reverse=True )
        resorted_neighbors_y = sorted( neighbors , key=lambda k: [k[1][1], k[1][0], k[1][2]], reverse=True )
        resorted_neighbors_z = sorted( neighbors , key=lambda k: [k[1][2], k[1][0], k[1][1]], reverse=True )
        #print rearranged neighbor indices                                                                                                                                                                                                            
        my_var["list_%s"%index+"_x"] = [n[0] for n in resorted_neighbors_x]
        my_var["list_%s"%index+"_y"] = [n[0] for n in resorted_neighbors_y]
        my_var["list_%s"%index+"_z"] = [n[0] for n in resorted_neighbors_z]
    
    for i in range(data.particles.count):
        
        my_var["make_remove_px_%s"%i] = []
        my_var["make_remove_py_%s"%i] = []
        my_var["make_remove_pz_%s"%i] = []
        my_var["make_remove_nx_%s"%i] = []
        my_var["make_remove_ny_%s"%i] = []
        my_var["make_remove_nz_%s"%i] = []
        
        my_var["make_list_%s"%i+"_x"] = my_var["list_%s"%i+"_x"].copy()
        my_var["make_remove_px_%s"%i].insert(0,my_var["make_list_%s"%i+"_x"].pop(0))
        my_var["make_remove_nx_%s"%i].insert(0,my_var["make_list_%s"%i+"_x"].pop(4))
        
        my_var["make_list_%s"%i+"_y"] = my_var["list_%s"%i+"_y"].copy()
        my_var["make_remove_py_%s"%i].insert(0,my_var["make_list_%s"%i+"_y"].pop(0))
        my_var["make_remove_ny_%s"%i].insert(0,my_var["make_list_%s"%i+"_y"].pop(4))
        
        my_var["make_list_%s"%i+"_z"] = my_var["list_%s"%i+"_z"].copy()
        my_var["make_remove_pz_%s"%i].insert(0,my_var["make_list_%s"%i+"_z"].pop(0))
        my_var["make_remove_nz_%s"%i].insert(0,my_var["make_list_%s"%i+"_z"].pop(4))
        
        my_var["neigh_list_%s"%i] = []
        my_var["neigh_list_%s"%i] = my_var["make_remove_px_%s"%i] + my_var["make_remove_py_%s"%i] + my_var["make_remove_pz_%s"%i] + my_var["make_remove_nx_%s"%i] + my_var["make_remove_ny_%s"%i] + my_var["make_remove_nz_%s"%i]
        my_var["neigh_list_%s"%i].insert(0,i)
        
    for j in range(data.particles.count):
        # +x direction
        my_var["list_%s"%j+"_px"] = my_var["neigh_list_%s"%j].copy()
        remove_px = my_var["list_%s"%j+"_px"].pop(1)
        remove_nx = my_var["list_%s"%j+"_px"].pop(3)
        #find next atom along +x direction
        my_var["list_%s"%j+"_npx"] = my_var["neigh_list_%s"%remove_px].copy()
        remove_cen_px = my_var["list_%s"%j+"_npx"].pop(0)
        remove_npx = my_var["list_%s"%j+"_npx"].pop(0)
        remove_nnx = my_var["list_%s"%j+"_npx"].pop(2)
        
        my_var["tetragonal_id_%s"%j+"_px"] = my_var["list_%s"%j+"_px"].copy() + my_var["list_%s"%j+"_npx"].copy()
        
        # ID +Y +Z -Y -Z +Y +Z -Y -Z

        #+Y DIRECTION
        my_var["list_%s"%j+"_py"] = my_var["neigh_list_%s"%j].copy()
        remove_py = my_var["list_%s"%j+"_py"].pop(2)
        remove_ny = my_var["list_%s"%j+"_py"].pop(4)
        #find next atom along +y direction
        my_var["list_%s"%j+"_npy"] = my_var["neigh_list_%s"%remove_py].copy()
        remove_cen_py = my_var["list_%s"%j+"_npy"].pop(0)
        remove_npy = my_var["list_%s"%j+"_npy"].pop(1)
        remove_nny = my_var["list_%s"%j+"_npy"].pop(3)
        
        my_var["tetragonal_id_%s"%j+"_py"] = my_var["list_%s"%j+"_py"].copy() + my_var["list_%s"%j+"_npy"].copy()
                
        #+z direction
        my_var["list_%s"%j+"_pz"] = my_var["neigh_list_%s"%j].copy()
        remove_pz = my_var["list_%s"%j+"_pz"].pop(3)
        remove_nz = my_var["list_%s"%j+"_pz"].pop(5)
        #find next atom along +z direction
        my_var["list_%s"%j+"_npz"] = my_var["neigh_list_%s"%remove_pz].copy()
        remove_cen_pz = my_var["list_%s"%j+"_npz"].pop(0)
        remove_npz = my_var["list_%s"%j+"_npz"].pop(2)
        remove_nnz = my_var["list_%s"%j+"_npz"].pop(4)
        my_var["tetragonal_id_%s"%j+"_pz"] = my_var["list_%s"%j+"_pz"].copy() + my_var["list_%s"%j+"_npz"].copy()

        
    for k in range(data.particles.count):
        
    #+x direciton three vector
        #take +Y -Y +Z NEXT +Z
        
        my_var["tetragonal_id_%s_copy"%k+"_px"] = my_var["tetragonal_id_%s"%k+"_px"].copy()
        
        #take id
        my_var["px_y_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_px"][1]
        my_var["px_ny_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_px"][3]
        my_var["px_z_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_px"][2]
        my_var["px_next_z_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_px"][6]
        
        #take position
        my_var["px_y_pos_%s"%k] = positions[my_var["px_y_%s"%k]]
        my_var["px_ny_pos_%s"%k] = positions[my_var["px_ny_%s"%k]]
        my_var["px_z_pos_%s"%k] = positions[my_var["px_z_%s"%k]]
        my_var["px_next_z_pos_%s"%k] = positions[my_var["px_next_z_%s"%k]]#STANDARD POINT
        
        #take three vector
        my_var["px_y_vector_%s"%k] = my_var["px_y_pos_%s"%k] - my_var["px_z_pos_%s"%k]

        if(my_var["px_y_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["px_y_vector_%s"%k][0] = my_var["px_y_vector_%s"%k][0] - pbc_x
        elif(my_var["px_y_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["px_y_vector_%s"%k][0] = my_var["px_y_vector_%s"%k][0] + pbc_x
        if(my_var["px_y_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["px_y_vector_%s"%k][1] = my_var["px_y_vector_%s"%k][1] - pbc_y
        elif(my_var["px_y_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["px_y_vector_%s"%k][1] = my_var["px_y_vector_%s"%k][1] + pbc_y
        if(my_var["px_y_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["px_y_vector_%s"%k][2] = my_var["px_y_vector_%s"%k][2] - pbc_z
        elif(my_var["px_y_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["px_y_vector_%s"%k][2] = my_var["px_y_vector_%s"%k][2] + pbc_z

        my_var["px_ny_vector_%s"%k] = my_var["px_ny_pos_%s"%k] - my_var["px_z_pos_%s"%k]

        if(my_var["px_ny_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["px_ny_vector_%s"%k][0] = my_var["px_ny_vector_%s"%k][0] - pbc_x
        elif(my_var["px_ny_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["px_ny_vector_%s"%k][0] = my_var["px_ny_vector_%s"%k][0] + pbc_x
        if(my_var["px_ny_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["px_ny_vector_%s"%k][1] = my_var["px_ny_vector_%s"%k][1] - pbc_y
        elif(my_var["px_ny_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["px_ny_vector_%s"%k][1] = my_var["px_ny_vector_%s"%k][1] + pbc_y
        if(my_var["px_ny_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["px_ny_vector_%s"%k][2] = my_var["px_ny_vector_%s"%k][2] - pbc_z
        elif(my_var["px_ny_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["px_ny_vector_%s"%k][2] = my_var["px_ny_vector_%s"%k][2] + pbc_z

        my_var["px_z_vector_%s"%k] = my_var["px_next_z_pos_%s"%k] - my_var["px_z_pos_%s"%k]

        if(my_var["px_z_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["px_z_vector_%s"%k][0] = my_var["px_z_vector_%s"%k][0] - pbc_x
        elif(my_var["px_z_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["px_z_vector_%s"%k][0] = my_var["px_z_vector_%s"%k][0] + pbc_x
        if(my_var["px_z_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["px_z_vector_%s"%k][1] = my_var["px_z_vector_%s"%k][1] - pbc_y
        elif(my_var["px_z_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["px_z_vector_%s"%k][1] = my_var["px_z_vector_%s"%k][1] + pbc_y
        if(my_var["px_z_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["px_z_vector_%s"%k][2] = my_var["px_z_vector_%s"%k][2] - pbc_z
        elif(my_var["px_z_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["px_z_vector_%s"%k][2] = my_var["px_z_vector_%s"%k][2] + pbc_z
        
        #TAKE VECTOR MATRIX
        my_var["px_y_matrix_%s"%k] = np.array(my_var["px_y_vector_%s"%k])
        my_var["px_ny_matrix_%s"%k] = np.array(my_var["px_ny_vector_%s"%k])
        my_var["px_z_matrix_%s"%k] = np.array(my_var["px_z_vector_%s"%k])
        my_var["px_matrix_%s"%k] = np.column_stack([my_var["px_ny_matrix_%s"%k],my_var["px_y_matrix_%s"%k],my_var["px_z_matrix_%s"%k]])


        #X


        #take F&& and Ft
        my_var["px_f_%s"%k] = np.dot(my_var["px_matrix_%s"%k],my_var["px_matrix_inv_ref%s"%k])
        my_var["px_ft_%s"%k] = np.transpose(my_var["px_f_%s"%k])
        
        #u ftf
        my_var["px_u2_%s"%k] = np.dot(my_var["px_ft_%s"%k],my_var["px_f_%s"%k])
        my_var["px_u_%s"%k] = sl.sqrtm(my_var["px_u2_%s"%k])

        #calculate point
        my_var["px_15_%s"%k] = np.zeros((1,15),dtype = np.float32)
        my_var["px_m1_score%s"%k] = 0
        my_var["px_m2_score%s"%k] = 0
        my_var["px_m3_score%s"%k] = 0
        my_var["px_m4_score%s"%k] = 0
        my_var["px_o1_score%s"%k] = 0
        my_var["px_o2_score%s"%k] = 0
        my_var["px_a_score%s"%k] = 0
        my_var["px_highest_score%s"%k] = 0
        my_var["px_variant%s"%k] = 0


        #A
        if my_var["px_u_%s"%k][0][0] - 1 >= tol:
            my_var["px_15_%s"%k][0,0] = 1
        #B
        if my_var["px_u_%s"%k][1][1] - 1 >= tol :
            my_var["px_15_%s"%k][0,1] = 1
        #C
        if my_var["px_u_%s"%k][2][2] - 1 >= tol :
            my_var["px_15_%s"%k][0,2] = 1
        #D
        if my_var["px_u_%s"%k][0][0] - 1 < ntol :
            my_var["px_15_%s"%k][0,3] = 1
        #E
        if my_var["px_u_%s"%k][1][1] - 1 < ntol :
            my_var["px_15_%s"%k][0,4] = 1
        #F
        if my_var["px_u_%s"%k][2][2] - 1 < ntol :
            my_var["px_15_%s"%k][0,5] = 1
        #G
        if my_var["px_u_%s"%k][0][1] >= tol :
            my_var["px_15_%s"%k][0,6] = 1
        #H
        if my_var["px_u_%s"%k][0][2] >= tol :
            my_var["px_15_%s"%k][0,7] = 1
        #I
        if my_var["px_u_%s"%k][1][2] >= tol :
            my_var["px_15_%s"%k][0,8] = 1
        #J
        if my_var["px_u_%s"%k][0][1] <= ntol :
            my_var["px_15_%s"%k][0,9] = 1
        #K
        if my_var["px_u_%s"%k][0][2] <= ntol :
            my_var["px_15_%s"%k][0,10] = 1
        #L
        if my_var["px_u_%s"%k][1][2] <= ntol :
            my_var["px_15_%s"%k][0,11] = 1
        #M
        if abs(abs(my_var["px_u_%s"%k][0][1]) - abs(my_var["px_u_%s"%k][0][2])) <= 2*tol:
            my_var["px_15_%s"%k][0,12] = 1
        #N
        if abs(abs(my_var["px_u_%s"%k][0][1]) - abs(my_var["px_u_%s"%k][1][2])) <= 2*tol:
            my_var["px_15_%s"%k][0,13] = 1
        #O
        if abs(abs(my_var["px_u_%s"%k][0][2]) - abs(my_var["px_u_%s"%k][1][2])) <= 2*tol:
            my_var["px_15_%s"%k][0,14] = 1

        #classify
        my_var["px_15_same_m1%s"%k] = my_var["px_15_%s"%k] == M1
        my_var["px_m1_score%s"%k] = np.sum(my_var["px_15_same_m1%s"%k])

        my_var["px_15_same_m2%s"%k] = my_var["px_15_%s"%k] == M2
        my_var["px_m2_score%s"%k] = np.sum(my_var["px_15_same_m2%s"%k])

        my_var["px_15_same_m3%s"%k] = my_var["px_15_%s"%k] == M3
        my_var["px_m3_score%s"%k] = np.sum(my_var["px_15_same_m3%s"%k])

        my_var["px_15_same_m4%s"%k] = my_var["px_15_%s"%k] == M4
        my_var["px_m4_score%s"%k] = np.sum(my_var["px_15_same_m4%s"%k])

        my_var["px_15_same_o1%s"%k] = my_var["px_15_%s"%k] == O1
        my_var["px_o1_score%s"%k] = np.sum(my_var["px_15_same_o1%s"%k])

        my_var["px_15_same_o2%s"%k] = my_var["px_15_%s"%k] == O2
        my_var["px_o2_score%s"%k] = np.sum(my_var["px_15_same_o2%s"%k])

        my_var["px_15_same_a%s"%k] = my_var["px_15_%s"%k] == A
        my_var["px_a_score%s"%k] = np.sum(my_var["px_15_same_a%s"%k])

        my_var["px_highest_score%s"%k] = max(my_var["px_m1_score%s"%k],my_var["px_m2_score%s"%k],my_var["px_m3_score%s"%k],my_var["px_m4_score%s"%k],my_var["px_o1_score%s"%k],my_var["px_o2_score%s"%k],my_var["px_a_score%s"%k])

        my_var["px_score_list%s"%k] = [my_var["px_m1_score%s"%k],my_var["px_m2_score%s"%k],my_var["px_m3_score%s"%k],my_var["px_m4_score%s"%k],my_var["px_o1_score%s"%k],my_var["px_o2_score%s"%k],my_var["px_a_score%s"%k]]
        #determine variants
        if my_var["px_highest_score%s"%k] == my_var["px_m1_score%s"%k]:
            my_var["px_variant%s"%k] = "M1"
        elif my_var["px_highest_score%s"%k] == my_var["px_m2_score%s"%k]:
            my_var["px_variant%s"%k] = "M2"
        elif my_var["px_highest_score%s"%k] == my_var["px_m3_score%s"%k]:
            my_var["px_variant%s"%k] = "M3"
        elif my_var["px_highest_score%s"%k] == my_var["px_m4_score%s"%k]:
            my_var["px_variant%s"%k] = "M4"
        elif my_var["px_highest_score%s"%k] == my_var["px_o1_score%s"%k]:
            my_var["px_variant%s"%k] = "O1"
        elif my_var["px_highest_score%s"%k] == my_var["px_o2_score%s"%k]:
            my_var["px_variant%s"%k] = "O2"
        elif my_var["px_highest_score%s"%k] == my_var["px_a_score%s"%k]:
            my_var["px_variant%s"%k] = "A"

        my_var["px_unknown%s"%k] = np.sum(my_var["px_score_list%s"%k] == my_var["px_highest_score%s"%k])

        if my_var["px_unknown%s"%k] != 1:
            my_var["px_variant%s"%k] = "unknown"

#+y direciton three vector
        #take +x -x +Z NEXT +Z
        
        my_var["tetragonal_id_%s_copy"%k+"_py"] = my_var["tetragonal_id_%s"%k+"_py"].copy()
        
        #take id
        my_var["py_x_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_py"][1]
        my_var["py_nx_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_py"][3]
        my_var["py_z_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_py"][2]
        my_var["py_next_z_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_py"][6]
        
        #take position
        my_var["py_x_pos_%s"%k] = positions[my_var["py_x_%s"%k]]
        my_var["py_nx_pos_%s"%k] = positions[my_var["py_nx_%s"%k]]
        my_var["py_z_pos_%s"%k] = positions[my_var["py_z_%s"%k]]
        my_var["py_next_z_pos_%s"%k] = positions[my_var["py_next_z_%s"%k]]#STANDARD POINT
        
        #take three vector
        my_var["py_x_vector_%s"%k] = my_var["py_x_pos_%s"%k] - my_var["py_z_pos_%s"%k]

        if(my_var["py_x_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["py_x_vector_%s"%k][0] = my_var["py_x_vector_%s"%k][0] - pbc_x
        elif(my_var["py_x_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["py_x_vector_%s"%k][0] = my_var["py_x_vector_%s"%k][0] + pbc_x
        if(my_var["py_x_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["py_x_vector_%s"%k][1] = my_var["py_x_vector_%s"%k][1] - pbc_y
        elif(my_var["py_x_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["py_x_vector_%s"%k][1] = my_var["py_x_vector_%s"%k][1] + pbc_y
        if(my_var["py_x_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["py_x_vector_%s"%k][2] = my_var["py_x_vector_%s"%k][2] - pbc_z
        elif(my_var["py_x_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["py_x_vector_%s"%k][2] = my_var["py_x_vector_%s"%k][2] + pbc_z

        my_var["py_nx_vector_%s"%k] = my_var["py_nx_pos_%s"%k] - my_var["py_z_pos_%s"%k]

        if(my_var["py_nx_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["py_nx_vector_%s"%k][0] = my_var["py_nx_vector_%s"%k][0] - pbc_x
        elif(my_var["py_nx_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["py_nx_vector_%s"%k][0] = my_var["py_nx_vector_%s"%k][0] + pbc_x
        if(my_var["py_nx_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["py_nx_vector_%s"%k][1] = my_var["py_nx_vector_%s"%k][1] - pbc_y
        elif(my_var["py_nx_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["py_nx_vector_%s"%k][1] = my_var["py_nx_vector_%s"%k][1] + pbc_y
        if(my_var["py_nx_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["py_nx_vector_%s"%k][2] = my_var["py_nx_vector_%s"%k][2] - pbc_z
        elif(my_var["py_nx_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["py_nx_vector_%s"%k][2] = my_var["py_nx_vector_%s"%k][2] + pbc_z

        my_var["py_z_vector_%s"%k] = my_var["py_next_z_pos_%s"%k] - my_var["py_z_pos_%s"%k]

        if(my_var["py_z_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["py_z_vector_%s"%k][0] = my_var["py_z_vector_%s"%k][0] - pbc_x
        elif(my_var["py_z_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["py_z_vector_%s"%k][0] = my_var["py_z_vector_%s"%k][0] + pbc_x
        if(my_var["py_z_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["py_z_vector_%s"%k][1] = my_var["py_z_vector_%s"%k][1] - pbc_y
        elif(my_var["py_z_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["py_z_vector_%s"%k][1] = my_var["py_z_vector_%s"%k][1] + pbc_y
        if(my_var["py_z_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["py_z_vector_%s"%k][2] = my_var["py_z_vector_%s"%k][2] - pbc_z
        elif(my_var["py_z_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["py_z_vector_%s"%k][2] = my_var["py_z_vector_%s"%k][2] + pbc_z
        
        #TAKE VECTOR MATRIX
        my_var["py_x_matrix_%s"%k] = np.array(my_var["py_x_vector_%s"%k])
        my_var["py_nx_matrix_%s"%k] = np.array(my_var["py_nx_vector_%s"%k])
        my_var["py_z_matrix_%s"%k] = np.array(my_var["py_z_vector_%s"%k])
        my_var["py_matrix_%s"%k] = np.column_stack([my_var["py_x_matrix_%s"%k],my_var["py_nx_matrix_%s"%k],my_var["py_z_matrix_%s"%k]])


        #Y


        #take F&& and Ft
        my_var["py_f_%s"%k] = np.dot(my_var["py_matrix_%s"%k],my_var["py_matrix_inv_ref%s"%k])
        my_var["py_ft_%s"%k] = np.transpose(my_var["py_f_%s"%k])
        
        #u ftf
        my_var["py_u2_%s"%k] = np.dot(my_var["py_ft_%s"%k],my_var["py_f_%s"%k])
        my_var["py_u_%s"%k] = sl.sqrtm(my_var["py_u2_%s"%k])

        #calculate point
        my_var["py_15_%s"%k] = np.zeros((1,15),dtype = np.float32)
        my_var["py_m5_score%s"%k] = 0
        my_var["py_m6_score%s"%k] = 0
        my_var["py_m7_score%s"%k] = 0
        my_var["py_m8_score%s"%k] = 0
        my_var["py_o3_score%s"%k] = 0
        my_var["py_o4_score%s"%k] = 0
        my_var["py_a_score%s"%k] = 0
        my_var["py_highest_score%s"%k] = 0
        my_var["py_variant%s"%k] = 0

        #A
        if my_var["py_u_%s"%k][0][0] - 1 >= tol :
            my_var["py_15_%s"%k][0,0] = 1
        #B
        if my_var["py_u_%s"%k][1][1] - 1 >= tol :
            my_var["py_15_%s"%k][0,1] = 1
        #C
        if my_var["py_u_%s"%k][2][2] - 1 >= tol :
            my_var["py_15_%s"%k][0,2] = 1
        #D
        if my_var["py_u_%s"%k][0][0] - 1 < ntol :
            my_var["py_15_%s"%k][0,3] = 1
        #E
        if my_var["py_u_%s"%k][1][1] - 1 < ntol :
            my_var["py_15_%s"%k][0,4] = 1
        #F
        if my_var["py_u_%s"%k][2][2] - 1 < ntol :
            my_var["py_15_%s"%k][0,5] = 1
        #G
        if my_var["py_u_%s"%k][0][1] >= tol :
            my_var["py_15_%s"%k][0,6] = 1
        #H
        if my_var["py_u_%s"%k][0][2] >= tol :
            my_var["py_15_%s"%k][0,7] = 1
        #I
        if my_var["py_u_%s"%k][1][2] >= tol :
            my_var["py_15_%s"%k][0,8] = 1
        #J
        if my_var["py_u_%s"%k][0][1] <= ntol :
            my_var["py_15_%s"%k][0,9] = 1
        #K
        if my_var["py_u_%s"%k][0][2] <= ntol :
            my_var["py_15_%s"%k][0,10] = 1
        #L
        if my_var["py_u_%s"%k][1][2] <= ntol :
            my_var["py_15_%s"%k][0,11] = 1
        #M
        if abs(abs(my_var["py_u_%s"%k][0][1]) - abs(my_var["py_u_%s"%k][0][2])) <= 2*tol:
            my_var["py_15_%s"%k][0,12] = 1
        #N
        if abs(abs(my_var["py_u_%s"%k][0][1]) - abs(my_var["py_u_%s"%k][1][2])) <= 2*tol:
            my_var["py_15_%s"%k][0,13] = 1
        #O
        if abs(abs(my_var["py_u_%s"%k][0][2]) - abs(my_var["py_u_%s"%k][1][2])) <= 2*tol:
            my_var["py_15_%s"%k][0,14] = 1

        #classify
        my_var["py_15_same_m5%s"%k] = my_var["py_15_%s"%k] == M5
        my_var["py_m5_score%s"%k] = np.sum(my_var["py_15_same_m5%s"%k])

        my_var["py_15_same_m6%s"%k] = my_var["py_15_%s"%k] == M6
        my_var["py_m6_score%s"%k] = np.sum(my_var["py_15_same_m6%s"%k])

        my_var["py_15_same_m7%s"%k] = my_var["py_15_%s"%k] == M7
        my_var["py_m7_score%s"%k] = np.sum(my_var["py_15_same_m7%s"%k])

        my_var["py_15_same_m8%s"%k] = my_var["py_15_%s"%k] == M8
        my_var["py_m8_score%s"%k] = np.sum(my_var["py_15_same_m8%s"%k])

        my_var["py_15_same_o3%s"%k] = my_var["py_15_%s"%k] == O3
        my_var["py_o3_score%s"%k] = np.sum(my_var["py_15_same_o3%s"%k])

        my_var["py_15_same_o4%s"%k] = my_var["py_15_%s"%k] == O4
        my_var["py_o4_score%s"%k] = np.sum(my_var["py_15_same_o4%s"%k])

        my_var["py_15_same_a%s"%k] = my_var["py_15_%s"%k] == A
        my_var["py_a_score%s"%k] = np.sum(my_var["py_15_same_a%s"%k])

        my_var["py_highest_score%s"%k] = max(my_var["py_m5_score%s"%k],my_var["py_m6_score%s"%k],my_var["py_m7_score%s"%k],my_var["py_m8_score%s"%k],my_var["py_o3_score%s"%k],my_var["py_o4_score%s"%k],my_var["py_a_score%s"%k])

        my_var["py_score_list%s"%k] = [my_var["py_m5_score%s"%k],my_var["py_m6_score%s"%k],my_var["py_m7_score%s"%k],my_var["py_m8_score%s"%k],my_var["py_o3_score%s"%k],my_var["py_o4_score%s"%k],my_var["py_a_score%s"%k]]
        #determine variants
        if my_var["py_highest_score%s"%k] == my_var["py_m5_score%s"%k]:
            my_var["py_variant%s"%k] = "M5"
        elif my_var["py_highest_score%s"%k] == my_var["py_m6_score%s"%k]:
            my_var["py_variant%s"%k] = "M6"
        elif my_var["py_highest_score%s"%k] == my_var["py_m7_score%s"%k]:
            my_var["py_variant%s"%k] = "M7"
        elif my_var["py_highest_score%s"%k] == my_var["py_m8_score%s"%k]:
            my_var["py_variant%s"%k] = "M8"
        elif my_var["py_highest_score%s"%k] == my_var["py_o3_score%s"%k]:
            my_var["py_variant%s"%k] = "O3"
        elif my_var["py_highest_score%s"%k] == my_var["py_o4_score%s"%k]:
            my_var["py_variant%s"%k] = "O4"
        elif my_var["py_highest_score%s"%k] == my_var["py_a_score%s"%k]:
            my_var["py_variant%s"%k] = "A"

        my_var["py_unknown%s"%k] = np.sum(my_var["py_score_list%s"%k] == my_var["py_highest_score%s"%k])
        
        if my_var["py_unknown%s"%k] != 1:
            my_var["py_variant%s"%k] = "unknown"

#+z direciton three vector
        #take +X -X +Y NEXT +Y
        
        my_var["tetragonal_id_%s_copy"%k+"_pz"] = my_var["tetragonal_id_%s"%k+"_pz"].copy()
        
        #take id
        my_var["pz_x_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_pz"][1]
        my_var["pz_nx_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_pz"][3]
        my_var["pz_y_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_pz"][2]
        my_var["pz_next_y_%s"%k] = my_var["tetragonal_id_%s_copy"%k+"_pz"][6]
        
        #take position
        my_var["pz_x_pos_%s"%k] = positions[my_var["pz_x_%s"%k]]
        my_var["pz_nx_pos_%s"%k] = positions[my_var["pz_nx_%s"%k]]
        my_var["pz_y_pos_%s"%k] = positions[my_var["pz_y_%s"%k]]
        my_var["pz_next_y_pos_%s"%k] = positions[my_var["pz_next_y_%s"%k]]#STANDARD POINT
        
        #take three vector
        my_var["pz_x_vector_%s"%k] = my_var["pz_x_pos_%s"%k] - my_var["pz_y_pos_%s"%k]

        if(my_var["pz_x_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["pz_x_vector_%s"%k][0] = my_var["pz_x_vector_%s"%k][0] - pbc_x
        elif(my_var["pz_x_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["pz_x_vector_%s"%k][0] = my_var["pz_x_vector_%s"%k][0] + pbc_x
        if(my_var["pz_x_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["pz_x_vector_%s"%k][1] = my_var["pz_x_vector_%s"%k][1] - pbc_y
        elif(my_var["pz_x_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["pz_x_vector_%s"%k][1] = my_var["pz_x_vector_%s"%k][1] + pbc_y
        if(my_var["pz_x_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["pz_x_vector_%s"%k][2] = my_var["pz_x_vector_%s"%k][2] - pbc_z
        elif(my_var["pz_x_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["pz_x_vector_%s"%k][2] = my_var["pz_x_vector_%s"%k][2] + pbc_z

        my_var["pz_nx_vector_%s"%k] = my_var["pz_nx_pos_%s"%k] - my_var["pz_y_pos_%s"%k]

        if(my_var["pz_nx_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["pz_nx_vector_%s"%k][0] = my_var["pz_nx_vector_%s"%k][0] - pbc_x
        elif(my_var["pz_nx_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["pz_nx_vector_%s"%k][0] = my_var["pz_nx_vector_%s"%k][0] + pbc_x
        if(my_var["pz_nx_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["pz_nx_vector_%s"%k][1] = my_var["pz_nx_vector_%s"%k][1] - pbc_y
        elif(my_var["pz_nx_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["pz_nx_vector_%s"%k][1] = my_var["pz_nx_vector_%s"%k][1] + pbc_y
        if(my_var["pz_nx_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["pz_nx_vector_%s"%k][2] = my_var["pz_nx_vector_%s"%k][2] - pbc_z
        elif(my_var["pz_nx_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["pz_nx_vector_%s"%k][2] = my_var["pz_nx_vector_%s"%k][2] + pbc_z

        my_var["pz_y_vector_%s"%k] = my_var["pz_next_y_pos_%s"%k] - my_var["pz_y_pos_%s"%k]

        if(my_var["pz_y_vector_%s"%k][0] >= pbc_x*0.9):
            my_var["pz_y_vector_%s"%k][0] = my_var["pz_y_vector_%s"%k][0] - pbc_x
        elif(my_var["pz_y_vector_%s"%k][0] <= pbc_x*0.9*-1):
            my_var["pz_y_vector_%s"%k][0] = my_var["pz_y_vector_%s"%k][0] + pbc_x
        if(my_var["pz_y_vector_%s"%k][1] >= pbc_y*0.9):
            my_var["pz_y_vector_%s"%k][1] = my_var["pz_y_vector_%s"%k][1] - pbc_y
        elif(my_var["pz_y_vector_%s"%k][1] <= pbc_y*0.9*-1):
            my_var["pz_y_vector_%s"%k][1] = my_var["pz_y_vector_%s"%k][1] + pbc_y
        if(my_var["pz_y_vector_%s"%k][2] >= pbc_z*0.9):
            my_var["pz_y_vector_%s"%k][2] = my_var["pz_y_vector_%s"%k][2] - pbc_z
        elif(my_var["pz_y_vector_%s"%k][2] <= pbc_z*0.9*-1):
            my_var["pz_y_vector_%s"%k][2] = my_var["pz_y_vector_%s"%k][2] + pbc_z
        
        #TAKE VECTOR MATRIX
        my_var["pz_x_matrix_%s"%k] = np.array(my_var["pz_x_vector_%s"%k])
        my_var["pz_nx_matrix_%s"%k] = np.array(my_var["pz_nx_vector_%s"%k])
        my_var["pz_y_matrix_%s"%k] = np.array(my_var["pz_y_vector_%s"%k])
        my_var["pz_matrix_%s"%k] = np.column_stack([my_var["pz_x_matrix_%s"%k],my_var["pz_nx_matrix_%s"%k],my_var["pz_y_matrix_%s"%k]])

        #Z


        
        #take F&& and Ft
        my_var["pz_f_%s"%k] = np.dot(my_var["pz_matrix_%s"%k],my_var["pz_matrix_inv_ref%s"%k])
        my_var["pz_ft_%s"%k] = np.transpose(my_var["pz_f_%s"%k])
        
        #u ftf
        my_var["pz_u2_%s"%k] = np.dot(my_var["pz_ft_%s"%k],my_var["pz_f_%s"%k])
        my_var["pz_u_%s"%k] = sl.sqrtm(my_var["pz_u2_%s"%k])

        #calculate point
        my_var["pz_15_%s"%k] = np.zeros((1,15),dtype = np.float32)
        my_var["pz_m9_score%s"%k] = 0
        my_var["pz_m10_score%s"%k] = 0
        my_var["pz_m11_score%s"%k] = 0
        my_var["pz_m12_score%s"%k] = 0
        my_var["pz_o5_score%s"%k] = 0
        my_var["pz_o6_score%s"%k] = 0
        my_var["pz_a_score%s"%k] = 0
        my_var["pz_highest_score%s"%k] = 0
        my_var["pz_variant%s"%k] = 0

        #A
        if my_var["pz_u_%s"%k][0][0] - 1 >= tol :
            my_var["pz_15_%s"%k][0,0] = 1
        #B
        if my_var["pz_u_%s"%k][1][1] - 1 >= tol :
            my_var["pz_15_%s"%k][0,1] = 1
        #C
        if my_var["pz_u_%s"%k][2][2] - 1 >= tol :
            my_var["pz_15_%s"%k][0,2] = 1
        #D
        if my_var["pz_u_%s"%k][0][0] - 1 < ntol :
            my_var["pz_15_%s"%k][0,3] = 1
        #E
        if my_var["pz_u_%s"%k][1][1] - 1 < ntol :
            my_var["pz_15_%s"%k][0,4] = 1
        #F
        if my_var["pz_u_%s"%k][2][2] - 1 < ntol :
            my_var["pz_15_%s"%k][0,5] = 1
        #G
        if my_var["pz_u_%s"%k][0][1] >= tol :
            my_var["pz_15_%s"%k][0,6] = 1
        #H
        if my_var["pz_u_%s"%k][0][2] >= tol :
            my_var["pz_15_%s"%k][0,7] = 1
        #I
        if my_var["pz_u_%s"%k][1][2] >= tol :
            my_var["pz_15_%s"%k][0,8] = 1
        #J
        if my_var["pz_u_%s"%k][0][1] <= ntol :
            my_var["pz_15_%s"%k][0,9] = 1
        #K
        if my_var["pz_u_%s"%k][0][2] <= ntol :
            my_var["pz_15_%s"%k][0,10] = 1
        #L
        if my_var["pz_u_%s"%k][1][2] <= ntol :
            my_var["pz_15_%s"%k][0,11] = 1
        #M
        if abs(abs(my_var["pz_u_%s"%k][0][1]) - abs(my_var["pz_u_%s"%k][0][2])) <= 2*tol:
            my_var["pz_15_%s"%k][0,12] = 1
        #N
        if abs(abs(my_var["pz_u_%s"%k][0][1]) - abs(my_var["pz_u_%s"%k][1][2])) <= 2*tol:
            my_var["pz_15_%s"%k][0,13] = 1
        #O
        if abs(abs(my_var["pz_u_%s"%k][0][2]) - abs(my_var["pz_u_%s"%k][1][2])) <= 2*tol:
            my_var["pz_15_%s"%k][0,14] = 1

        #classify
        my_var["pz_15_same_m9%s"%k] = my_var["pz_15_%s"%k] == M9
        my_var["pz_m9_score%s"%k] = np.sum(my_var["pz_15_same_m9%s"%k])

        my_var["pz_15_same_m10%s"%k] = my_var["pz_15_%s"%k] == M10
        my_var["pz_m10_score%s"%k] = np.sum(my_var["pz_15_same_m10%s"%k])

        my_var["pz_15_same_m11%s"%k] = my_var["pz_15_%s"%k] == M11
        my_var["pz_m11_score%s"%k] = np.sum(my_var["pz_15_same_m11%s"%k])

        my_var["pz_15_same_m12%s"%k] = my_var["pz_15_%s"%k] == M12
        my_var["pz_m12_score%s"%k] = np.sum(my_var["pz_15_same_m12%s"%k])

        my_var["pz_15_same_o5%s"%k] = my_var["pz_15_%s"%k] == O5
        my_var["pz_o5_score%s"%k] = np.sum(my_var["pz_15_same_o5%s"%k])

        my_var["pz_15_same_o6%s"%k] = my_var["pz_15_%s"%k] == O6
        my_var["pz_o6_score%s"%k] = np.sum(my_var["pz_15_same_o6%s"%k])

        my_var["pz_15_same_a%s"%k] = my_var["pz_15_%s"%k] == A
        my_var["pz_a_score%s"%k] = np.sum(my_var["pz_15_same_a%s"%k])

        my_var["pz_highest_score%s"%k] = max(my_var["pz_m9_score%s"%k],my_var["pz_m10_score%s"%k],my_var["pz_m11_score%s"%k],my_var["pz_m12_score%s"%k],my_var["pz_o5_score%s"%k],my_var["pz_o6_score%s"%k],my_var["pz_a_score%s"%k])

        my_var["pz_score_list%s"%k] = [my_var["pz_m9_score%s"%k],my_var["pz_m10_score%s"%k],my_var["pz_m11_score%s"%k],my_var["pz_m12_score%s"%k],my_var["pz_o5_score%s"%k],my_var["pz_o6_score%s"%k],my_var["pz_a_score%s"%k]]
        #determine variants
        if my_var["pz_highest_score%s"%k] == my_var["pz_m9_score%s"%k]:
            my_var["pz_variant%s"%k] = "M9"
        elif my_var["pz_highest_score%s"%k] == my_var["pz_m10_score%s"%k]:
            my_var["pz_variant%s"%k] = "M10"
        elif my_var["pz_highest_score%s"%k] == my_var["pz_m11_score%s"%k]:
            my_var["pz_variant%s"%k] = "M11"
        elif my_var["pz_highest_score%s"%k] == my_var["pz_m12_score%s"%k]:
            my_var["pz_variant%s"%k] = "M12"
        elif my_var["pz_highest_score%s"%k] == my_var["pz_o5_score%s"%k]:
            my_var["pz_variant%s"%k] = "O5"
        elif my_var["pz_highest_score%s"%k] == my_var["pz_o6_score%s"%k]:
            my_var["pz_variant%s"%k] = "O6"
        elif my_var["pz_highest_score%s"%k] == my_var["pz_a_score%s"%k]:
            my_var["pz_variant%s"%k] = "A"

        my_var["pz_unknown%s"%k] = np.sum(my_var["pz_score_list%s"%k] == my_var["pz_highest_score%s"%k])
        
        if my_var["pz_unknown%s"%k] != 1:
            my_var["pz_variant%s"%k] = "unknown"



        #Determine final variant
        my_var["variant%s"%k] = 0
        my_var["variant_color%s"%k] = 0
        my_var["highest_score_list%s"%k] = [my_var["px_highest_score%s"%k],my_var["py_highest_score%s"%k],my_var["pz_highest_score%s"%k]]
        my_var["highest_score%s"%k] = max(my_var["highest_score_list%s"%k])
        my_var["unknown%s"%k] = np.sum(my_var["highest_score_list%s"%k] == my_var["highest_score%s"%k])
        


        if my_var["highest_score%s"%k] == my_var["pz_highest_score%s"%k]:
            my_var["variant%s"%k] = my_var["pz_variant%s"%k]

        elif my_var["highest_score%s"%k] == my_var["px_highest_score%s"%k]:
            my_var["variant%s"%k] = my_var["px_variant%s"%k]
        elif my_var["highest_score%s"%k] == my_var["py_highest_score%s"%k]:
            my_var["variant%s"%k] = my_var["py_variant%s"%k]

        
        if my_var["variant%s"%k] != "A" and my_var["unknown%s"%k] != 1:
            my_var["variant%s"%k] = "unknown"


        

        #color
        if my_var["variant%s"%k] == "M1":
            my_var["variant_color%s"%k] = 0
        elif my_var["variant%s"%k] == "M2":
            my_var["variant_color%s"%k] = 1
        elif my_var["variant%s"%k] == "M3":
            my_var["variant_color%s"%k] = 2
        elif my_var["variant%s"%k] == "M4":
            my_var["variant_color%s"%k] = 3
        elif my_var["variant%s"%k] == "M5":
            my_var["variant_color%s"%k] = 4
        elif my_var["variant%s"%k] == "M6":
            my_var["variant_color%s"%k] = 5
        elif my_var["variant%s"%k] == "M7":
            my_var["variant_color%s"%k] = 6
        elif my_var["variant%s"%k] == "M8":
            my_var["variant_color%s"%k] = 7
        elif my_var["variant%s"%k] == "M9":
            my_var["variant_color%s"%k] = 8
        elif my_var["variant%s"%k] == "M10":
            my_var["variant_color%s"%k] = 9
        elif my_var["variant%s"%k] == "M11":
            my_var["variant_color%s"%k] = 10
        elif my_var["variant%s"%k] == "M12":
            my_var["variant_color%s"%k] = 11
        elif my_var["variant%s"%k] == "O1":
            my_var["variant_color%s"%k] = 12
        elif my_var["variant%s"%k] == "O2":
            my_var["variant_color%s"%k] = 13
        elif my_var["variant%s"%k] == "O3":
            my_var["variant_color%s"%k] = 14
        elif my_var["variant%s"%k] == "O4":
            my_var["variant_color%s"%k] = 15
        elif my_var["variant%s"%k] == "O5":
            my_var["variant_color%s"%k] = 16
        elif my_var["variant%s"%k] == "O6":
            my_var["variant_color%s"%k] = 17
        elif my_var["variant%s"%k] == "A":
            my_var["variant_color%s"%k] = 18
        elif my_var["variant%s"%k] == "unknown":
            my_var["variant_color%s"%k] = 19

        property_list.append(my_var["variant_color%s"%k])

    data.particles_.create_property('Variant', data=property_list)

node.modifiers.append(PythonScriptModifier(function = compute_myproperty))

for frame_index in range(node.source.num_frames):    
    data = node.compute(frame_index)
    f = open('E:/MD_PostProcess/data/10nm_tensile/result/post.%s'%frame_index,'w')
    export_file(node, 'E:/MD_PostProcess/data/10nm_tensile/result/post.%s'%frame_index, 'imd',columns = ['Particle Identifier','Particle Type','Position','Variant'],frame = frame_index)
    f.close()


