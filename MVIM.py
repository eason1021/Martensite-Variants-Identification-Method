from ovito.io import import_file
from ovito.io import export_file
from ovito.modifiers import *
from ovito.data import *
from ovito.data import NearestNeighborFinder
import numpy as np
import sys
import scipy.linalg as sl
import math
import datetime
#FIXME Adjust the scale
starttime = datetime.datetime.now()

#File path
file_path = '/data2/average.*'
ref_path = '/data2/large_model.dump'

#parameter area
#tolerance
tol = 0.006
ntol = tol*-1
my_var = {}
switch = False
#periodic boundary setting
pbc_setting = (True, True, False)
#rotation angle
rotation_angle = 90
# 0 = x-axis, 1 = y-axis, 2 = z-axis 
rotation_axis = 0

#input_data
node = import_file(file_path)
node_1 = import_file(file_path)
node_1.modifiers.append(SelectParticleTypeModifier(property='Particle Type', types={2}))
node_1.modifiers.append(DeleteSelectedParticlesModifier())
node_2 = import_file(file_path)
node_2.modifiers.append(SelectParticleTypeModifier(property='Particle Type', types={1}))
node_2.modifiers.append(DeleteSelectedParticlesModifier())

node_ref = import_file(ref_path)
node_ref_1 = import_file(ref_path)
node_ref_1.modifiers.append(SelectParticleTypeModifier(property='Particle Type', types={2}))
node_ref_1.modifiers.append(DeleteSelectedParticlesModifier())
node_ref_2 = import_file(ref_path)
node_ref_2.modifiers.append(SelectParticleTypeModifier(property='Particle Type', types={1}))
node_ref_2.modifiers.append(DeleteSelectedParticlesModifier())

def rotate(theta, axis, data, pbc_setting, bool):
	if bool == True :
		rotate_data = data.clone()
		theta = np.deg2rad(theta)  # time-dependent angle of rotation
		if axis == 0:
			tm = [[1, 0, 0, 0],
				[0 , np.cos(theta), -np.sin(theta), 0],
				[0, np.sin(theta), np.cos(theta), 0]]
		elif axis == 1:
			tm = [[np.cos(theta), 0, np.sin(theta), 0],
				[0,  1, 0, 0],
				[-np.sin(theta), 0, np.cos(theta), 0]]
		elif axis == 2:
			tm = [[np.cos(theta), -np.sin(theta), 0, 0],
				[np.sin(theta),  np.cos(theta), 0, 0],
				[            0,              0, 1, 0]]
		# Execute AffineTransformationModifier as a sub-operation:
		rotate_data.apply(AffineTransformationModifier(transformation = tm))
		if rotation_axis == 0 :
			if theta > 0 :
				rotate_data.cell_[1][1] = data.cell_[2][2] * -1
				rotate_data.cell_[2][2] = data.cell_[1][1]
			elif theta < 0 :
				rotate_data.cell_[1][1] = data.cell_[2][2]
				rotate_data.cell_[2][2] = data.cell_[1][1] * -1
			rotate_data.cell_[0][1] = 0
			rotate_data.cell_[0][2] = 0
			rotate_data.cell_[1][0] = 0
			rotate_data.cell_[1][2] = 0
			rotate_data.cell_[2][0] = 0
			rotate_data.cell_[2][1] = 0
			rotate_data.cell_.pbc = pbc_setting
		elif rotation_axis == 1 :
			if theta > 0:
				rotate_data.cell_[0][0] = data.cell_[2][2]
				rotate_data.cell_[2][2] = data.cell_[0][0] * -1
			elif theta < 0:
				rotate_data.cell_[0][0] = data.cell_[2][2] * -1
				rotate_data.cell_[2][2] = data.cell_[0][0]
			rotate_data.cell_[0][1] = 0
			rotate_data.cell_[0][2] = 0
			rotate_data.cell_[1][0] = 0
			rotate_data.cell_[1][2] = 0
			rotate_data.cell_[2][0] = 0
			rotate_data.cell_[2][1] = 0
			rotate_data.cell_.pbc = pbc_setting
		elif rotation_axis == 2 :
			if theta > 0:
				rotate_data.cell_[0][0] = data.cell_[1][1] * -1
				rotate_data.cell_[1][1] = data.cell_[0][0] 
			elif theta < 0:
				rotate_data.cell_[0][0] = data.cell_[1][1] 
				rotate_data.cell_[1][1] = data.cell_[0][0] * -1
			rotate_data.cell_[0][1] = 0
			rotate_data.cell_[0][2] = 0
			rotate_data.cell_[1][0] = 0
			rotate_data.cell_[1][2] = 0
			rotate_data.cell_[2][0] = 0
			rotate_data.cell_[2][1] = 0
			rotate_data.cell_.pbc = pbc_setting

	else :
		rotate_data = data.clone()
	return rotate_data

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

#R direction
R1 = np.array([0,0,0,0,0,0,0,0,0,1,1,1,1,1,1])
R2 = np.array([0,0,0,0,0,0,1,0,1,0,1,0,1,1,1])
R3 = np.array([0,0,0,0,0,0,0,1,1,1,0,0,1,1,1])
R4 = np.array([0,0,0,0,0,0,1,1,0,0,0,1,1,1,1])

def compute_ref_vector(data, theta, axis, pbc_setting, switch):
	ref_list = []
	N = 6
	rotate_data = rotate(theta, axis, data, pbc_setting, switch)
	finder = NearestNeighborFinder(N, rotate_data)
	ptype = data.particles['Particle Type']
	position = data.particles['Position']

	#neighbors[n]第幾近
	#neighbors[n][0]index
	#neighbors[n][1]座標
	#Loop over all input particles:
	for index in range(data.particles.count):
		neighbors = [ (neigh.index, neigh.delta) for neigh in finder.find(index) ]
		neigh_list = [0]*6
		resorted_neighbors_x_ref = sorted( neighbors , key=lambda k: [k[1][0], k[1][1], k[1][2]], reverse=True )
		resorted_neighbors_y_ref = sorted( neighbors , key=lambda k: [k[1][1], k[1][0], k[1][2]], reverse=True )
		resorted_neighbors_z_ref = sorted( neighbors , key=lambda k: [k[1][2], k[1][0], k[1][1]], reverse=True )
		resorted_neighbors_nx_ref = sorted( neighbors , key=lambda k: [k[1][0], k[1][1], k[1][2]], reverse=False )
		resorted_neighbors_ny_ref = sorted( neighbors , key=lambda k: [k[1][1], k[1][0], k[1][2]], reverse=False )
		resorted_neighbors_nz_ref = sorted( neighbors , key=lambda k: [k[1][2], k[1][0], k[1][1]], reverse=False )
		neigh_list[0] = resorted_neighbors_x_ref[0]
		neigh_list[1] = resorted_neighbors_y_ref[0]
		neigh_list[2] = resorted_neighbors_z_ref[0]
		neigh_list[3] = resorted_neighbors_nx_ref[0]
		neigh_list[4] = resorted_neighbors_ny_ref[0]
		neigh_list[5] = resorted_neighbors_nz_ref[0]
		my_var["neigh_list%s"%index] = neigh_list

	for neigh_sort in range(data.particles.count):
		phase_index_list = [0]*8
		#+X
		phase_index_list[0] = my_var["neigh_list%s"%neigh_sort][0]
		#+Y
		phase_index_list[1] = my_var["neigh_list%s"%neigh_sort][1]
		#+Z
		phase_index_list[2] = my_var["neigh_list%s"%neigh_sort][2]
		#-X
		phase_index_list[3] = my_var["neigh_list%s"%neigh_sort][3]
		#-Y
		phase_index_list[4] = my_var["neigh_list%s"%neigh_sort][4]
		#-Z
		phase_index_list[5] = my_var["neigh_list%s"%neigh_sort][5]
		#+X+Z
		phase_index_list[6] = (my_var["neigh_list%s"%my_var['neigh_list%s'%neigh_sort][0][0]][2][0],tuple(((np.array(my_var['neigh_list%s'%my_var['neigh_list%s'%neigh_sort][0][0]][2][1]) + np.array((my_var['neigh_list%s'%neigh_sort][0][1]))).tolist())))
		#+Y+Z
		phase_index_list[7] = (my_var["neigh_list%s"%my_var['neigh_list%s'%neigh_sort][1][0]][2][0],tuple(((np.array(my_var['neigh_list%s'%my_var['neigh_list%s'%neigh_sort][1][0]][2][1]) + np.array((my_var['neigh_list%s'%neigh_sort][1][1]))).tolist())))
		my_var["phase_index_list%s"%neigh_sort] = phase_index_list

	for transform_matrix in range(data.particles.count):
		#X Direction 
		phase_vector_x = [0]*3
		#+Y - +Z
		phase_vector_x[0] = np.array(my_var["phase_index_list%s"%transform_matrix][1][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#-Y - +Z
		phase_vector_x[1] = np.array(my_var["phase_index_list%s"%transform_matrix][4][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#+X+Z - +Z
		phase_vector_x[2] = np.array(my_var["phase_index_list%s"%transform_matrix][6][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		phase_vector_x = np.array(phase_vector_x)
		phase_vector_x  = np.transpose(phase_vector_x)
		phase_vector_x = np.linalg.pinv(phase_vector_x)

		#Y Direction 
		phase_vector_y = [0]*3
		#+X - +Z
		phase_vector_y[0] = np.array(my_var["phase_index_list%s"%transform_matrix][0][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#-X - +Z
		phase_vector_y[1] = np.array(my_var["phase_index_list%s"%transform_matrix][3][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#+Y+Z - +Z
		phase_vector_y[2] = np.array(my_var["phase_index_list%s"%transform_matrix][7][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		phase_vector_y = np.array(phase_vector_y)
		phase_vector_y  = np.transpose(phase_vector_y)
		phase_vector_y = np.linalg.pinv(phase_vector_y)

		#Z Direction 
		phase_vector_z = [0]*3
		#+X - +Y
		phase_vector_z[0] = np.array(my_var["phase_index_list%s"%transform_matrix][0][1])-np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		#-X - +Y
		phase_vector_z[1] = np.array(my_var["phase_index_list%s"%transform_matrix][3][1])-np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		#+Z+Y - +Y
		phase_vector_z[2] = np.array(my_var["phase_index_list%s"%transform_matrix][7][1])-np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		phase_vector_z = np.array(phase_vector_z)
		phase_vector_z  = np.transpose(phase_vector_z)
		phase_vector_z = np.linalg.pinv(phase_vector_z)

		#R Direction
		phase_vector_r = [0]*3
		#+X
		phase_vector_r[0] = np.array(my_var["phase_index_list%s"%transform_matrix][0][1])
		#+Y
		phase_vector_r[1] = np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		#+Z
		phase_vector_r[2] = np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		phase_vector_r = np.array(phase_vector_r)
		phase_vector_r  = np.transpose(phase_vector_r)
		phase_vector_r = np.linalg.pinv(phase_vector_r)

		vector_tuple = tuple((phase_vector_x,phase_vector_y,phase_vector_z,phase_vector_r))
		ref_list.append(vector_tuple)

	return ref_list

def mvim(frame, data ,ref, theta, axis, pbc_setting, switch):
	property_list = []
	score_list = []
	N = 6
	rotate_data = rotate(theta, axis, data, pbc_setting, switch)
	finder = NearestNeighborFinder(N, rotate_data)
	ptype = data.particles['Particle Type']
	position = data.particles['Position']
	#neighbors[n]第幾近
	#neighbors[n][0]index
	#neighbors[n][1]座標
	#Loop over all input particles:
	for index in range(data.particles.count):
		neighbors = [ (neigh.index, neigh.delta) for neigh in finder.find(index) ]
		neigh_list = [0]*6
		resorted_neighbors_x_ref = sorted( neighbors , key=lambda k: [k[1][0], k[1][1], k[1][2]], reverse=True )
		resorted_neighbors_y_ref = sorted( neighbors , key=lambda k: [k[1][1], k[1][0], k[1][2]], reverse=True )
		resorted_neighbors_z_ref = sorted( neighbors , key=lambda k: [k[1][2], k[1][0], k[1][1]], reverse=True )
		resorted_neighbors_nx_ref = sorted( neighbors , key=lambda k: [k[1][0], k[1][1], k[1][2]], reverse=False )
		resorted_neighbors_ny_ref = sorted( neighbors , key=lambda k: [k[1][1], k[1][0], k[1][2]], reverse=False )
		resorted_neighbors_nz_ref = sorted( neighbors , key=lambda k: [k[1][2], k[1][0], k[1][1]], reverse=False )
		neigh_list[0] = resorted_neighbors_x_ref[0]
		neigh_list[1] = resorted_neighbors_y_ref[0]
		neigh_list[2] = resorted_neighbors_z_ref[0]
		neigh_list[3] = resorted_neighbors_nx_ref[0]
		neigh_list[4] = resorted_neighbors_ny_ref[0]
		neigh_list[5] = resorted_neighbors_nz_ref[0]
		my_var["neigh_list%s"%index] = neigh_list

	for neigh_sort in range(data.particles.count):
		phase_index_list = [0]*8
		#+X
		phase_index_list[0] = my_var["neigh_list%s"%neigh_sort][0]
		#+Y
		phase_index_list[1] = my_var["neigh_list%s"%neigh_sort][1]
		#+Z
		phase_index_list[2] = my_var["neigh_list%s"%neigh_sort][2]
		#-X
		phase_index_list[3] = my_var["neigh_list%s"%neigh_sort][3]
		#-Y
		phase_index_list[4] = my_var["neigh_list%s"%neigh_sort][4]
		#-Z
		phase_index_list[5] = my_var["neigh_list%s"%neigh_sort][5]
		#+X+Z
		phase_index_list[6] = (my_var["neigh_list%s"%my_var['neigh_list%s'%neigh_sort][0][0]][2][0],tuple(((np.array(my_var['neigh_list%s'%my_var['neigh_list%s'%neigh_sort][0][0]][2][1]) + np.array((my_var['neigh_list%s'%neigh_sort][0][1]))).tolist())))
		#+Y+Z
		phase_index_list[7] = (my_var["neigh_list%s"%my_var['neigh_list%s'%neigh_sort][1][0]][2][0],tuple(((np.array(my_var['neigh_list%s'%my_var['neigh_list%s'%neigh_sort][1][0]][2][1]) + np.array((my_var['neigh_list%s'%neigh_sort][1][1]))).tolist())))
		my_var["phase_index_list%s"%neigh_sort] = phase_index_list

	for transform_matrix in range(data.particles.count):
		#X Direction 
		phase_vector_x = [0]*3
		#+Y - +Z
		phase_vector_x[0] = np.array(my_var["phase_index_list%s"%transform_matrix][1][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#-Y - +Z
		phase_vector_x[1] = np.array(my_var["phase_index_list%s"%transform_matrix][4][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#+X+Z - +Z
		phase_vector_x[2] = np.array(my_var["phase_index_list%s"%transform_matrix][6][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		phase_vector_x = np.array(phase_vector_x)
		phase_vector_x  = np.transpose(phase_vector_x)
		#compute F
		m_o_f_x = np.dot(phase_vector_x,ref[transform_matrix][0])
		m_o_f_t_x = np.transpose(m_o_f_x)
		#compute U
		m_o_u2_x = np.dot(m_o_f_t_x,m_o_f_x)
		m_o_u_x = sl.sqrtm(m_o_u2_x)

		score_15_x = np.zeros((1,15),dtype = np.float32)
		m1_score = 0
		m2_score = 0
		m3_score = 0
		m4_score = 0
		o1_score = 0
		o2_score = 0
		a_score_x = 0
		highest_score_x = 0
		variant_x = 0

		if m_o_u_x[0][0] - 1 >= tol:
			score_15_x[0,0] = 1
		if m_o_u_x[1][1] - 1 >= tol :
			score_15_x[0,1] = 1
		if m_o_u_x[2][2] - 1 >= tol :
			score_15_x[0,2] = 1
		if m_o_u_x[0][0] - 1 < ntol :
			score_15_x[0,3] = 1
		if m_o_u_x[1][1] - 1 < ntol :
			score_15_x[0,4] = 1
		if m_o_u_x[2][2] - 1 < ntol :
			score_15_x[0,5] = 1
		if m_o_u_x[0][1] >= tol :
			score_15_x[0,6] = 1
		if m_o_u_x[0][2] >= tol :
			score_15_x[0,7] = 1
		if m_o_u_x[1][2] >= tol :
			score_15_x[0,8] = 1
		if m_o_u_x[0][1] <= ntol :
			score_15_x[0,9] = 1
		if m_o_u_x[0][2] <= ntol :
			score_15_x[0,10] = 1
		if m_o_u_x[1][2] <= ntol :
			score_15_x[0,11] = 1
		if abs(abs(m_o_u_x[0][1]) - abs(m_o_u_x[0][2])) <= 2*tol:
			score_15_x[0,12] = 1
		if abs(abs(m_o_u_x[0][1]) - abs(m_o_u_x[1][2])) <= 2*tol:
			score_15_x[0,13] = 1
		if abs(abs(m_o_u_x[0][2]) - abs(m_o_u_x[1][2])) <= 2*tol:
			score_15_x[0,14] = 1

		#classify
		same_m1 = score_15_x == M1
		m1_score = np.sum(same_m1)
		same_m2 = score_15_x == M2
		m2_score = np.sum(same_m2)
		same_m3 = score_15_x == M3
		m3_score = np.sum(same_m3)
		same_m4 = score_15_x == M4
		m4_score = np.sum(same_m4)
		same_o1 = score_15_x == O1
		o1_score = np.sum(same_o1)
		same_o2 = score_15_x == O2
		o2_score = np.sum(same_o2)
		same_a_x = score_15_x == A
		a_score_x = np.sum(same_a_x)
		highest_score_x = max(m1_score,m2_score,m3_score,m4_score,o1_score,o2_score,a_score_x)
		score_list_x = [m1_score,m2_score,m3_score,m4_score,o1_score,o2_score,a_score_x]

		#determine variants
		if highest_score_x == m1_score:
			variant_x = "M1"
		elif highest_score_x == m2_score:
			variant_x = "M2"
		elif highest_score_x == m3_score:
			variant_x = "M3"
		elif highest_score_x == m4_score:
			variant_x = "M4"
		elif highest_score_x == o1_score:
			variant_x = "O1"
		elif highest_score_x == o2_score:
			variant_x = "O2"
		elif highest_score_x == a_score_x:
			variant_x = "A"
		unknown = np.sum(score_list_x == highest_score_x)
		if unknown != 1:
			variant_x = "unknown"

		#Y Direction 
		phase_vector_y = [0]*3
		#+X - +Z
		phase_vector_y[0] = np.array(my_var["phase_index_list%s"%transform_matrix][0][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#-X - +Z
		phase_vector_y[1] = np.array(my_var["phase_index_list%s"%transform_matrix][3][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		#+Y+Z - +Z
		phase_vector_y[2] = np.array(my_var["phase_index_list%s"%transform_matrix][7][1])-np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		phase_vector_y = np.array(phase_vector_y)
		phase_vector_y  = np.transpose(phase_vector_y)
		#compute F
		m_o_f_y = np.dot(phase_vector_y,ref[transform_matrix][1])
		m_o_f_t_y = np.transpose(m_o_f_y)
		#compute U
		m_o_u2_y = np.dot(m_o_f_t_y,m_o_f_y)
		m_o_u_y = sl.sqrtm(m_o_u2_y)

		score_15_y = np.zeros((1,15),dtype = np.float32)
		m5_score = 0
		m6_score = 0
		m7_score = 0
		m8_score = 0
		o3_score = 0
		o4_score = 0
		a_score_y = 0
		highest_score_y = 0
		variant_y = 0

		if m_o_u_y[0][0] - 1 >= tol:
			score_15_y[0,0] = 1
		if m_o_u_y[1][1] - 1 >= tol :
			score_15_y[0,1] = 1
		if m_o_u_y[2][2] - 1 >= tol :
			score_15_y[0,2] = 1
		if m_o_u_y[0][0] - 1 < ntol :
			score_15_y[0,3] = 1
		if m_o_u_y[1][1] - 1 < ntol :
			score_15_y[0,4] = 1
		if m_o_u_y[2][2] - 1 < ntol :
			score_15_y[0,5] = 1
		if m_o_u_y[0][1] >= tol :
			score_15_y[0,6] = 1
		if m_o_u_y[0][2] >= tol :
			score_15_y[0,7] = 1
		if m_o_u_y[1][2] >= tol :
			score_15_y[0,8] = 1
		if m_o_u_y[0][1] <= ntol :
			score_15_y[0,9] = 1
		if m_o_u_y[0][2] <= ntol :
			score_15_y[0,10] = 1
		if m_o_u_y[1][2] <= ntol :
			score_15_y[0,11] = 1
		if abs(abs(m_o_u_y[0][1]) - abs(m_o_u_y[0][2])) <= 2*tol:
			score_15_y[0,12] = 1
		if abs(abs(m_o_u_y[0][1]) - abs(m_o_u_y[1][2])) <= 2*tol:
			score_15_y[0,13] = 1
		if abs(abs(m_o_u_y[0][2]) - abs(m_o_u_y[1][2])) <= 2*tol:
			score_15_y[0,14] = 1

		#classify
		same_m5 = score_15_y == M5
		m5_score = np.sum(same_m5)
		same_m6 = score_15_y == M6
		m6_score = np.sum(same_m6)
		same_m7 = score_15_y == M7
		m7_score = np.sum(same_m7)
		same_m8 = score_15_y == M8
		m8_score = np.sum(same_m8)
		same_o3 = score_15_y == O3
		o3_score = np.sum(same_o3)
		same_o4 = score_15_y == O4
		o4_score = np.sum(same_o4)
		same_a_y = score_15_y == A
		a_score_y = np.sum(same_a_y)
		highest_score_y = max(m5_score,m6_score,m7_score,m8_score,o3_score,o4_score,a_score_y)
		score_list_y = [m5_score,m6_score,m7_score,m8_score,o3_score,o4_score,a_score_y]
		
		#determine variant_ys
		if highest_score_y == m5_score:
			variant_y = "M5"
		elif highest_score_y == m6_score:
			variant_y = "M6"
		elif highest_score_y == m7_score:
			variant_y = "M7"
		elif highest_score_y == m8_score:
			variant_y = "M8"
		elif highest_score_y == o3_score:
			variant_y = "O3"
		elif highest_score_y == o4_score:
			variant_y = "O4"
		elif highest_score_y == a_score_y:
			variant_y = "A"
		unknown = np.sum(score_list_y == highest_score_y)
		if unknown != 1:
			variant_y = "unknown"

		#Z Direction 
		phase_vector_z = [0]*3
		#+X - +Y
		phase_vector_z[0] = np.array(my_var["phase_index_list%s"%transform_matrix][0][1])-np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		#-X - +Y
		phase_vector_z[1] = np.array(my_var["phase_index_list%s"%transform_matrix][3][1])-np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		#+Z+Y - +Y
		phase_vector_z[2] = np.array(my_var["phase_index_list%s"%transform_matrix][7][1])-np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		phase_vector_z = np.array(phase_vector_z)
		phase_vector_z  = np.transpose(phase_vector_z)
		#compute F
		m_o_f_z = np.dot(phase_vector_z,ref[transform_matrix][2])
		m_o_f_t_z = np.transpose(m_o_f_z)
		#compute U
		m_o_u2_z = np.dot(m_o_f_t_z,m_o_f_z)
		m_o_u_z = sl.sqrtm(m_o_u2_z)

		score_15_z = np.zeros((1,15),dtype = np.float32)
		m9_score = 0
		m10_score = 0
		m11_score = 0
		m12_score = 0
		o5_score = 0
		o6_score = 0
		a_score_z = 0
		highest_score_z = 0
		variant_z = 0

		if m_o_u_z[0][0] - 1 >= tol:
			score_15_z[0,0] = 1
		if m_o_u_z[1][1] - 1 >= tol :
			score_15_z[0,1] = 1
		if m_o_u_z[2][2] - 1 >= tol :
			score_15_z[0,2] = 1
		if m_o_u_z[0][0] - 1 < ntol :
			score_15_z[0,3] = 1
		if m_o_u_z[1][1] - 1 < ntol :
			score_15_z[0,4] = 1
		if m_o_u_z[2][2] - 1 < ntol :
			score_15_z[0,5] = 1
		if m_o_u_z[0][1] >= tol :
			score_15_z[0,6] = 1
		if m_o_u_z[0][2] >= tol :
			score_15_z[0,7] = 1
		if m_o_u_z[1][2] >= tol :
			score_15_z[0,8] = 1
		if m_o_u_z[0][1] <= ntol :
			score_15_z[0,9] = 1
		if m_o_u_z[0][2] <= ntol :
			score_15_z[0,10] = 1
		if m_o_u_z[1][2] <= ntol :
			score_15_z[0,11] = 1
		if abs(abs(m_o_u_z[0][1]) - abs(m_o_u_z[0][2])) <= 2*tol:
			score_15_z[0,12] = 1
		if abs(abs(m_o_u_z[0][1]) - abs(m_o_u_z[1][2])) <= 2*tol:
			score_15_z[0,13] = 1
		if abs(abs(m_o_u_z[0][2]) - abs(m_o_u_z[1][2])) <= 2*tol:
			score_15_z[0,14] = 1

		#classify
		same_m9 = score_15_z == M9
		m9_score = np.sum(same_m9)
		same_m10 = score_15_z == M10
		m10_score = np.sum(same_m10)
		same_m11 = score_15_z == M11
		m11_score = np.sum(same_m11)
		same_m12 = score_15_z == M12
		m12_score = np.sum(same_m12)
		same_o5 = score_15_z == O5
		o5_score = np.sum(same_o5)
		same_o6 = score_15_z == O6
		o6_score = np.sum(same_o6)
		same_a_z = score_15_z == A
		a_score_z = np.sum(same_a_z)
		highest_score_z = max(m9_score,m10_score,m11_score,m12_score,o5_score,o6_score,a_score_z)
		score_list_z = [m9_score,m10_score,m11_score,m12_score,o5_score,o6_score,a_score_z]
		
		#determine variant_zs
		if highest_score_z == m9_score:
			variant_z = "M9"
		elif highest_score_z == m10_score:
			variant_z = "M10"
		elif highest_score_z == m11_score:
			variant_z = "M11"
		elif highest_score_z == m12_score:
			variant_z = "M12"
		elif highest_score_z == o5_score:
			variant_z = "O5"
		elif highest_score_z == o6_score:
			variant_z = "O6"
		elif highest_score_z == a_score_z:
			variant_z = "A"
		unknown = np.sum(score_list_z == highest_score_z)
		if unknown != 1:
			variant_z = "unknown"

		#R Direction
		phase_vector_r = [0]*3
		#+X
		phase_vector_r[0] = np.array(my_var["phase_index_list%s"%transform_matrix][0][1])
		#+Y
		phase_vector_r[1] = np.array(my_var["phase_index_list%s"%transform_matrix][1][1])
		#+Z
		phase_vector_r[2] = np.array(my_var["phase_index_list%s"%transform_matrix][2][1])
		phase_vector_r = np.array(phase_vector_r)
		phase_vector_r  = np.transpose(phase_vector_r)
		#compute F
		r_f = np.dot(phase_vector_r,ref[transform_matrix][3])
		r_f_t = np.transpose(r_f)
		#compute U
		r_u2_z = np.dot(r_f_t,r_f)
		r_u_z = sl.sqrtm(r_u2_z)
		score_15_r = np.zeros((1,15),dtype = np.float32)
		a_r_score = 0
		r1_score = 0
		r2_score = 0
		r3_score = 0
		r4_score = 0
		highest_score_r = 0
		variant_r = 0

		if r_u_z[0][0] - 1 >= tol:
			score_15_r[0,0] = 1
		if r_u_z[1][1] - 1 >= tol :
			score_15_r[0,1] = 1
		if r_u_z[2][2] - 1 >= tol :
			score_15_r[0,2] = 1
		if r_u_z[0][0] - 1 < ntol :
			score_15_r[0,3] = 1
		if r_u_z[1][1] - 1 < ntol :
			score_15_r[0,4] = 1
		if r_u_z[2][2] - 1 < ntol :
			score_15_r[0,5] = 1
		if r_u_z[0][1] >= tol :
			score_15_r[0,6] = 1
		if r_u_z[0][2] >= tol :
			score_15_r[0,7] = 1
		if r_u_z[1][2] >= tol :
			score_15_r[0,8] = 1
		if r_u_z[0][1] <= ntol :
			score_15_r[0,9] = 1
		if r_u_z[0][2] <= ntol :
			score_15_r[0,10] = 1
		if r_u_z[1][2] <= ntol :
			score_15_r[0,11] = 1
		if abs(abs(r_u_z[0][1]) - abs(r_u_z[0][2])) <= 2*tol:
			score_15_r[0,12] = 1
		if abs(abs(r_u_z[0][1]) - abs(r_u_z[1][2])) <= 2*tol:
			score_15_r[0,13] = 1
		if abs(abs(r_u_z[0][2]) - abs(r_u_z[1][2])) <= 2*tol:
			score_15_r[0,14] = 1

		same_a = score_15_r == A
		a_r_score = np.sum(same_a)
		same_r1 = score_15_r == R1
		r1_score = np.sum(same_r1)
		same_r2 = score_15_r == R2
		r2_score = np.sum(same_r2)
		same_r3 = score_15_r == R3
		r3_score = np.sum(same_r3)
		same_r4 = score_15_r == R4
		r4_score = np.sum(same_r4)

		highest_score_r = max(a_r_score, r1_score, r2_score, r3_score, r4_score)
		score_list_r = [a_r_score, r1_score, r2_score, r3_score, r4_score]

		#determine variants
		if highest_score_r == a_r_score:
			variant_r = "A"
		elif highest_score_r == r1_score:
			variant_r = 'R1'
		elif highest_score_r == r2_score:
			variant_r = 'R2'
		elif highest_score_r == r3_score:
			variant_r = 'R3'
		elif highest_score_r == r4_score:
			variant_r = 'R4'

		unknown = np.sum(score_list_r == highest_score_r)
		if unknown != 1:
			variant_r = "unknown"

		#Determine final variant
		variant = 0
		variant_color = 0
		highest_score_list = [highest_score_x,highest_score_y,highest_score_z,highest_score_r]
		highest_score = max(highest_score_list)
		unknown = np.sum(highest_score_list == highest_score)

		if highest_score == highest_score_x:
			variant = variant_x
		elif highest_score == highest_score_z:
			variant = variant_z
		elif highest_score == highest_score_y:
			variant = variant_y
		elif highest_score == highest_score_r:
			variant = variant_r

		if variant != "A" and unknown != 1:
			variant = "unknown"
		if highest_score < 9 :
			variant = "unknown"
		#color
		if variant == "M1":
			variant_color = 0
		elif variant == "M2":
			variant_color = 1
		elif variant == "M3":
			variant_color = 2
		elif variant == "M4":
			variant_color = 3
		elif variant == "M5":
			variant_color = 4
		elif variant == "M6":
			variant_color = 5
		elif variant == "M7":
			variant_color = 6
		elif variant == "M8":
			variant_color = 7
		elif variant == "M9":
			variant_color = 8
		elif variant == "M10":
			variant_color = 9
		elif variant == "M11":
			variant_color = 10
		elif variant == "M12":
			variant_color = 11
		elif variant == "O1":
			variant_color = 12
		elif variant == "O2":
			variant_color = 13
		elif variant == "O3":
			variant_color = 14
		elif variant == "O4":
			variant_color = 15
		elif variant == "O5":
			variant_color = 16
		elif variant == "O6":
			variant_color = 17
		elif variant == "R1":
			variant_color = 18
		elif variant == "R2":
			variant_color = 19
		elif variant == "R3":
			variant_color = 20
		elif variant == "R4":
			variant_color = 21
		elif variant == "A":
			variant_color = 22
		elif variant == "unknown":
			variant_color = 23
		property_list.append(variant_color)
		score_list.append(highest_score)
	
	
	return property_list, score_list

data_ref_1 = node_ref_1.compute()
data_ref_2 = node_ref_2.compute()
ref_list_type1 = compute_ref_vector(data_ref_1, rotation_angle, rotation_axis, pbc_setting, switch)
ref_list_type2 = compute_ref_vector(data_ref_2, rotation_angle, rotation_axis, pbc_setting, switch)

if switch == True:
	data_ref = node_ref.compute()
	rotated_data_ref = rotate(rotation_angle, rotation_axis, data_ref, pbc_setting, switch)
	export_file(rotated_data_ref, '/data2/paper_oo/set/ref_y.dump', 'lammps/dump',columns = ['Particle Identifier','Particle Type','Position.X','Position.Y','Position.Z'])

for frame_index in range(0, node.source.num_frames, 1):
	data = node.compute(frame_index)
	data_1 = node_1.compute(frame_index)
	data_2 = node_2.compute(frame_index)
	mvim_1 = mvim(frame_index, data_1, ref_list_type1, rotation_angle, rotation_axis, pbc_setting, switch)
	mvim_type1 = mvim_1[0]
	mvim_score1 = mvim_1[1]
	mvim_2 = mvim(frame_index, data_2, ref_list_type2, rotation_angle, rotation_axis, pbc_setting, switch)
	mvim_type2 = mvim_2[0]
	mvim_score2 = mvim_2[1]
	mvim_type = mvim_type1 + mvim_type2
	mvim_score = mvim_score1 + mvim_score2
	rotated_data = rotate(rotation_angle, rotation_axis, data, pbc_setting, switch)
	rotated_data.particles_.create_property('Variant', data = mvim_type)
	rotated_data.particles_.create_property('Score', data = mvim_score)
	
	export_file(rotated_data, '/data2/post_50_y_%s.dump'%frame_index, 'lammps/dump',columns = ['Particle Identifier','Particle Type','Position.X','Position.Y','Position.Z','Variant','Score'],frame = frame_index)

endtime = datetime.datetime.now()
print (endtime - starttime)