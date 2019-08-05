from ovito.io import import_file
from ovito.io import export_file
from ovito.modifiers import *
from ovito.data import *
import numpy as np
import os

# inputfile
node = import_file('/Users/eason/Desktop/U/c')
node.add_to_scene()

my_var = {}

for frame_index in range(node.source.num_frames):
    node.modifiers.append(SelectParticleTypeModifier(property='Particle Type', types={2}))
    node.modifiers.append(DeleteSelectedParticlesModifier())
    data = node.compute(frame_index)
    positions = data.particles_['Position_']
    my_var["position_frame%s"%frame_index] = positions[:]


position_total = np.zeros(my_var["position_frame0"].shape)

count = 0
series = 0
num_average = 50
num_average_reci = 1/num_average

for i in range(node.source.num_frames):
    count = count + 1
    position_total = position_total + my_var["position_frame%s"%i]

    if count%num_average == 0:
        my_var["average_position%s"%series] = position_total.copy()
        my_var["average%s"%series] = my_var["average_position%s"%series] * [num_average_reci]

        series = series + 1
        position_total = np.zeros(my_var["position_frame0"].shape)
        count = 0

for j in range(series):
    def compute_myproperty(frame, data):
        data.particles_.create_property('Position', data = my_var["average%s"%j])
    node.modifiers.append(PythonScriptModifier(function = compute_myproperty))
 
    data = node.compute()

    f = open('/Users/eason/Desktop/U/average.%s'%j,'w')
    export_file(node, '/Users/eason/Desktop/U/average.%s'%j, 'imd',columns = ['Particle Identifier','Particle Type','Position'],frame = j)
    f.close()














