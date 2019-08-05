from ovito.io import import_file
import numpy as np
import matplotlib.pyplot as plt

# Load input simulation file.
node = import_file('/Users/eason/Downloads/drive-download-20190123T185018Z-001/post.*')
my_var = {}
M1list = []
M2list = []
M3list = []
M4list = []
M5list = []
M6list = []
M7list = []
M8list = []
M9list = []
M10list = []
M11list = []
M12list = []
O1list = []
O2list = []
O3list = []
O4list = []
O5list = []
O6list = []
Alist = []
unknownlist = []
for frame_index in range(node.source.num_frames):
	data = node.compute(frame_index)
	variant = data.particles['Variant']
	variant2 = list(variant)
	total = data.particles.count
	my_var["total%s"%frame_index] = data.particles.count
	my_var["M1%s"%frame_index] = variant2.count(0)/total
	my_var["M2%s"%frame_index] = variant2.count(1)/total
	my_var["M3%s"%frame_index] = variant2.count(2)/total
	my_var["M4%s"%frame_index] = variant2.count(3)/total
	my_var["M5%s"%frame_index] = variant2.count(4)/total
	my_var["M6%s"%frame_index] = variant2.count(5)/total
	my_var["M7%s"%frame_index] = variant2.count(6)/total
	my_var["M8%s"%frame_index] = variant2.count(7)/total
	my_var["M9%s"%frame_index] = variant2.count(8)/total
	my_var["M10%s"%frame_index] = variant2.count(9)/total
	my_var["M11%s"%frame_index] = variant2.count(10)/total
	my_var["M12%s"%frame_index] = variant2.count(11)/total
	my_var["O1%s"%frame_index] = variant2.count(12)/total
	my_var["O2%s"%frame_index] = variant2.count(13)/total
	my_var["O3%s"%frame_index] = variant2.count(14)/total
	my_var["O4%s"%frame_index] = variant2.count(15)/total
	my_var["O5%s"%frame_index] = variant2.count(16)/total
	my_var["O6%s"%frame_index] = variant2.count(17)/total
	my_var["A%s"%frame_index] = variant2.count(18)/total
	my_var["unknown%s"%frame_index] = variant2.count(19)/total
for j in range(node.source.num_frames):
	M1list.append(my_var["M1%s"%j])
	M2list.append(my_var["M2%s"%j])
	M3list.append(my_var["M3%s"%j])
	M4list.append(my_var["M4%s"%j])
	M5list.append(my_var["M5%s"%j])
	M6list.append(my_var["M6%s"%j])
	M7list.append(my_var["M7%s"%j])
	M8list.append(my_var["M8%s"%j])
	M9list.append(my_var["M9%s"%j])
	M10list.append(my_var["M10%s"%j])
	M11list.append(my_var["M11%s"%j])
	M12list.append(my_var["M12%s"%j])
	O1list.append(my_var["O1%s"%j])
	O2list.append(my_var["O2%s"%j])
	O3list.append(my_var["O3%s"%j])
	O4list.append(my_var["O4%s"%j])
	O5list.append(my_var["O5%s"%j])
	O6list.append(my_var["O6%s"%j])
	Alist.append(my_var["A%s"%j])
	unknownlist.append(my_var["unknown%s"%j])
x = np.linspace(0,1200,46)
#plt.plot(x,M1list,color = '#DEB887',label = 'M1')
#plt.plot(x,M2list,color = '#CD853F',label = 'M2')
#plt.plot(x,M3list,color = '#A0522D',label = 'M3')
#plt.plot(x,M4list,color = '#6C320A',label = 'M4')
#plt.plot(x,M5list,color = '#FFFF00',label = 'M5')
#plt.plot(x,M6list,color = '#FFD700',label = 'M6')
#plt.plot(x,M7list,color = '#EEDD82',label = 'M7')
#plt.plot(x,M8list,color = '#DAA520',label = 'M8')
plt.plot(x,M9list,color = '#92D050',label = 'M9')
plt.plot(x,M10list,color = '#385723',label = 'M10')
plt.plot(x,M11list,color = '#9AFF9A',label = 'M11')
plt.plot(x,M12list,color = '#548B54',label = 'M12')
#plt.plot(x,O1list,color = '#836FFF',label = 'O1')
#plt.plot(x,O2list,color = '#473C8B',label = 'O2')
#plt.plot(x,O3list,color = '#0000FF',label = 'O3')
#plt.plot(x,O4list,color = '#00008B',label = 'O4')
plt.plot(x,O5list,color = '#4876FF',label = 'O5')
plt.plot(x,O6list,color = '#27408B',label = 'O6')
plt.plot(x,Alist,color = '#C00000',label = 'A')
plt.plot(x,unknownlist,color = '#1E1E1E',label = 'Other')
plt.xlabel('time(ps)', fontsize = 32)
plt.ylabel('propotion', fontsize = 32)
plt.xticks([0.0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200], fontsize = 32)
plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize = 32)
plt.legend(loc='best', fontsize = 28)
plt.show()