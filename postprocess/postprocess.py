from __future__ import division
import numpy as numpy
import pylab
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from scipy import *
from math import *
import sqlitedict

s0 = 270.0e6 #######

def plot(ax, sol, elems, nodes):
    for i in xrange(len(sol)):
        sol_loc = sol[i]

        ind_node_1 = elems[i][0]
        ind_node_2 = elems[i][1]

        if (1 == 0):

            x1 = nodes[ind_node_1][0]
            y1 = nodes[ind_node_1][1]
            z1 = nodes[ind_node_1][2]

            x2 = nodes[ind_node_2][0]
            y2 = nodes[ind_node_2][1]
            z2 = nodes[ind_node_2][2]

            x = [x1, x2]
            y = [y1, y2]
            z = [z1, z2]

        else:

            x1 = nodes[ind_node_1][0]
            y1 = nodes[ind_node_1][1]

            x2 = nodes[ind_node_2][0]
            y2 = nodes[ind_node_2][1]

            x = [x1, x2]
            y = [y1, y2]
        #print '++++++',z

        # ax.plot(x,y,z,'-ro',
        #         #c=colorVal)
        #         markersize=0.1)

        if (1 == 0):
            ax.plot(x,y,z,'-k',
                    #c=colorVal)
                    linewidth=sol_loc)
        else:

            ax.plot(x,y,'-k',
                    #c=colorVal)
                    linewidth=sol_loc)

        i += 1


    ax.view_init(elev=90., azim=-90)
    ax.grid(b=False)
    ax.set_aspect(1)




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

nodes = loadtxt('../output/data_nodes.dat')
elems = loadtxt('../output/data_elems.dat')


db = sqlitedict.SqliteDict('data.db', 'iterations')
print( list( db.keys() ) )




index = 0


for case_name, case_data in db.iteritems():
    if "metadata" in case_name or "derivs" in case_name:
        continue # don't plot these cases

    area = case_data['Unknowns']['areas']
    stress = case_data['Unknowns']['stress']

    if 1 and index % 10 == 0:
        max_abs_area = max(area)
        dimless_area = area/max_abs_area

        plot(ax, dimless_area, elems, nodes)
        plt.axis('off')
        pylab.savefig('tmp/area%03i.png'%(index))

        abs_stress = abs(stress)
        max_abs_stress = max(abs_stress)
        dimless_stress = abs_stress/max_abs_stress
        for i in xrange(len(dimless_stress)):
            if dimless_stress[i]>0.999:
                dimless_stress[i] = 5.0

        plot(ax, dimless_stress, elems, nodes)
        plt.axis('off')
        pylab.savefig('tmp/stress%03i.png'%(index))
        print index

    print index
    index += 1

