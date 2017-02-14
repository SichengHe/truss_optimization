# Run script for the continuous 3-D wing-truss optimization problem
# Dr. John T. Hwang
# Sicheng He
# June, 2016

from __future__ import division
import numpy
import time
import fractions as fractions

from openmdao.api import IndepVarComp, Problem, Group, ScipyOptimizer, SqliteRecorder, pyOptSparseDriver
from component import SysDispAug, SysDisplacements, SysCompliance, SysVolume, SysStress, SysKS, EulerBucklingKS
from utils import setup_problem, writeBDF, getLength
from openmdao.devtools.partition_tree_n2 import view_tree


# other param (not subject to change...)
u1 = 0.005
u2 = 0.995

def case_preprocess(ind_case):

    # cantilevered truss case
    if (ind_case==1):

        pattern_length = 1.0

        N_pattern = 2

        # nodes
        nodes = []
        for i in xrange(N_pattern+1):
            for j in xrange(2):
                nodes_loc = [i*pattern_length,j*-pattern_length]
                nodes.append(nodes_loc)
        nodes = numpy.array(nodes)

        print(nodes)


        #
        elements = []
        for i in xrange(N_pattern):
            loc_ind = i*2
            elements.append([loc_ind+1,loc_ind+3])
            elements.append([loc_ind,  loc_ind+3])
            elements.append([loc_ind+2,loc_ind+3])
            elements.append([loc_ind+1,loc_ind+2])
            elements.append([loc_ind,  loc_ind+2])


        elements = numpy.array(elements)

        #
        cons = []
        cons.append(0)
        cons.append(1)
        cons = numpy.array(cons)

        #
        forces = numpy.zeros((N_pattern*2+2,2))
        forces[-1,1] = 1e5

        return [nodes, elements, cons, forces]

    # Michell structure
    elif (ind_case==2):

        n_node_x = 8
        n_node_y = 5
        n_node_conn_x = n_node_x-1
        n_node_conn_y = n_node_y-1

        ind_ld_pt_x = n_node_x-1
        ind_ld_pt_y = n_node_y/2
        ind_ld_pt = ind_ld_pt_x*n_node_y+ind_ld_pt_y

        N_node = n_node_x*n_node_y

        # node

        nodes = []

        for i in xrange(n_node_x):
            for j in xrange(n_node_y):
                coord = [i*1.0, j*-1.0, 0.0]
                nodes.append(coord)

        nodes = numpy.array(nodes)

        # elements

        elements = []

        for i in xrange(n_node_x):
            for j in xrange(n_node_y):

                ind_node_low = i*n_node_y+j

                dn_x = n_node_x-i
                dn_y = n_node_y-j

                dn_x = min(dn_x,n_node_conn_x+1)
                dn_y = min(dn_y,n_node_conn_y+1)

                #print(dn_x,dn_y)

                for ii in xrange(dn_x):
                    for jj in xrange(dn_y):

                        if (fractions.gcd(ii,jj)!=1): # exclude the overlap bars
                            pass
                        elif (i==0 and ii==0): # exclude the overconstrained bar
                            pass
                        elif (jj==0 and ii!=1):
                            pass
                        elif (ii==0 and jj!=1):
                            pass
                        else:
                            ind_node_up_x = i+ii
                            ind_node_up_y = j+jj

                            #print(ind_node_up_x,ind_node_up_y)

                            ind_node_up = ind_node_up_x*n_node_y+ind_node_up_y

                            elem_loc = [ind_node_low,ind_node_up]

                            elements.append(elem_loc)

                dn_y = j+1
                dn_y = min(dn_y,n_node_conn_y+1)

                for ii in xrange(dn_x):
                    for jj in xrange(dn_y):

                        if (fractions.gcd(ii,jj)!=1): # exclude the overlap bars
                            pass
                        elif (jj==0):
                            pass
                        elif (ii==0):
                            pass
                        else:
                            ind_node_up_x = i+ii
                            ind_node_up_y = j-jj

                            ind_node_up = ind_node_up_x*n_node_y+ind_node_up_y

                            elem_loc = [ind_node_low,ind_node_up]

                            elements.append(elem_loc)

        elements = numpy.array(elements)

        cons = []

        for i in xrange(n_node_y):
            cons.append(i)

        cons = numpy.array(cons)

        forces = numpy.zeros([N_node,3])
        forces[ind_ld_pt,1] = 4.0*1e6



        return [nodes, elements, cons, forces]


    # wing case
    elif (ind_case==3):


        geom_file = '../CRM_AVL/wing_coarse.avl'
        upp_file = '../airfoils/rae2822_upp.txt'
        low_file = '../airfoils/rae2822_low.txt'
        results_file = '../CRM_AVL/results_coarse.txt'

        factor = 1

        xyz, nodes, elements, cons, forces, forcesArray = setup_problem(u1, u2, geom_file, upp_file, low_file, results_file,factor)
        forces /= 1.e0

        return [nodes, elements, cons, forces]

    



def scaling(ind_case, scale_x, scale_mass, E, s0, area_low, area_up):

    # scaled
    E = E/scale_x*scale_mass
    s0 = s0/scale_x*scale_mass


    # set up the case

    [nodes, elements, cons, forces] = case_preprocess(ind_case)


    # scaled load and geometry
    forces = forces*scale_mass*scale_x
    nodes = nodes*scale_x

    elem_len = getLength(nodes,elements)



    # scaled area bd
    area_low = area_low * scale_x**2
    area_up = area_up * scale_x**2

    return E, s0, nodes, elements, cons, forces, area_low, area_up, elem_len




# case to runpip install git+http://github.com/OpenMDAO/OpenMDAO.git@master
ind_case = 1

# scaling factor
scale_x = 1e2
scale_mass = 1e-6

# Youngs modulus and yield stress
E = 69.0e9
s0 = 270.0e6

# area bd
area_low = 0.0001
area_up = 0.1


E, s0, nodes, elements, cons, forces, area_low, area_up, elem_len = scaling(ind_case, scale_x, scale_mass, E, s0, area_low, area_up)



# save the files
numpy.savetxt('../output/data_elems.dat', elements)
numpy.savetxt('../output/data_nodes.dat', nodes)
numpy.savetxt('../output/data_constraints.dat', cons)
numpy.savetxt('../output/data_forces.dat', forces)


# set up OpenMDAO systems and solve the optimzation problem
root = Group()
root.add('comp_areas',
         IndepVarComp([('areas', 4.e-3 * numpy.ones(elements.shape[0]))]),
         promotes=['*'])
root.add('sys_disp_aug',
         SysDispAug(nodes, elements, forces.flatten(), cons, E),
         promotes=['*'])
root.add('sys_displacements',
         SysDisplacements(nodes, cons),
         promotes=['*'])
root.add('sys_compliance',
         SysCompliance(nodes, forces),
         promotes=['*'])
root.add('sys_volume',
         SysVolume(elements, elem_len),
         promotes=['*'])
root.add('sys_stress',
         SysStress(nodes, elements, E, s0),
         promotes=['*'])
root.add('sys_ks',
         SysKS(elements, 1.0),
         promotes=['*'])
# root.add('bkl_ks',
#          EulerBucklingKS(E, elements, nodes, cons),
#          promotes=['*'])

prob = Problem()
prob.root = root
#prob.root.deriv_options['type'] = 'fd'
prob.setup()

t0 = time.time()
prob.run()
t1 = time.time()

#print t1-t0

nodes0 = nodes
nodes1 = nodes + prob['disp']

#writeBDF('jig.bdf', nodes0, elements+1)
#writeBDF('deflected.bdf', nodes1, elements+1)

if 0:
    prob.check_partial_derivatives(compact_print=True)
    exit()

prob.driver = pyOptSparseDriver()
prob.driver.options['optimizer'] = "SNOPT"
prob.driver.opt_settings = {'Major optimality tolerance': 1.0e-7,
                            'Major feasibility tolerance': 1.0e-7,
                            'Iterations limit': int(1e8),
}

prob.driver.add_desvar('areas',lower=area_low, upper=area_up, scaler=1e0) # test
prob.driver.add_objective('volume', scaler=1e0)
prob.driver.add_constraint('minstress', upper=0.)
prob.driver.add_constraint('maxstress', upper=0.)
#\prob.driver.add_constraint('neg_stress_plus_buckling_con', upper=0.)

# setup data recording
recorder = SqliteRecorder('postprocess/data.db')
recorder.options['record_params'] = True
recorder.options['record_metadata'] = True
prob.driver.add_recorder(recorder)

prob.setup()
view_tree(prob, outfile="partition_tree/aerostruct.html", show_browser=True)
prob.run()

# print("prob['areas']",prob['areas'])
print("prob['stress']",prob['stress'])
# print("prob['disp_aug']",prob['disp_aug'])
# print("prob['disp']",prob['disp'])
# print("neg_stress_plus_buckling_con",prob['neg_stress_plus_buckling_con'])

#writeBDF('optimized.bdf', nodes+prob['disp'], elements+1)
