# OpenMDAO components for the continuous 3-D wing-truss optimization problem
# Dr. John T. Hwang
# Sicheng He
# June, 2016

from __future__ import division
import numpy
import scipy.sparse
import time
from scipy.sparse.linalg import gmres, LinearOperator, aslinearoperator

from openmdao.api import Component

import lib

from utils import getLength, getkFactor


class SysDispAug(Component):

    def __init__(self, nodes, elements, loads, cons, E):
        '''
        nodes : ndarray[num_nodes, 3], float
        elements : ndarray[num_elems, 2], int
        loads : ndarray[num_nodes, 3], float
        cons : ndarray[num_cons], int
        E : float
        '''
        super(SysDispAug, self).__init__()

        num_nodes = nodes.shape[0]
        num_elem = elements.shape[0]
        num_cons = cons.shape[0]

        self.add_param('areas', val=numpy.zeros(num_elem))

        if (1 == 0):

            self.add_state('disp_aug', val=numpy.zeros(3*num_nodes + 3*num_cons))

        else:

            self.add_state('disp_aug', val=numpy.zeros(2*num_nodes + 2*num_cons))


        self.nodes = nodes
        self.elements = elements
        self.loads = loads
        self.cons = cons

        self.E = E


        if (1 == 0):

            self.rhs = numpy.zeros(3*num_nodes + 3*num_cons)
            self.rhs[:3*num_nodes] = loads

        else:

            self.rhs = numpy.zeros(2*num_nodes + 2*num_cons)
            self.rhs[:2*num_nodes] = loads


        nodes = self.nodes
        elems = self.elements
        cons = self.cons
        E = self.E

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        num_cons = cons.shape[0]

        if (1 == 0):

            nnz = 36 * num_elems + 2 * 3 * num_cons

        else:

            nnz = 16 * num_elems + 2 * 2 * num_cons

        if (1 == 0):

            data, self.rows, self.cols = lib.getmtx(num_nodes, num_elems, num_cons, nnz,
                                          E, nodes, elems+1, numpy.ones(num_elems), cons+1)
        else:

            data, self.rows, self.cols = lib.getmtx2d(num_nodes, num_elems, num_cons, nnz,
                                          E, nodes, elems+1, numpy.ones(num_elems), cons+1)


        # if (1 == 0):

        #     nnz = 36 * num_elems

        #     out = lib.getresder(num_nodes, num_elems, nnz, E, nodes, elems+1)
        #     self.data2, self.rows2, self.cols2, self.ind_aug2 = out


        # else:

        #     nnz = 16 * num_elems

        #     out = lib.getresder(num_nodes, num_elems, nnz, E, nodes, elems+1)
        #     self.data2, self.rows2, self.cols2, self.ind_aug2 = out

#        self.deriv_options['type'] = 'cs'
#        self.deriv_options['form'] = 'central'


    def get_matrix(self, areas):

        nodes = self.nodes
        elems = self.elements
        cons = self.cons
        E = self.E

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        num_cons = cons.shape[0]

        if (1 == 0):

            nnz = 36 * num_elems + 2 * 3 * num_cons

            data = lib.getmtx2(num_nodes, num_elems, nnz,
                           E, nodes, elems+1, areas)

            size = 3 * num_nodes + 3 * num_cons

        else:

            nnz = 16 * num_elems + 2 * 2 * num_cons

            data = lib.getmtx2d2(num_nodes, num_elems, nnz,
                           E, nodes, elems+1, areas)

            size = 2 * num_nodes + 2 * num_cons


        mat = scipy.sparse.csc_matrix((data, (self.rows, self.cols)),
                                      shape=(size, size))

        return mat

    def apply_nonlinear(self, params, unknowns, resids):
        self.mat = self.get_matrix(params['areas'])

        resids['disp_aug'] = self.mat.dot(unknowns['disp_aug']) - self.rhs

    def callback(self, res):
        #print self.counter, numpy.linalg.norm(res)
        self.counter += 1

    def gmres(self, mat, rhs, lu):
        if 1:
            size = len(rhs)
            A = aslinearoperator(mat)
            M = LinearOperator((size, size), dtype=float, matvec=lu.solve)
            self.counter = 0
            sol, info = gmres(A, rhs, M=M, maxiter=10,
                              #callback=self.callback,
                              tol=1e-12)
            return sol
        else:
            return lu.solve(rhs)

    def solve_nonlinear(self, params, unknowns, resids):
        self.mat = self.get_matrix(params['areas'])
        self.lu = scipy.sparse.linalg.splu(self.mat)

        #unknowns['disp_aug'] = self.lu.solve(self.rhs)
        unknowns['disp_aug'] = self.gmres(self.mat, self.rhs, self.lu)

    def linearize(self, params, unknowns, resids):
        self.lu_T = scipy.sparse.linalg.splu(self.mat.T)

        nodes = self.nodes
        elems = self.elements
        cons = self.cons
        disp_aug = unknowns['disp_aug']
        E = self.E

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        num_cons = cons.shape[0]
        num_aug = disp_aug.shape[0]

        if (1 == 0):
            nnz = 36 * num_elems

            data, rows, cols = lib.getresder(num_nodes, num_elems, num_aug, nnz,
                                            E, nodes, elems+1, disp_aug)

            size = 3 * num_nodes + 3 * num_cons

        else:
            nnz = 16 * num_elems

            data, rows, cols = lib.getresder2d(num_nodes, num_elems, num_aug, nnz,
                                            E, nodes, elems+1, disp_aug)

            size = 2 * num_nodes + 2 * num_cons


        mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                      shape=(size, num_elems))

        jac = {} # self.alloc_jacobian()
        jac['disp_aug', 'disp_aug'] = self.mat
        jac['disp_aug', 'areas'] = mat

        return jac

    def solve_linear(self, dumat, drmat, vois, mode=None):
        if mode == 'fwd':
            sol_vec, rhs_vec = self.dumat, self.drmat
            t = 0
        else:
            sol_vec, rhs_vec = self.drmat, self.dumat
            t = 1

        for voi in vois:
            if mode == 'fwd':
                sol_vec[voi].vec[:] = self.lu.solve(rhs_vec[voi].vec)
            else:
                #sol_vec[voi].vec[:] = self.lu_T.solve(
                sol_vec[voi].vec[:] = self.gmres(self.mat.T, rhs_vec[voi].vec,
                                                 self.lu_T)



class SysDisplacements(Component):
    """ Selects displacements from augmented vector """

    def __init__(self, nodes, cons):
        super(SysDisplacements, self).__init__()

        n = nodes.shape[0]

        if (1 == 0):
            size = 3 * n + 3 * cons.shape[0]
        else:
            size = 2 * n + 2 * cons.shape[0]

        self.n = n

        self.add_param('disp_aug', val=numpy.zeros((size)))

        if (1 == 0):
            self.add_output('disp', val=numpy.zeros((n, 3)))
        else:
            self.add_output('disp', val=numpy.zeros((n, 2)))


        #self.deriv_options['type'] = 'cs'
        #self.deriv_options['form'] = 'central'
        #self.deriv_options['extra_check_partials_form'] = "central"
        if (1 == 0):
            data = numpy.ones(3 * n)
            lins = numpy.arange(3 * n)
            self.mat = scipy.sparse.csc_matrix((data, (lins, lins)),
                                           shape=(3*n, size))
        else:
            data = numpy.ones(2 * n)
            lins = numpy.arange(2 * n)
            self.mat = scipy.sparse.csc_matrix((data, (lins, lins)),
                                           shape=(2*n, size))


    def solve_nonlinear(self, params, unknowns, resids):
        n = self.n
        if (1 == 0):
            unknowns['disp'] = numpy.array(params['disp_aug'][:3*n].reshape((n, 3)))
        else:
            unknowns['disp'] = numpy.array(params['disp_aug'][:2*n].reshape((n, 2)))


    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['disp', 'disp_aug'] = self.mat
        return jac



class SysCompliance(Component):

    def __init__(self, nodes, loads):
        super(SysCompliance, self).__init__()

        self.n = nodes.shape[0]
        self.loads = loads

        if (1 == 0):

            self.add_param('disp', val=numpy.zeros((self.n, 3)))

        else:

            self.add_param('disp', val=numpy.zeros((self.n, 2)))

        self.add_output('compliance', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['compliance'] = numpy.sum(params['disp'] * self.loads)

    def linearize(self, params, unknowns, resids):

        jac = {}

        if (1 == 0):
            jac['compliance', 'disp'] = self.loads.reshape((1, 3*self.n))
        else:
            jac['compliance', 'disp'] = self.loads.reshape((1, 2*self.n))

        return jac



class SysVolume(Component):

    def __init__(self, elems, lengths):
        super(SysVolume, self).__init__()

        self.n = elems.shape[0]
        self.lengths = lengths

        self.add_param('areas', val=numpy.zeros(self.n))
        self.add_output('volume', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['volume'] = numpy.sum(params['areas'] * self.lengths)

    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['volume', 'areas'] = self.lengths.reshape((1, self.n))
        return jac




class SysStress(Component):

    def __init__(self, nodes, elems, E, s0):
        super(SysStress, self).__init__()

        self.nodes = nodes
        self.elems = elems
        self.E = E
        self.s0 = s0

        if (1 == 0):

            self.add_param('disp', val=numpy.zeros((nodes.shape[0], 3)))

        else:

            self.add_param('disp', val=numpy.zeros((nodes.shape[0], 2)))

        self.add_output('stress', val=numpy.zeros(elems.shape[0]))

        nodes = self.nodes
        elems = self.elems
        E = self.E
        s0 = self.s0

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]

        if (1 == 0):

            nnz = 2 * 3 * num_elems

            data, rows, cols = lib.getstressder(num_nodes, num_elems, nnz,
                                            E, nodes, elems+1)
            data /= s0

            self.mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                               shape=(num_elems, 3 * num_nodes))

        else:
            
            nnz = 2 * 2 * num_elems

            data, rows, cols = lib.getstressder2d(num_nodes, num_elems, nnz,
                                            E, nodes, elems+1)
            data /= s0

            self.mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                               shape=(num_elems, 2 * num_nodes))

    def solve_nonlinear(self, params, unknowns, resids):
        nodes = self.nodes
        elems = self.elems
        E = self.E
        s0 = self.s0

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]

        #unknowns['stress'] = lib.getstresses(num_nodes, num_elems,
        #                                     E, nodes, elems+1, params['disp']) / s0
        unknowns['stress'] = self.mat.dot(params['disp'].flatten())

    def linearize(self, params, unknowns, resids):
        '''
        nodes = self.nodes
        elems = self.elems
        E = self.E
        s0 = self.s0

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        nnz = 2 * 3 * num_elems

        data, rows, cols = lib.getstressder(num_nodes, num_elems, nnz,
                                            E, nodes, elems+1)
        data /= s0
        mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                      shape=(num_elems, 3 * num_nodes))
        '''

        jacs = {}
        jacs['stress', 'disp'] = self.mat

        return jacs



class SysKS(Component):
    """ Aggregates failure constraints from the structure """

    def __init__(self, elems, s0, rho=10**8):
        super(SysKS, self).__init__()

        self.n = elems.shape[0]
        self.s0 = s0
        self.rho = rho

        self.add_param('stress', val=numpy.zeros(self.n))
        self.add_output('minstress', val=0.)
        self.add_output('maxstress', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        s0 = self.s0
        rho = self.rho
        stress = params['stress']

        fmin = -s0 - stress
        fmax =  stress - s0

        maxfmin = numpy.max(fmin)
        maxfmax = numpy.max(fmax)

        unknowns['minstress'] = maxfmin + 1 / rho * \
                                numpy.log(numpy.sum(numpy.exp(rho*(fmin - maxfmin))))
        unknowns['maxstress'] = maxfmax + 1 / rho * \
                                numpy.log(numpy.sum(numpy.exp(rho*(fmax - maxfmax))))

    def linearize(self, params, unknowns, resids):
        n = self.n
        s0 = self.s0
        rho = self.rho
        stress = params['stress']

        fmin = -s0 - stress
        fmax =  stress - s0
        sgnmin = -1.0
        sgnmax =  1.0

        indmin = numpy.argmax(fmin)
        indmax = numpy.argmax(fmax)
        maxfmin = fmin[indmin]
        maxfmax = fmax[indmax]

        dmaxfmin_ds = numpy.zeros(self.n)
        dmaxfmax_ds = numpy.zeros(self.n)
        dmaxfmin_ds[indmin] = sgnmin * 1.0
        dmaxfmax_ds[indmax] = sgnmax * 1.0

        datamin = dmaxfmin_ds + 1/rho * 1.0 /\
                  numpy.sum(numpy.exp(rho * (fmin - maxfmin))) * \
                  numpy.exp(rho * (fmin - maxfmin)) * (sgnmin * rho)
        datamax = dmaxfmax_ds + 1/rho * 1.0 /\
                  numpy.sum(numpy.exp(rho * (fmax - maxfmax))) * \
                  numpy.exp(rho * (fmax - maxfmax)) * (sgnmax * rho)
        datamin[indmin] -= 1/rho * sgnmin * rho
        datamax[indmax] -= 1/rho * sgnmax * rho

        jacs = {}
        jacs['minstress', 'stress'] = datamin.reshape((1, self.n))
        jacs['maxstress', 'stress'] = datamax.reshape((1, self.n))
        return jacs

class EulerBucklingKS(Component):

    """Euler buckling constraint"""
    def __init__(self, E, elems, nodes, cons, rho=10**8):
        super(EulerBucklingKS, self).__init__()

        self.E = E
        self.elems = elems
        self.nodes = nodes
        self.cons = cons
        self.rho = rho

        self.n = elems.shape[0]

        self.add_param('areas', val=numpy.zeros(self.n))
        self.add_param('stress', val=numpy.zeros(self.n))
        self.add_output('neg_stress_plus_buckling_con',val=0.0)

        self.elem_len = getLength(nodes,elems) # get the length of bars
        k_list =  getkFactor(elems,cons)

        self.eq_len = numpy.zeros(self.n)
        for i in xrange(self.n):
            self.eq_len[i] = self.elem_len[i] * k_list[i]

    def solve_nonlinear(self, params, unknowns, resids):

        rho = self.rho

        area = params['areas']
        stress = params['stress']

        buckling_con = numpy.zeros(self.n)
        for i in xrange(self.n):
            buckling_con[i] = -stress[i]-numpy.pi*self.E*area[i]/(self.eq_len[i]**2)

        buckling_con_max = numpy.max(buckling_con)

        unknowns['neg_stress_plus_buckling_con'] = buckling_con_max + 1 / rho * \
                                numpy.log(numpy.sum(numpy.exp(rho*(buckling_con - buckling_con_max))))

    def linearize(self, params, unknowns, resids):

        n = self.n
        rho = self.rho
        E = self.E
        eq_len = self.eq_len

        area = params['areas']
        stress = params['stress']

        # get the maximum con index
        f = numpy.zeros(n)
        for i in xrange(n):
            f[i] = -stress[i]-numpy.pi*E/eq_len[i]**2*area[i]

        indmax = numpy.argmax(f)

        # get the normalized weight
        weight_array = numpy.zeros(n)
        for i in xrange(n):
            weight_array[i] = numpy.exp(rho*(f[i]-f[indmax]))
        weight_sum = numpy.sum(weight_array)
        for i in xrange(n):
            # normalize the weight
            weight_array[i] /= weight_sum

        # get the derivative wrt area for KS function
        dKS_darea = numpy.zeros(n)
        for i in xrange(n):
            dfi_darea = numpy.zeros(n)
            dfi_darea[i] = -numpy.pi*E/eq_len[i]**2

            dKS_darea += weight_array[i]*dfi_darea

        # get the derivative wrt stress for KS function
        dKS_dstress = numpy.zeros(n)
        for i in xrange(n):
            dfi_dstress = numpy.zeros(n)
            dfi_dstress[i] = -1.0

            dKS_dstress += weight_array[i]*dfi_dstress


        jacs = {}
        jacs['neg_stress_plus_buckling_con', 'areas'] = dKS_darea.reshape((1, self.n))
        jacs['neg_stress_plus_buckling_con', 'stress'] = dKS_dstress.reshape((1, self.n))
        return jacs
