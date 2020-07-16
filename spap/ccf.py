'''
Author:
Dr. Chuanxun Su

Email:
suchuanxun@163.cn / scx@calypso.cn
'''
import numpy as np
from ase import Atoms
from math import sqrt, pi, exp
from ase.data import atomic_numbers
# from numba import jit


def struc2ccf(struc, r_cut_off, r_vector, apw=60.0, ftype='CCF'):
    rho = len(struc.numbers) / struc.get_volume()
    ccf = d2ccf(cal_inter_atomic_d(struc, r_cut_off), r_cut_off, r_vector, apw, ftype, rho)
    nspec = get_nspec(struc)
    lccf=len(r_vector)
    for i in range(1, nspec + 1):
        for j in range(i, nspec + 1):
            pairt = str(i) + '_' + str(j)
            if not pairt in ccf:
                ccf[pairt] = np.zeros(lccf,float)
    return ccf


def cal_ccf_d(ccf1, ccf2):
    ccf_wf1 = {}
    lenth = len(ccf1)
    for key in ccf1.keys():
        ccf_wf1[key] = np.sum(ccf1[key])
    sum_ccf_wf = sum(ccf_wf1.values())
    if abs(sum_ccf_wf) < 1.0e-20:
        for key in ccf_wf1.keys():
            ccf_wf1[key] = 1.0 / lenth
    else:
        for key in ccf_wf1.keys():
            ccf_wf1[key] = ccf_wf1[key] / sum_ccf_wf
    ccf_wf2 = {}
    for key in ccf2.keys():
        ccf_wf2[key] = np.sum(ccf2[key])
    sum_ccf_wf = sum(ccf_wf2.values())
    if abs(sum_ccf_wf) < 1.0e-20:
        for key in ccf_wf2.keys():
            ccf_wf2[key] = 1.0 / lenth
    else:
        for key in ccf_wf2.keys():
            ccf_wf2[key] = ccf_wf2[key] / sum_ccf_wf
    for key in ccf_wf1.keys():
        ccf_wf1[key] = (ccf_wf1[key] + ccf_wf2[key]) / 2.0
    struc_d = 0.0
    for key in ccf1.keys():
        struc_d += ccf_wf1[key] * pearson_cc(ccf1[key], ccf2[key])
    return 1 - struc_d


# @jit(nopython=True)
def cal_inter_atomic_d(struc, r_cut_off):
    distances = {}
    # square_rcut=r_cut_off*r_cut_off
    i_range = cell_range(struc.cell, struc.pbc, r_cut_off)
    natoms = len(struc.numbers)
    ele_tag = element_tag(struc.numbers)
    # square_delta=0.01
    d_delta = 0.001
    # d_delta2=0.001
    # n_atom_dict=count_atoms_dict(struc.numbers)
    pair_tags = []
    l_d_empty = int(r_cut_off / d_delta + 0.5) + 1
    # d_empty=np.zeros(l_d_empty,dtype=int)
    # d_empty = l_d_empty * [0]

    n_species = len(ele_tag)
    for i1 in range(n_species):
        for i2 in range(n_species):
            if i1 < i2:
                pair_tags.append(str(i1 + 1) + '_' + str(i2 + 1))
            else:
                pair_tags.append(str(i2 + 1) + '_' + str(i1 + 1))
                # distances[pair_tags[-1]] = d_empty.copy()
                # distances[pair_tags[-1]] = d_empty[:]
                distances[pair_tags[-1]] = l_d_empty * [0]
    pair_tags = [[x for x in pair_tags[i * n_species:(i + 1) * n_species]] for i in range(n_species)]
    transa = -i_range[0] * struc.cell[0]
    for ia in range(-i_range[0] + 1, i_range[0] + 1):
        transa = np.row_stack((transa, ia * struc.cell[0]))
    transb = -i_range[1] * struc.cell[1]
    for ib in range(-i_range[1] + 1, i_range[1] + 1):
        transb = np.row_stack((transb, ib * struc.cell[1]))
    transc = -i_range[2] * struc.cell[2]
    for ic in range(-i_range[2] + 1, i_range[2] + 1):
        transc = np.row_stack((transc, ic * struc.cell[2]))
    # temp_d = d_empty.copy()
    # temp_d = d_empty[:]
    temp_d = l_d_empty * [0]
    for ia in range(-i_range[0], i_range[0] + 1):
        for ib in range(-i_range[1], i_range[1] + 1):
            for ic in range(-i_range[2], i_range[2] + 1):
                for i1 in range(natoms):
                    for i2 in range(i1 + 1, natoms):
                        d = sqrt(np.sum(
                            np.square(struc.positions[i1] - struc.positions[i2] + transc[ic + i_range[2]] +
                                      transb[ib + i_range[1]] + transa[ia + i_range[0]])))
                        if d < r_cut_off:
                            distances[pair_tags[ele_tag[struc.numbers[i1]] - 1][ele_tag[struc.numbers[i2]] - 1]][
                                -int((r_cut_off - d) / d_delta + 0.5) - 1] += 1
            for ic in range(1, i_range[2] + 1):
                d = sqrt(np.sum(np.square(transc[ic + i_range[2]] + transb[ib + i_range[1]] + transa[ia + i_range[0]])))
                if d < r_cut_off:
                    temp_d[-int((r_cut_off - d) / d_delta + 0.5) - 1] += 1
        for ib in range(1, i_range[1] + 1):
            d = sqrt(np.sum(np.square(transb[ib + i_range[1]] + transa[ia + i_range[0]])))
            if d < r_cut_off:
                temp_d[-int((r_cut_off - d) / d_delta + 0.5) - 1] += 1
    for ia in range(1, i_range[0] + 1):
        d = sqrt(np.sum(np.square(transa[ia + i_range[0]])))
        if d < r_cut_off:
            temp_d[-int((r_cut_off - d) / d_delta + 0.5) - 1] += 1
    for ele_n, n in count_atoms_dict(struc.numbers).items():
        for i, n_pair in enumerate(temp_d):
            if n_pair != 0:
                distances[pair_tags[ele_tag[ele_n] - 1][ele_tag[ele_n] - 1]][i] += n_pair * n
    final_d = {}
    for key1 in distances.keys():
        final_d[key1] = {}
        for i, n_pair in enumerate(distances[key1]):
            if n_pair != 0:
                final_d[key1][r_cut_off - (l_d_empty - i - 1) * d_delta] = float(n_pair) / natoms
    return final_d


# @jit()
def d2ccf(distances, r_cut_off, r_vector, a=60.0, ftype='CCF', rho=1.0):
    ccf = {}
    if ftype == 'RDF':
        norm = 2.0 * pi * rho
    for key1 in distances.keys():
        for key2 in distances[key1].keys():
            if key1 in ccf:
                if ftype == 'CCF':
                    ccf[key1] = ccf[key1] + gaussian_f(distances[key1][key2] * weight_f(key2, r_cut_off), key2,
                                                       r_vector, a)
                elif ftype == 'RDF':
                    ccf[key1] = ccf[key1] + gaussian_f(distances[key1][key2] / norm / key2 ** 2, key2, r_vector, a)
            else:
                if ftype == 'CCF':
                    ccf[key1] = gaussian_f(distances[key1][key2] * weight_f(key2, r_cut_off), key2, r_vector, a)
                elif ftype == 'RDF':
                    ccf[key1] = gaussian_f(distances[key1][key2] / norm / key2 ** 2, key2, r_vector, a)
    return ccf


def weight_f(r, r_cut_off):
    tail_r = 0.5
    if r < r_cut_off - tail_r:
        return 1.0
    else:
        return exp(-3.0 * (r - r_cut_off + tail_r) / tail_r)


# @jit(nopython=True)
def pearson_cc(x, y):
    # average1=np.mean(ccf1)
    # average2=np.mean(ccf2)
    # Send in mean(x) might be more efficient. Because we calculated sum(x) before.
    local1 = x - np.mean(x)
    local2 = y - np.mean(y)
    # sum1=np.sum(np.multiply(local1,local2))
    sum2 = np.sum(np.square(local1))
    sum3 = np.sum(np.square(local2))
    if sum2 + sum3 == 0.0:
        return 1.0
    elif (sum2 == 0.0) or (sum3 == 0.0):
        return 0.0
    return np.sum(np.multiply(local1, local2)) / sqrt(sum2 * sum3)


# @jit
def gaussian_f(weight_f, b, x, a=60.0):
    # a = 60.0
    # a=290
    # return weight_f * sqrt(a / pi) * np.exp(-a * (x - b)**2)
    return weight_f * sqrt(a / pi) * np.exp(-a * np.square(x - b))


def element_tag(numbers, irt=1):
    ele_n = {}
    ele_tag = {}
    for i in numbers:
        if i in ele_n:
            ele_n[i] += 1
        else:
            ele_n[i] = 1
    sotn = sorted(ele_n.items(), key=lambda x: x[1])
    for i, x in enumerate(sotn):
        if irt == 1:
            ele_tag[x[0]] = i + 1
        elif irt == 2:
            for symb in atomic_numbers.keys():
                if atomic_numbers[symb] == x[0]:
                    break
            ele_tag[x[0]] = [i + 1, symb]
        # Worth to try something like below
        # ele_tag[x[0]] = i
    return ele_tag


def get_nspec(struc):
    eles = {}
    for i in struc.numbers:
        if not i in eles:
            eles[i] = 1
    return len(eles)


def cell_range(cell, pbc, rcut):
    recipc_no2pi = Atoms(cell=cell).get_reciprocal_cell()
    i_range = []
    for i in range(3):
        if pbc[i] == True:
            i_range.append(int(rcut * ((np.sum(recipc_no2pi[i] ** 2)) ** 0.5) + 1.0e-6) + 1)
            # i_range.append(int(rcut * ((np.sum(recipc_no2pi[i] ** 2)) ** 0.5)) + 1)
        else:
            i_range.append(0)
    return i_range
    # return [int(rcut * ((np.sum(recipc_no2pi[i] ** 2)) ** 0.5)) + 1 for i in range(3)]


def count_atoms_dict(numbers):
    ctype = {}
    for i in numbers:
        if i in ctype:
            ctype[i] += 1
        else:
            ctype[i] = 1
    return ctype
