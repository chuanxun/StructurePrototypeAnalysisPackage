#!/usr/bin/python3
'''
Author:
Dr. Chuanxun Su
State Key Lab of Superhard Materials, Jilin University, Changchun, China

Email:
suchuanxun@163.cn / scx@calypso.cn
'''
from spap import run_spap
from ase.io import read

# strul=read('XDATCAR',index=':')
ccf=run_spap(
    r_cut_off=9.0,
    #ilat=1,
    # structure_list=strul[0:299],
    ftype='CCF',
    readf='8_157.cif',
    # readf='XDATCAR6',
    # readf='CONTCAR_3',
    i_mode=5,
    apw=60.0,
    # index='6933:6953:',
    # ccf_step=0.04,
    lplt=True,
)