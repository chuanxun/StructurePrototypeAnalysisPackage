#!/usr/bin/python3
'''
Author:
Dr. Chuanxun Su
State Key Lab of Superhard Materials, Jilin University, Changchun, China

Email:
suchuanxun@163.cn / scx@calypso.cn
'''
from spap import analyze


analyze(
    # total_struc=150,
    l_cif=True,
    # r_cut_off=5.5,
    # symprec=0.02,
    # l_poscar=True,
    # threshold=0.04,
    # l_comp=False,
    # ilat=1,
    # work_dir='C:\\Users\\null\\Documents\\share\\wks\\Examples\\1_example\\results',
    work_dir='./example/results',
    # i_mode=2,
)