#!/usr/bin/python3
'''
Author:
Dr. Chuanxun Su
Key Laboratory of Quantum Information, University of Science and Technology of China, Hefei 230026, China

Email:
sucx@ustc.edu.cn / suchuanxun@163.cn
'''
from spap import run_spap
import datetime


start_time = datetime.datetime.now()
print('Start time: ' + datetime.datetime.strftime(start_time, '%Y-%m-%d %H:%M:%S') + '\n')

run_spap(
    l_cif=True,
    # l_db=True,
    r_cut_off=9.0,
    symprec=0.01,
    e_range=666,
    # l_poscar=True,
    # threshold=0.075,
    # l_comp=False,
    ilat=2,
    iag=2,
    work_dir=r'D:\data\Examples\6_example\results',
    i_mode=1,
    # lplt=True,
    # nrdcmf=True,
    # nsm=True,
    # nfu=True,
    # lprm=True,
)

end_time = datetime.datetime.now()
print(('\nEnd time: ' + datetime.datetime.strftime(end_time, '%Y-%m-%d %H:%M:%S')))
print('Duration: ' + str(end_time - start_time))
timelog = open('timelog.dat', 'a')
timelog.write('Duration: ' + str(end_time - start_time) + '\n')
timelog.close()