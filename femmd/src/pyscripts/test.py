import os
import femmd as md
import femmd_module as fm
import math
import numpy as np
import matplotlib.pyplot as plt
import time
import mdmesh
import Plotter
import pdb


import timeit


def test():
    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))
    L = np.array([1.0,2.0])
    md.system_init(L)
    p=md.system_get_walls_pos()
    print(p[0][0], p[1][1])

test()