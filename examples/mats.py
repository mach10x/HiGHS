import highspy
import numpy as np
import time

h = highspy.Highs()
nk = 6
inf = highspy.kHighsInf
nvars = 10
for k in range(nk):
    h.clear()
    h.setOptionValue('output_flag', False)
    tt0 = time.perf_counter()
    for i in range(nvars):
        h.addVar(0, inf)
        h.changeColCost(i, 1)
    row_nnz = 2
    index = np.array([0, 1])
    value = np.array([1, 1], dtype=np.double)
    for i in range(nvars - 1):
        index[0] = i
        index[1] = i+1
        h.addRow(1, inf, row_nnz, index, value)
    tt_build = time.perf_counter() - tt0
    tt0 = time.perf_counter()
    h.run()
    tt_run = time.perf_counter() - tt0
    tt0 = time.perf_counter()
    solution = h.getSolution()
    tt_solution = time.perf_counter() - tt0
    print('N = {0:6d}; build = {1:6f}; run = {2:6f}; solution = {3:6f}'.format(nvars, tt_build, tt_run, tt_solution))
    nvars *= 10


#h.writeModel("")
