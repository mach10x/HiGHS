import highspy

h = highspy.Highs()
x = {}
nvars=10
for i in range(nvars):
    x[i] = h.addVar(obj=1, name=f"x_{i}")

for i in range(nvars - 1):
    h.addConstr(x[i] + x[i + 1] >= 1, name=f"c_{i}")

#h.minimize()
h.run()

#h.writeModel("")
