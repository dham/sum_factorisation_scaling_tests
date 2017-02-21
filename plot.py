import pylab as plt
import json
import numpy as np

dataset = json.load(file("data.json", "r"))

scaling = {"matrix": lambda dim: 2 * dim + 1,
           "action": lambda dim: dim + 1}


for dim in (2, 3):
    for operation, data in dataset[str(dim)].iteritems():

        plt.figure()

        for form, times in data.iteritems():

            degrees = sorted(times.keys())

            t = [min(times[d]) for d in degrees]

            degrees = np.array(degrees, dtype=float)

            plt.loglog(degrees, t, "*")

            reference = 8
            modeltimes = (degrees / reference) ** scaling[operation](dim) * t[reference]

        plt.loglog(degrees, modeltimes)
        plt.title("%s in %dD" % (operation, dim))
        plt.legend(data.keys() + ["O(p**%d)" % scaling[operation](dim)], loc="lower right")
        plt.xlabel("Polynomial degree")
        plt.ylabel("Time per cell (s)")
    
plt.show()

