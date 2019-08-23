
"""Simple matshow() example."""
import matplotlib.pyplot as plt
import numpy as np





input = np.loadtxt("ZMatrixNoReorder.txt", delimiter=' ')

# Display matrix
plt.matshow(input)

plt.show()