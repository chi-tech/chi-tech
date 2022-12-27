import matplotlib.pyplot as plt
import numpy as np

xsfile = open("../ChiTest/xs_graphite_pure.cxs", "r")
# xsfile = open("./xs_2g_mat1up.cxs", "r")

num_groups: int = 0
num_moments: int = 0
T = np.zeros((num_groups, num_groups))
S = np.zeros((num_groups, num_groups))
S1 = np.zeros((num_groups, num_groups))

lines = xsfile.readlines()
line_count = 0
for line in lines:
    words = line.split()
    if len(words) == 0:
        continue

    if words[0] == "NUM_GROUPS":
        num_groups = int(words[1])
        T = np.zeros((num_groups, num_groups))
        S = np.zeros((num_groups, num_groups))
        S1 = np.zeros((num_groups, num_groups))

    if words[0] == "SIGMA_T_BEGIN":
        k = line_count+1
        while True:
            block_words = lines[k].split()
            if len(block_words) >= 1:
                if block_words[0] == "SIGMA_T_END":
                    break
            if len(block_words) > 1:
                if block_words[0].find("//") < 0:
                    g = int(block_words[0])
                    val = float(block_words[1])
                    T[g, g] = val
            k += 1

    if words[0] == "M_GPRIME_G_VAL":
        m = int(words[1])
        gp = int(words[2])
        g = int(words[3])
        val = float(words[4])

        if m == 0:
            S[g, gp] = val
        if m == 1:
            S1[g, gp] = val

    line_count += 1

xsfile.close()


def PowerIteration(A):
    y = np.ones(num_groups)
    eigv = 0.0

    Ay = np.matmul(A, y)
    eigv = np.dot(y, Ay)
    y = Ay/np.sum(Ay)
    if eigv < 0.0:
        y = y*-1.0

    converged = False
    it_counter = 0
    while not converged and it_counter < 1000:
        Ay = np.matmul(A, y)
        eigv = np.dot(y, Ay)
        y = Ay/np.sum(Ay)
        if eigv < 0.0:
            y = y*-1.0

        it_counter += 1

    print("Check: ", np.linalg.norm(np.matmul(A,y)-eigv*y))

    y = y/np.sum(y)

    return eigv, y

# ===================================== Do scattering stuff
# A = np.copy(T)
# B = np.zeros((num_groups, num_groups))
# # for g in range(0, num_groups):
# #     A[g, g] = T[g, g] - S[g, g]
# #     for gp in range(0, num_groups):
# #         if gp != g:
# #             B[g, gp] = S[g, gp]
# for g in range(63, num_groups):
#     A[g, g] = T[g, g]
#     for gp in range(0, num_groups):
#         B[g, gp] = S[g, gp]
#
# Ainv = np.linalg.inv(A)
# C = np.matmul(Ainv, B)
#
# eigv, phi = PowerIteration(C)
# print("EigenValue", eigv)
# for g in range(0, num_groups):
#     print(g, phi[g])
#
# print(np.sum(phi))

# ===================================== Plot transfer matrix
maxS = np.max(S)
Sclipped = np.clip(S, 1.0e-9, maxS)

Atest = np.log10(Sclipped) + 10.0

plt.figure()
plt.title('Transfer matrix')
plt.imshow(Atest)
plt.show()

# =====================================


A = T-S;
q = np.zeros(num_groups)
q[0] =1;
spectrum = np.linalg.solve(A,q)
print(spectrum)



