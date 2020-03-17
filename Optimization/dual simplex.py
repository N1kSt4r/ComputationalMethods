import numpy as np
import itertools

AT = np.array([
    [8, 25, 1, 0, 0],
    [8, 5, 0, 1, 0],
    [1, 5, 0, 0, 1]
]).T
b = [800, 640, 145]
c = np.array([80, 70, 0, 0, 0])
limitations = [[3, 78], [0, 28], [0, 776], [0, 616], [0, 142]]
rows, columns = len(AT[0]), len(AT)

def dual_task(AT, b, c, limitations):
    new_AT = np.vstack((AT.T, -np.identity(columns), np.identity(columns)))
    new_limitations = np.array(limitations)
    new_c = -np.append(np.append(b, -new_limitations.T[0]), new_limitations.T[1])
    new_b = [*np.append(c, np.zeros(columns - len(c)))]
    new_limitations = [[0, np.inf] for _ in range(columns * 3)]
    print(new_AT.T)
    print(new_b)
    print(new_c)
    print(new_limitations)
    return new_AT, new_b, new_c, new_limitations

base = None

for base_columns in itertools.combinations(range(columns), rows):
    if abs(np.linalg.det(AT[np.array(base_columns)])) > 1e-6:
        base = {*base_columns}
        break

print(' A:\n{}'.format(AT.T))
print(' base: {}'.format(base))
print(' c: {}'.format(c))
print(' limitation: {}\n'.format(limitations))

for i in range(1, 10):
    mask = np.zeros(columns, dtype=bool)
    mask[list(base)] = True
    J_b = np.arange(columns)[mask]
    J_n = np.arange(columns)[np.logical_not(mask)]
    AT_b = AT[J_b]
    AT_n = AT[J_n]
    c_b = c[J_b]
    u = np.linalg.solve(AT_b, c_b)
    delta_n = [c[j] - np.dot(u, AT[j]) for j in J_n]

    print('{} iteration'.format(i))
    print('\t A_base:')
    for string in AT_b.T:
        print('\t\t {}'.format(string))
    print('\t c_base: {}'.format(c_b))
    print('\t index_base: {}'.format(J_b))
    print('\t index_nobase: {}'.format(J_n))
    print('\t potential: {}'.format(u))
    print('\t delta_nobase: {}'.format(delta_n))

    x_n = [limitations[J_n[j]][1] if delta_n[j] > 0 else limitations[J_n[j]][0] for j in range(len(delta_n))]
    x_b = np.linalg.solve(AT_b.T, b - np.dot(AT_n.T, x_n))

    print('\t x_nobase: {}'.format(x_n))
    print(f'{AT_b.T} * x = {b} - {np.dot(AT_n.T, x_n)}')
    print('\t x_base: {}'.format(x_b))

    j_star = None
    for index, x in enumerate(x_b):
        if not limitations[J_b[index]][0] < x < limitations[J_b[index]][1]:
            j_star = index

    if j_star is None:
        print('\n optimality criterion is satisfied')
        x = np.zeros(columns)
        x[J_b] = x_b
        x[J_n] = x_n
        print(f' optimal plan: {x}')
        print(f' optimal target function: {np.dot(x, c)}')
        break

    p_u = np.zeros(len(J_b))
    if abs(x_b[j_star] - limitations[J_b[j_star]][0]) < abs(x_b[j_star] - limitations[J_b[j_star]][1]):
        p_u[j_star] = -np.sign(x_b[j_star] - limitations[J_b[j_star]][0])
    else:
        p_u[j_star] = -np.sign(x_b[j_star] - limitations[J_b[j_star]][1])
    p_u = np.linalg.solve(AT_b, p_u)

    p_n = [-np.dot(AT[j], p_u) for j in J_n]
    print(f'\t p_u: {p_u}')
    print(f'\t p_nobase: {p_n}')

    sigma = [-delta_n[j] / p_n[j] if delta_n[j] * p_n[j] < 0 else np.inf for j in range(len(J_n))]
    print(f'\t sigma: {sigma}')
    j_0 = J_n[np.argmin(sigma)]

    print(f'\t j*, j0: {j_star}, {j_0}')
    base -= {J_b[j_star]}
    base.add(j_0)
    print(f'\t new base:', *base, end='\n\n')
