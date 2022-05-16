def generate_homogenous_symmetric_polinomials_powers(n_vars, max_pow = None):

    max_pow = max_pow if max_pow else n_vars
    max_num = 10 ** n_vars
    combinations = []

    for num in range(round(max_num)):
        algs = [int(x) for x in str(num)]
        summ = sum(algs)
        if summ <= max_pow and summ > 0:
            combinations.append((algs, num, summ))

    combinations.sort(key = lambda x: (x[2], -x[1]))
    return [ [0 for x in range(n_vars - len(combination[0]))] + combination[0] for combination in combinations]