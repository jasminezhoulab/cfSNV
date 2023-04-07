# https://github.com/llimllib/jenks-python

import json
from pprint import pprint as pp

inf = float('inf')

def jenks_init_matrices(data, n_classes):
    #where N = len(data) and K = len(n_classes)+1
    k = n_classes + 1
    #An N x K matrix of optimal lower class limits
    lower_class_limits =    [[0.] * k for i in range(len(data))]
    #An N x K matrix of optimal variance combinations
    #the first row is 0s, the rest of the rows are [0, inf*K]
    variance_combinations = [[0.] * k] + \
                            [[0.] + [inf] * k
                             for i in range(len(data))]
    #the original code sets lc[0][i] to 1... figure out what this row is actually for!
    for i in range(1, n_classes+1):
        lower_class_limits[1][i] = 1.
    return lower_class_limits, variance_combinations




def jenks_matrices(data, n_classes):
    lower_class_limits, variance_combinations = jenks_init_matrices(data, n_classes)
    variance = 0.0
    for l in range(1, len(data)):
        sum = 0.0
        sum_squares = 0.0
        w = 0.0
        for m in range(l):
            # `III` originally
            lower_class_limit = l - m + 1
            val = data[lower_class_limit-1]
            # here we're estimating variance for each potential classing
            # of the data, for each potential number of classes. `w`
            # is the number of data points considered so far.
            w += 1
            # increase the current sum and sum-of-squares
            sum += val
            sum_squares += val * val
            # the variance at this point in the sequence is the difference
            # between the sum of squares and the total x 2, over the number
            # of samples.
            variance = sum_squares - (sum * sum) / w
            i4 = lower_class_limit - 1
            if i4 != 0:
                for j in range(2, n_classes+1):
                    if variance_combinations[l][j] >= (variance + variance_combinations[i4][j - 1]):
                        lower_class_limits[l][j] = lower_class_limit
                        variance_combinations[l][j] = variance + variance_combinations[i4][j - 1]
        lower_class_limits[l][1] = 1.
        variance_combinations[l][1] = variance
    return lower_class_limits, variance_combinations

def get_jenks_breaks(data, lower_class_limits, n_classes):
    k = len(data) - 1
    kclass = [0.] * (n_classes+1)
    countNum = n_classes
    kclass[0] = data[0]
    kclass[n_classes] = data[len(data) - 1]
    while countNum > 1:
        elt = int(lower_class_limits[k][countNum] - 2)
        kclass[countNum - 1] = data[elt]
        k = int(lower_class_limits[k][countNum] - 1)
        countNum -= 1
    return kclass

def separate_data(data, boundaries):
    clusters = []
    cluster = []
    boundaries.pop(0)
    i=0
    for d in data:
        if d <= boundaries[i]:
           cluster.append(d)
        else:
            clusters.append(cluster)
            cluster = [d]
            i += 1
    clusters.append(cluster)
    return clusters

def jenks(data, n_classes):
    if n_classes > len(data): return ###nothing?
    data.sort()
    lower_class_limits, _ = jenks_matrices(data, n_classes)
    boundaries = get_jenks_breaks(data, lower_class_limits, n_classes)
    return separate_data(data, boundaries), boundaries


def main():
    pp(jenks([1,2,3,4,5,6,7,8,9,10], 3))






