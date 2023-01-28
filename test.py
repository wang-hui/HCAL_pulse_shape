#0.018729, 0.660921, 0.156969, 0.0469796, 0.0213707

import numpy as np

shapes = [
    np.array([0.02, 0.60, 0.20, 0.10, 0.05, 0.03]),
    np.array([0.01, 0.30, 0.19, 0.01, 0.30, 0.19]),
    np.array([0.00, 0.01, 0.09, 0.60, 0.20, 0.10]),
]

for shape in shapes:
    print "=========", sum(shape) , "========="
    
    filter_ascend = np.array([-1, 1])
    filter_descend = np.array([1, -1])
    
    fMap_ascend = np.zeros(shape.size - 1)
    fMap_descend = np.zeros(shape.size - 1)
    
    for i in range(shape.size - 1):
        temp = shape[i:i+2]
        fMap_ascend[i] = np.dot(temp, filter_ascend)
        fMap_descend[i] = np.dot(temp, filter_descend)
    
    print fMap_ascend
    print fMap_descend
    
    filter3 = np.array([1, -1, -1, -1, -1])
    filter4 = np.array([-1, 1, 0.4, 0.2, 0.1])
    
    print np.dot(fMap_ascend, filter3)
    print np.dot(fMap_descend, filter4)
