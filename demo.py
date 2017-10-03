from laplace import LaplaceTransform 
from math import sin

# the following four functions were taken from a Laplace Transform table
def sin_lt(s):
    return 1.0 / (s ** 2 + 1.0 ** 2) + 0.0j


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    order = [2, 5, 10, 20, 50, 100]
    domain = [0.01 * s for s in range(1000)]
    plt.plot(domain, [sin_lt(s) for s in domain])
    for i in order:
        func = LaplaceTransform(sin, order = i)
        plt.plot(domain, [func(s) for s in domain])
    
    plt.legend(["Exact"] + ["order = " + str(i) for i in order])
    plt.xlabel("s")
    plt.ylabel("F(s)")
    plt.ylim((-1.0, 1.0))
        
    plt.show()

    

