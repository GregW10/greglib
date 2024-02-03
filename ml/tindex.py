import numpy as np
import math
import random as rand


def offset(shape: tuple, *args):
    if len(shape) != len(args):
        raise ValueError("Mismatch.")
    o = args[-1]
    for counter, index in enumerate(args[:-1], start=1):
        o += index*math.prod(shape[counter:])
    return o


def main():
    d1 = 9
    d2 = 8
    d3 = 7
    d4 = 6
    d5 = 19
    d6 = 72
    shape = (d1, d2, d3, d4, d5, d6)
    arr = np.zeros(shape=shape, dtype=int)
    counter = 0
    for i in range(d1):
        for j in range(d2):
            for k in range(d3):
                for l in range(d4):
                    for m in range(d5):
                        for n in range(d6):
                            arr[i][j][k][l][m][n] = counter
                            counter += 1
    nrand_idx = 100
    for i in range(nrand_idx):
        index = (rand.randint(a=0, b=d1 - 1),
                 rand.randint(a=0, b=d2 - 1),
                 rand.randint(a=0, b=d3 - 1),
                 rand.randint(a=0, b=d4 - 1),
                 rand.randint(a=0, b=d5 - 1),
                 rand.randint(a=0, b=d6 - 1))
        print(f"True offset: {arr[index]}")
        print(f"Calculated offset: {offset(shape, *index)}\n")
    print(f"o = {offset((9, 7, 6), 4, 3, 2)}")


if __name__ == "__main__":
    main()
