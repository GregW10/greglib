import numpy as np


def main():
    arr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 8, 7, 6, 5, 4, 3, 2, 1])
    arr = arr.reshape((2, 2, 2, 2))
    print(arr)
    arr = np.array([])
    print(arr.shape)


if __name__ == "__main__":
    main()
