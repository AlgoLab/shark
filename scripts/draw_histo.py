import sys

import matplotlib.pyplot as plt
import numpy as np

def main():
    data = []
    for n in open(sys.argv[1]).readlines():
        n = int(n.strip("\n"))
        data.append(n)

    print("Minimum coverage:", min(data))
    plt.hist(data, bins=50, color='blue', alpha=0.75)
    plt.xlabel('Coverage')
    plt.ylabel('# Transcripts')
    plt.show()

if __name__ == '__main__':
    main()
