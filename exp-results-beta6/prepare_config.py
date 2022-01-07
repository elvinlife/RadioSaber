import argparse
import sys

def main():
    pass


if __name__ == "__main__":
    args = sys.argv[1:]
    print(args)
    num_ues = args[0]
    num_slices = args[1]
    f = open("config", "w")
    f.write(str(0)+" "+str(num_slices))
    f.write("\n")
    slice_weight = 1/int(num_slices)
    weights = [slice_weight for _ in range(int(num_slices))]
    print(weights)
    f.write(" ".join(map(str, weights)))
    f.write("\n")
    ues = [num_ues for _ in range(int(num_slices))]
    f.write(" ".join(map(str, ues)))
    f.write("\n")