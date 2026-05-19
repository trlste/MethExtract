import numpy as np
import pyBigWig
import argparse

parser = argparse.ArgumentParser(description='Get breakpoints from bigwig files and compute the distance matrix')
parser.add_argument('-i','--input_file',help='Bigwig input file',required=True)
parser.add_argument('-c','--chromosome',default="2",help='Chromosome name',required=False)
parser.add_argument('-s','--start_pos',type=int,default=86780000 ,help='Start position',required=False)
parser.add_argument('-e','--end_pos',type=int,default=86870000,help='End position',required=False)
parser.add_argument('-o','--output_file',help='Output file .csv of distance matrix.',required=True)

args = parser.parse_args()

def extract_intervals(bw_file, chr, start_pos=86780000, end_pos=86870000):
    bw = pyBigWig.open(bw_file)
    chr_name = "chr" + chr
    #chroms = bw.chroms(chr_name)
    #print(f"Chromosome {chr}: {chroms}.")

    #values = bw.values(chr_name, start_pos, end_pos, numpy=True)
    intervals = bw.intervals(chr_name, start_pos, end_pos)

    bw.close()
    return intervals

def extract_breakpoints(intervals):
    breakpoints = []
    for i in range(len(intervals)-1):
        end_pos = intervals[i][1]
        breakpoints.append((end_pos))

    return breakpoints

def distance_matrix(breakpoints):
    positions = sorted(breakpoints)
    pos = np.array(positions)
    return np.abs(pos[:, None] - pos[None, :])

def main():
    intervals = extract_intervals(args.input_file, args.chromosome, args.start_pos, args.end_pos)
    breakpoints = extract_breakpoints(intervals)
    dm = distance_matrix(breakpoints)
    np.savetxt(args.output_file, dm, delimiter=",")

main()