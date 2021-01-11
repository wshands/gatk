import argparse

from matplotlib import pyplot as plt
import pandas as pd


def plot_depth(args):
    mt = pd.read_csv(args.depth_file, sep='\t', index_col=[0, 1, 2], header=0)
    sample_depth = pd.read_csv(args.mean_depth_file, sep='\t', index_col=0, header=None).transpose()

    mt = mt.div(sample_depth.iloc[0])

    out_path = args.out_name + ".depth.png"
    figure_width = 6
    figure_height = 4

    fig = plt.figure(figsize=(figure_width, figure_height))

    linewidth = 0.5

    mt[mt.keys()[1:]].plot(legend=False, rot=90, color='r', marker='.', linewidth=0, markersize=1, alpha=0.05)

    mt['HG00096'].plot(legend=False, rot=90, color='b', marker='.', linewidth=linewidth, markersize=1)
    mt[mt.keys()[1:]].median(axis=1).plot(legend=False, rot=90, color='k', marker='.', linewidth=linewidth, markersize=1)
    plt.ylim([0, 2])

    plt.savefig(out_path)
    plt.close(fig)


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--depth-file', help='Depth file path', required=True)
    parser.add_argument('--mean-depth-file', help='Mean sample depth file path', required=True)
    parser.add_argument('--out-name', help='Output base name', required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    plot_depth(args)


if __name__ == "__main__":
    main()
