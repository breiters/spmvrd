import matplotlib.pyplot as plt
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="table")
    parser.add_argument("-f", "--file", type=str,
                        default=None, help="csv file")
    args = parser.parse_args()

    fig, ax = plt.subplots(1, 1, figsize=(16, 9))

    df = pd.read_csv(args.file)
    df['mindist'] = df['mindist'] * 256 # convert to bytes
    df = df[df['shared'] == 0]
    df = df[df['cache_id'] == 0]

    plt.axvline(x = 64 * 1024, color = 'b', label = 'L1')
    plt.axvline(x = 8 * 1024 * 1024, color = 'b', label = 'L2')

    df.plot.bar(x='mindist', y='count', rot=45)

    ax.set_xlabel('Reuse Distance [Bytes]')
    ax.set_ylabel('Reference Count')

    plt.savefig(args.file + '.png')
