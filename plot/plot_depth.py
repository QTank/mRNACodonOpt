import matplotlib.pyplot as plt
import time, csv


def plot_gates(max_depth, mean_depth, output_file):
    x = range(6, 20)
    plt.plot(x, max_depth[0], marker='o', label='maximum on one-hot')
    plt.plot(x, max_depth[1], marker='*', label='maximum on dense')
    plt.plot(x, mean_depth[0], marker='h', label='mean on one-hot')
    plt.plot(x, mean_depth[1], marker='v', label='mean on dense')
    plt.xlabel("The length of each fragment of protein")
    plt.ylabel("Depth of Circuits")
    plt.title("Depth of circuits with one-hot and dense encoding")
    plt.legend()
    # plt.show()
    plt.savefig(output_file)


def read_depth(file_name):
    mean_depth = []
    max_depth = []

    with open(file_name, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar="|")
        for row in reader:
            mean_depth.append(float(row['Mean_Depth']))
            max_depth.append(float(row['Max_Depth']))

    return mean_depth, max_depth

encode_type = "dense"
dense_mean_depth, dense_max_depth = read_depth(f"../output/{encode_type}_depth_gate.csv")
encode_type = "one_hot"
one_hot_mean_depth, one_hot_max_depth = read_depth(f"../output/{encode_type}_depth_gate.csv")
plot_gates([one_hot_max_depth, dense_max_depth], [one_hot_mean_depth, dense_mean_depth], f"../output/figure_depth.png")


