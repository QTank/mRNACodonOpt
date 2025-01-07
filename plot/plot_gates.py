import matplotlib.pyplot as plt
import time, csv


def plot_gates(max_depth, mean_depth, output_file):
    x = range(6, 20)
    plt.plot(x, max_depth[0], marker='o', label='maximum on one-hot')
    plt.plot(x, max_depth[1], marker='*', label='maximum on dense')
    plt.plot(x, mean_depth[0], marker='h', label='mean on one-hot')
    plt.plot(x, mean_depth[1], marker='v', label='mean on dense')
    plt.xlabel("Length of each fragment")
    plt.ylabel("Gates of circuits")
    # plt.title("Total Gates of circuits with one-hot and dense encoding")
    plt.legend()
    # plt.show()
    plt.savefig(output_file)


def read_gates(file_name):
    mean_gate = []
    max_gate = []

    with open(file_name, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar="|")
        for row in reader:
            mean_gate.append(float(row['Mean_Gate']))
            max_gate.append(float(row['Max_Gate']))

    return mean_gate, max_gate

encode_type = "dense"
dense_mean_gate, dense_max_gate = read_gates(f"../output/{encode_type}_depth_gate.csv")
encode_type = "one_hot"
one_hot_mean_gate, one_hot_max_gate = read_gates(f"../output/{encode_type}_depth_gate.csv")
plot_gates([one_hot_max_gate, dense_max_gate], [one_hot_mean_gate, dense_mean_gate], f"../output/figure_gates.png")


