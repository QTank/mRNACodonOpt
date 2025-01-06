import matplotlib.pyplot as plt
import numpy as np
import time, util


def plot(opt_data_name, exact_data_name, one_hot_data_name, n, output_file):
	opt_data = util.read_data(opt_data_name)
	exact_data = util.read_data(exact_data_name)
	one_hot_data = util.read_data(one_hot_data_name)
	x, y, one_hot = [], [], []

	for sequence in opt_data:
		x.append(float(exact_data[sequence]['value']))
		y.append(float(opt_data[sequence]['value']))
		one_hot.append(float(one_hot_data[sequence]['value']))

	plt.scatter(x, y, marker="*", label="optimal value on VQE")
	# plt.scatter(x, one_hot, marker="h")
	max_value = max(one_hot)
	max_value = round(max_value)
	plt.plot(range(max_value), range(max_value), color='r', label="exact value")
	plt.xlabel("Exact minimal value")
	plt.ylabel("Optimal minimal value")
	plt.title("Fragment with length of %i" % n)
	plt.legend()
	#plt.show()
	plt.savefig(output_file)


def plot_dense(opt_data_name, exact_data_name, n, output_file):
	opt_data = util.read_data(opt_data_name)
	exact_data = util.read_data(exact_data_name)
	x, y, one_hot = [], [], []

	for sequence in opt_data:
		x.append(float(exact_data[sequence]['value']))
		y.append(float(opt_data[sequence]['value']))

	plt.scatter(x, y, marker="*", label="optimal value on VQE")
	# plt.scatter(x, one_hot, marker="h")
	max_value = max(x)
	max_value = round(max_value)
	plt.plot(range(max_value), range(max_value), color='r', label="exact value")
	plt.xlabel("Exact minimal value")
	plt.ylabel("Optimal minimal value")
	plt.title("Each amino acid sequence with length of %i" % n)
	plt.legend()
	# plt.show()
	plt.savefig(output_file)


split_len = 8
path = f"../output/{split_len}"
dense = "dense"
exact = "exact"
one_hot = "one_hot"
time_stamp = time.strftime("%m-%d_%H-%M-%S")

output_file_name = f"{path}/P0DTC2_{time_stamp}_presentation"
# opt_file_name = f"{path}/P0DTC2_04-19_12-04-21.csv"
dense_file_name = f"{path}/{dense}/P0DTC2_04-16_18-12-59.csv"
# exact_file_name = f"{path}/{exact}/P0DTC2_04-19_13-13-13.csv"
exact_file_name = f"{path}/{exact}/P0DTC2_04-18_16-04-03.csv"
# one_hot_file_name = f"{path}/{one_hot}/P0DTC2_04-19_10-06-14.csv"
# plot(opt_file_name, exact_file_name, one_hot_file_name, split_len, output_file_name)

plot_dense(dense_file_name, exact_file_name, split_len, output_file_name)



