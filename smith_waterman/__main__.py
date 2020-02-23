from smith_waterman import algs
import numpy as np
import seaborn as sns
import pandas as pd

def read_pairs(filename):
	pairs = []
	with open(filename, 'r') as pairs_file:
		for line in pairs_file.read().splitlines():
			pairs.append(tuple(map(lambda x: read_fasta(x), line.split())))
	return pairs

def read_fasta(filename):
	seq = ""
	for line in open(filename, 'r').read().splitlines():
		if line[0] is not '>':
			seq += line
	return seq

def calc_all_scores(pairs, matrix, gap_start, gap_extend):
	scores = []
	#calculate true and false scores with specified open/extend
	for pair in pairs:
		score = algs.score(*(pair), score_matrix = matrix, 
			gap_start = gap_start, gap_extend = gap_extend)
		scores.append((pair, score))
	return scores

def calc_all_aligns(pairs, matrix, gap_start, gap_extend, filename):
	scores = []
	#calculate true and false scores with specified open/extend
	for pair in pairs:
		score = algs.score(*(pair), score_matrix = matrix, 
			gap_start = gap_start, gap_extend = gap_extend)
		scores.append((pair, score))
	#with open(filename, 'w') as f:

	return scores

def best_fp_rate(pos, neg, matrix = "BLOSUM50"):
	fp_rates = pd.DataFrame(index=range(1,21), columns=range(1,6), dtype=np.double)

	for gap_start in range(1,21):
		for gap_extend in range(1,6):
			true_scores = []

			#calculate true and false scores with specified open/extend
			true_scores = calc_all_scores(pos, matrix, 
					gap_start = gap_start, gap_extend = gap_extend)
			
			#get cutoff st true positive rate is 0.7
			true_scores = sorted([x[1] for x in true_scores])
			cutoff = true_scores[int(len(true_scores)*0.3)]

			false_scores =  calc_all_scores(neg, matrix, 
					gap_start = gap_start, gap_extend = gap_extend)

			fps = 0
			for pair, score in false_scores:
				if score > cutoff:
					fps += 1

			fp_rates[gap_extend][gap_start] = fps/len(neg)

	ax = sns.heatmap(fp_rates, annot=True, annot_kws={"size": 7})
	ax.set(xlabel='Gap extend', ylabel='Gap start')
	fig = ax.get_figure()
	fig.tight_layout()
	fig.savefig("fp_rates.png", dpi=200)

	print(fp_rates)

def compare_matrices(pos, neg):
	TSs = []
	FSs = []
	matrices = ["BLOSUM50", "BLOSUM62", "MATIO", "PAM100", "PAM250"]
	for matrix in matrices:
		pos_scores = calc_all_scores(pos, matrix, 11, 1)
		TSs.append([x[1] for x in pos_scores])
		neg_scores = calc_all_scores(neg, matrix, 11, 1)
		FSs.append([x[1] for x in neg_scores])
	algs.roc(TSs, FSs, matrices, save = True, filename = 'all_matrices_roc.png')

def compare_normalized(pos, neg, matrix):
	TSs = []
	FSs = []
	pos_scores = calc_all_scores(pos, matrix, 11, 1)
	TSs.append([x[1] for x in pos_scores])
	TSs.append([score/min(len(pair[0]), len(pair[1])) for pair,score in pos_scores])
	
	neg_scores = calc_all_scores(neg, matrix, 11, 1)
	FSs.append([x[1] for x in neg_scores])
	FSs.append([score/min(len(pair[0]), len(pair[1])) for pair,score in neg_scores])

	algs.roc(TSs, FSs, ["Raw", "Normalized"], save = True, filename = matrix + "_normalization.png")


def loss_function(true_scores, false_scores):
	"""
	Devise an optimization algorithm to modify the values in a 
	starting score matrix such as to maximize the following objective 
	function: sum of TP rates for FP rates of 0.0, 0.1, 0.2, and 0.3. 
	The maximum value for the objective function is 4.0 (where you are 
	getting perfect separation of positive and negative pairs even at
	the lowest false positive rate).
	"""

	#just to make sure
	true_scores = sorted(true_scores)
	false_scores = sorted(false_scores)

	score = 0

	for fp_rate in [0.0, 0.1, 0.2, 0.3]:
		TPs = 0.0
		cutoff = false_scores[int(len(false_scores)*(1-fp_rate))]+1
		for score in true_scores:
			if true_scores > cutoff:
				TPs += 1
		score += TPs/len(true_scores)
	return score




pos = read_pairs("Pospairs.txt")
neg = read_pairs("Negpairs.txt")

compare_normalized(pos, neg, "PAM250")

