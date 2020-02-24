from smith_waterman import algs
import numpy as np
import seaborn as sns
import pandas as pd
import copy
import os

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

	#speedup - read in score matrix here if its a string
	if isinstance(matrix, str):
		matrix = algs.get_scoring_matrix(matrix)

	#calculate true and false scores with specified open/extend
	for pair in pairs:
		score = algs.score(*(pair), score_matrix = matrix, 
			gap_start = gap_start, gap_extend = gap_extend)
		scores.append((pair, score))
	return scores

def calc_all_aligns(pairs, matrix, gap_start, gap_extend, filename):
	aligned = []

	#speedup - read in score matrix here if its a string
	if isinstance(matrix, str):
		matrix = algs.get_scoring_matrix(matrix)

	#calculate true and false scores with specified open/extend
	for pair in pairs:
		results = algs.align(*(pair), score_matrix = matrix, 
			gap_start = gap_start, gap_extend = gap_extend)
		aligned.append((results[1], results[2], results[0]))
	with open(filename, 'w') as f:
		for seq1, seq2, score in aligned:
			f.write(seq1 + '\n')
			f.write(seq2 + '\n')
			f.write(str(score) + '\n')
			f.write('\n')

def calc_fp_rate(true_scores, false_scores, cutoff):
	true_scores = sorted(true_scores)
	cutoff = true_scores[int(len(true_scores)*(1-cutoff))]

	fps = 0
	for score in false_scores:
		if score > cutoff:
			fps += 1

	return fps/len(false_scores)


def best_fp_rate(pos, neg, matrix = "BLOSUM50"):
	fp_rates = pd.DataFrame(index=range(1,21), columns=range(1,6), dtype=np.double)

	for gap_start in range(1,21):
		for gap_extend in range(1,6):
			true_scores = []

			#calculate true and false scores with specified open/extend
			true_scores = calc_all_scores(pos, matrix, 
					gap_start = gap_start, gap_extend = gap_extend)
			false_scores =  calc_all_scores(neg, matrix, 
					gap_start = gap_start, gap_extend = gap_extend)
			
			fp_rates[gap_extend][gap_start] = calc_fp_rate(
				[x[1] for x in true_scores], 
				[y[1] for y in false_scores], 0.7)

	ax = sns.heatmap(fp_rates, annot=True, annot_kws={"size": 7})
	ax.set(xlabel='Gap extend', ylabel='Gap start')
	fig = ax.get_figure()
	fig.tight_layout()
	fig.savefig("fp_rates.png", dpi=200)

	print(fp_rates)

def compare_matrices(pos, neg, 
	matrices = ["BLOSUM50", "BLOSUM62", "MATIO", "PAM100", "PAM250"],
	gen_roc = True,
	filename = "all_matrices_roc.png"):

	TSs = []
	FSs = []
	for matrix in matrices:
		pos_scores = calc_all_scores(pos, matrix, 11, 1)
		TSs.append([x[1] for x in pos_scores])
		neg_scores = calc_all_scores(neg, matrix, 11, 1)
		FSs.append([x[1] for x in neg_scores])
	if gen_roc:
		algs.roc(TSs, FSs, matrices, save = True, filename = filename)
	else:
		for ts, fs, matrix in zip(TSs, FSs, matrices):
			print(matrix, calc_fp_rate(ts, fs, 0.7))

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

	fp_score = 0

	for fp_rate in [0.0, 0.1, 0.2, 0.3]:
		TPs = 0.0
		cutoff = false_scores[int(len(false_scores)*(1-fp_rate))-1]
		for score in true_scores:
			if score > cutoff:
				TPs += 1
		fp_score += (TPs/len(true_scores))
	return fp_score

def check_symmetry(matrix):
	for aa1 in matrix:
		for aa2 in matrix[aa1]:
			assert(matrix[aa1][aa2] == matrix[aa2][aa1])

def mutate_matrix(starting_matrix, sigma, mut_prob = 1):
	new_matrix = {}

	#for each value, random sample around starting value
	for aa1 in starting_matrix:
		if aa1 not in new_matrix:
			new_matrix[aa1] = {}
		for aa2 in starting_matrix[aa1]:
			#ensure symmetry, check if this pair exists 
			if aa2 in new_matrix and aa1 in new_matrix[aa2]:
				new_matrix[aa1][aa2] = new_matrix[aa2][aa1]
			#if it doesn't, check what to do next
			else:
				#check if we should mutate
				mut = np.random.uniform(0,1)
				#if we should, select from normal distribution around starting value 
				if mut < mut_prob:
					new_matrix[aa1][aa2] = np.random.normal(starting_matrix[aa1][aa2], sigma)
				#keep value from start matrix
				else:
					new_matrix[aa1][aa2] = starting_matrix[aa1][aa2]
	check_symmetry(new_matrix)
	return new_matrix

def create_matrix_file(matrix, filename):
	aas = list(matrix.keys())
	with open(filename, 'w') as f:
		f.write(" ".join(aas) + '\n')
		for aa1 in aas:
			f.write(" ".join([str(matrix[aa1][aa2]) for aa2 in aas]) + '\n')


def genetic_loop(pos, neg, population, calculated_scores, generation, goal):
	"""
	Selection
    Crossover
    Mutation
    Compute fitness
	"""

	#score them and select the two most "fit"
	scores = {}
	for i, individual in enumerate(population):
		#check if we already have the score for this matrix, no need to recalculate
		#the parents that didn't change
		if calculated_scores[i] != -1:
			scores[i] = calculated_scores[i]
			continue
		
		#calculate scores
		true_scores = [x[1] for x in 
			calc_all_scores(pos, individual, 11, 1)]
		false_scores = [x[1] for x in 
			calc_all_scores(neg, individual, 11, 1)]

		#evaluate current 
		scores[i] = loss_function(true_scores, false_scores)

		#if this matrix scored better, we're done!
		if scores[i] > goal:
			return [individual], (i, scores[i])

	print("Gen " + str(generation), scores.values())
	#keep top 2
	pop_to_keep = sorted(scores, key = lambda x: scores[x], reverse = True)[:2]
	population = [population[i] for i in pop_to_keep] #get matrices
	new_scores = [scores[i] for i in pop_to_keep] #get corresponding scores

	#take the top 2 as the 'parents' for the next generation
	p1, p2 = copy.deepcopy(population[0]), copy.deepcopy(population[1])

	#crossover
	#randomly select n aas
	aas = list(p1.keys())
	switch_aas = np.random.choice(aas, np.random.randint(1, 23), replace = False)
	
	#switch those
	for aa1 in aas:
		for aa2 in p1[aa1]:
			if aa1 in switch_aas:
				#switch them, ensure symmetry
				temp1, temp2 = p1[aa1][aa2], p1[aa2][aa1]
				p1[aa1][aa2], p1[aa2][aa1] = p2[aa1][aa2], p2[aa2][aa1] 
				p2[aa1][aa2], p2[aa2][aa1] = temp1, temp2

	#ensure symmetry
	check_symmetry(p1)
	check_symmetry(p2)

	#mutation
	#mutate 10% of the scores slightly
	o1 = mutate_matrix(p1, 0.2, 0.1)
	o2 = mutate_matrix(p2, 0.2, 0.1)
	
	population += [o1, o2]
	new_scores += [-1, -1]

	return population, new_scores

def optimize_score_matrix(pos, neg, starting_matrix, goal = None):
	"""
	Genetic, starting from https://towardsdatascience.com/introduction-to-genetic-algorithms-
	including-example-code-e396e98d8bf3
	"""
	starting_matrix = algs.get_scoring_matrix(starting_matrix)

	#calculate the starting score 
	true_scores = [x[1] for x in calc_all_scores(pos, starting_matrix, 11, 1)]
	false_scores = [x[1] for x in calc_all_scores(neg, starting_matrix, 11, 1)]
	scores = [loss_function(true_scores, false_scores)]

	#if we dont have a specific goal, beat the starting matrix
	if goal is None:
		goal = scores[0]
	print("Starting score:", scores[0])

	population = [starting_matrix]

	#generate 3 random matrices for a total population of 4
	for i in range(3):
		population.append(mutate_matrix(population[i], 0.5))
		scores.append(-1)
	
	generation = 0
	best_matrix = None
	while best_matrix is None:
		population, scores = genetic_loop(pos, neg, population, scores, generation, goal)
		#check if we found the best
		if len(population) == 1:
			i, score = scores
			print ("Optimized: {0}_{1}, score: {2}".format(generation, i, score))
			best_matrix = population[0]
			create_matrix_file(best_matrix, "optimized_" + starting_matrix)
		generation += 1

	return best_matrix

if __name__ == "__main__":
	pos = read_pairs("Pospairs.txt")
	neg = read_pairs("Negpairs.txt")

	optimize_score_matrix(pos, neg, "BLOSUM50", 2.5)
