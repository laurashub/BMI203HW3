from smith_waterman import algs
import numpy as np
import seaborn as sns
import pandas as pd
import copy
import os

def read_pairs(filename):
	"""
	Read paires from pospairs or negpairs
	"""
	pairs = []
	with open(filename, 'r') as pairs_file:
		for line in pairs_file.read().splitlines():
			pairs.append(tuple(map(lambda x: read_fasta(x), line.split())))
	return pairs

def read_fasta(filename):
	"""
	Get sequence from FASTA file
	"""
	seq = ""
	for line in open(filename, 'r').read().splitlines():
		if line[0] is not '>':
			seq += line
	return seq

def read_alignment_scores(filename):
	"""
	get all the scores in filename, previously calculated scores
	"""
	scores = []
	with open(filename, 'r') as f:
		for line in f.read().splitlines():
			try:
				scores.append(float(line))
			except:
				pass
	return scores

def read_existing_aligns(filename):
	"""
	Get previously calculated alignments
	"""
	pairs = []
	with open(filename, 'r') as f:
		lines = f.read().splitlines()
	for i in range(0, len(lines), 4):
		pairs.append((lines[i], lines[i+1]))
	return pairs

def calc_all_scores(pairs, matrix, gap_start, gap_extend):
	"""
	Calculate alignment score for all pairs in pairs
	"""
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
	"""
	Calculate and save alignments and score for all pairs in pairs
	"""
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
	"""
	For a given tp cutoff, calculate the fp rate
	"""
	true_scores = sorted(true_scores)

	#get cutoff value
	cutoff = true_scores[int(len(true_scores)*(1-cutoff))]

	fps = 0
	for score in false_scores:
		if score > cutoff:
			fps += 1

	return fps/len(false_scores)


def best_fp_rate(pos, neg, matrix = "BLOSUM50"):
	"""
	For a range of gap starts and gap extends, plot fp rate as a heatmap
	"""
	fp_rates = pd.DataFrame(index=range(1,21), columns=range(1,6), dtype=np.double)

	#loop through gap extend and gap start values
	for gap_start in range(1,21):
		for gap_extend in range(1,6):
			true_scores = []

			#calculate true and false scores with specified open/extend
			true_scores = calc_all_scores(pos, matrix, 
					gap_start = gap_start, gap_extend = gap_extend)
			false_scores =  calc_all_scores(neg, matrix, 
					gap_start = gap_start, gap_extend = gap_extend)
			
			#calculate fp rate from those scores
			fp_rates[gap_extend][gap_start] = calc_fp_rate(
				[x[1] for x in true_scores], 
				[y[1] for y in false_scores], 0.7)

	#plot fp rates
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

	"""
	Takes in pos pairs, neg pairs, matrices, plots ROC curve with all provided matrices
	"""

	TSs = []
	FSs = []

	if prescored is None:
		prescored = []
		for m in matix:
			prescored.append(None)


	for i, matrix in enumerate(matrices):
		if prescored[i] is not None and prescored[0] == i:
			TSs.append(prescored[1][0])
			FSs.append(prescored[1][1])
			continue
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
	"""
	Calculate ROC curve for raw and normalized scores for a given matrix and pos/neg pairs
	"""
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
		#Get cutoff value
		TPs = 0.0
		cutoff = false_scores[int(len(false_scores)*(1-fp_rate))-1]
		#add up true positives
		for score in true_scores:
			if score > cutoff:
				TPs += 1
		fp_score += (TPs/len(true_scores))
	return fp_score

def check_symmetry(matrix):
	#make sure matrix is symmetrical
	for aa1 in matrix:
		for aa2 in matrix[aa1]:
			assert(matrix[aa1][aa2] == matrix[aa2][aa1])

def mutate_matrix(starting_matrix, sigma, mut_prob = 1):
	"""
	"Mutate" a few values in the matrix to introduce diversity into population
	"""
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
	#make sure its still symmetric
	check_symmetry(new_matrix)
	return new_matrix

def create_matrix_file(matrix, filename):
	"""
	Write matrix to a file that can be interpreted by algs.get_scoring_matrix
	"""
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
		
		#calculate scores on each of the alignments
		true_scores = [algs.score_existing_align(*(p), individual) for p in pos]
		false_scores = [algs.score_existing_align(*(n), individual) for n in neg]

		#evaluate current 
		scores[i] = loss_function(true_scores, false_scores)

		#if this matrix scored better, we're done!
		if scores[i] > goal:
			return [individual], (i, scores[i])

	#print("Gen " + str(generation), list(scores.values()))
	#keep top 2 as parents
	pop_to_keep = sorted(scores, key = lambda x: scores[x], reverse = True)[:2]
	population = [population[i] for i in pop_to_keep] #get matrices
	new_scores = [scores[i] for i in pop_to_keep] #get corresponding scores

	#take the top 2 as the 'parents' for the next generation
	p1, p2 = copy.deepcopy(population[0]), copy.deepcopy(population[1])

	#crossover
	#randomly select n aas
	aas = list(p1.keys())
	switch_aas = np.random.choice(aas, np.random.randint(1, 23), replace = False)
	
	#switch the values between the parents for those AAs
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
	o1 = mutate_matrix(p1, 0.5, 0.1)
	o2 = mutate_matrix(p2, 0.5, 0.1)
	
	#append the offspring to the new population, give placeholder scores
	population += [o1, o2]
	new_scores += [-1, -1]

	return population, new_scores

def optimize_score_matrix(pos, neg, starting_matrix, goal = None, max_gen = 100):
	"""
	Genetic, starting from https://towardsdatascience.com/introduction-to-genetic-algorithms-
	including-example-code-e396e98d8bf3
	"""
	starting_mat = algs.get_scoring_matrix(starting_matrix)

	#calculate the starting score 
	true_scores = [algs.score_existing_align(*(p), starting_mat) for p in pos]
	false_scores = [algs.score_existing_align(*(n), starting_mat) for n in neg]
	
	scores = [loss_function(true_scores, false_scores)]

	#if we dont have a specific goal, beat the starting matrix
	if goal is None:
		goal = scores[0]
	print("Starting score:", scores[0])

	population = [starting_mat]

	#generate 3 random matrices for a total population of 4
	for i in range(3):
		population.append(mutate_matrix(population[i], 0.5))
		scores.append(-1)
	
	generation = 0
	best_matrix = None
	matrix_filename = None
	while best_matrix is None:
		#get current generation and scores
		population, scores = genetic_loop(pos, neg, population, scores, generation, goal)
		#check if we found the best
		if len(population) == 1:
			#save matrix and return
			i, score = scores
			print ("Optimized: {0}_{1}, score: {2}".format(generation, i, score))
			best_matrix = population[0]
			create_matrix_file(best_matrix, "optimized_" + starting_matrix)
			matrix_filename = "optimized_" + starting_matrix
		#hit max without hitting goal
		if generation == max_gen:
			#get the best matrix in the current population
			best_score = -1
			for individual, score in zip(population, scores):
				if score > best_score:
					best_score = score
					best_matrix = individual
			create_matrix_file(best_matrix, "optimized_" + starting_matrix)
			matrix_filename = "optimized_" + starting_matrix
		generation += 1

	return best_matrix, matrix_filename

def full_optimization_run(starting_pos, starting_neg, starting_matrix, goal = 4, max_gen = 5000):
	
	#optimize starting
	pos_aligns = read_existing_aligns(starting_pos)
	neg_aligns = read_existing_aligns(starting_neg)


	#get best matrix and filename
	new_matrix, matrix_filename = optimize_score_matrix(pos_aligns, neg_aligns, starting_matrix, goal = goal, max_gen = max_gen)
	
	#generate scores in initial alignment
	starting_mat = algs.get_scoring_matrix(starting_matrix)

	#calculate the starting score 
	true_scores = [algs.score_existing_align(*(p), starting_mat) for p in pos_aligns]
	false_scores = [algs.score_existing_align(*(n), starting_mat) for n in neg_aligns]

	#calculate score on previous alignment
	precal_pos = [algs.score_existing_align(*(p), new_matrix) for p in pos_aligns]
	precal_neg = [algs.score_existing_align(*(n), new_matrix) for n in neg_aligns]


	#read in initial pairs for realign
	pos = read_pairs("Pospairs.txt")
	neg = read_pairs("Negpairs.txt")

	#save new aligns
	calc_all_aligns(pos, new_matrix, 11, 1, starting_matrix + "_optimized_pos_aligns.txt")
	calc_all_aligns(neg, new_matrix, 11, 1, starting_matrix + "_optimized_neg_aligns.txt")

	#read in and score new aligns for consistency in scoring
	pos_aligns_2 = read_existing_aligns(starting_matrix + "_optimized_pos_aligns.txt")
	neg_aligns_2 = read_existing_aligns(starting_matrix + "_optimized_neg_aligns.txt")

	precal_pos_2 = [algs.score_existing_align(*(p), new_matrix) for p in pos_aligns_2]
	precal_neg_2 = [algs.score_existing_align(*(n), new_matrix) for n in neg_aligns_2]

	prescored = [(0,(true_scores, false_scores)), (1,(precal_pos, precal_neg)), (2,())]

	TSs = [true_scores, precal_pos, precal_pos_2]
	FSs = [false_scores, precal_neg, precal_neg_2]

	#plot all new scores
	algs.roc(TSs, FSs, 
		[starting_matrix, starting_matrix + "_optimized", starting_matrix + "_optimized_realign"], 
		save = True, filename = starting_matrix + "_optimization.png")

def get_score_diff(matrix):
	"""
	get score on starting alignments using initial matrix and optimized matrix
	"""
	pos_start = read_existing_aligns("BLOSUM50_pos_aligns.txt")
	neg_start = read_existing_aligns('BLOSUM50_neg_aligns.txt')

	#matrices
	start_mat = algs.get_scoring_matrix(matrix)
	op_mat = algs.get_scoring_matrix("optimized_" + matrix)

	#initial scores
	tp_b50 = [algs.score_existing_align(*(p), start_mat) for p in pos_start]
	fp_b50 = [algs.score_existing_align(*(n), start_mat) for n in neg_start]

	#optimized scores
	tp_b50o = [algs.score_existing_align(*(p), op_mat) for p in pos_start]
	fp_b50o = [algs.score_existing_align(*(n), op_mat) for n in neg_start]

	all_b50 = tp_b50 + fp_b50
	all_b50o = tp_b50o + fp_b50o

	#calculate and print differences
	total_diff = sum([x[0] - x[1] for x in zip(all_b50o, all_b50)]) / len(all_b50)
	pos_diff = sum([x[0] - x[1] for x in zip(tp_b50o, tp_b50)]) / len(tp_b50)
	neg_diff = sum([x[0] - x[1] for x in zip(fp_b50o, fp_b50)]) / len(fp_b50)

	print(total_diff, pos_diff, neg_diff)

def get_len_diff(initial_pairs, realigned_pairs):
	"""
	Get the differences in length between initial alignment and realignment
	"""
	#initial aligns
	pos_start = read_existing_aligns(initial_pairs[0])
	neg_start = read_existing_aligns(initial_pairs[1])

	#realigned
	pos_realign = read_existing_aligns(realigned_pairs[0])
	neg_realign = read_existing_aligns(realigned_pairs[1])

	#get length of everything
	splen = [len(x[0]) for x in pos_start]
	sflen = [len(x[0]) for x in neg_start]
	salen = splen + sflen

	rplen = [len(x[0]) for x in pos_realign]
	rflen = [len(x[0]) for x in neg_realign]
	ralen = rplen + rflen

	#calculate differences
	total_diff = sum([x[0] - x[1] for x in zip(salen, ralen)]) / len(salen)
	pos_diff = sum([x[0] - x[1] for x in zip(splen, rplen)]) / len(rplen)
	neg_diff = sum([x[0] - x[1] for x in zip(sflen, rflen)]) / len(rflen)

	print(total_diff, pos_diff, neg_diff)


if __name__ == "__main__":
	pass
	#full_optimization_run('BLOSUM50_pos_aligns.txt', 'BLOSUM50_neg_aligns.txt', "MATIO", max_gen = 5000)

	#get_len_diff(('BLOSUM50_pos_aligns.txt', 'BLOSUM50_neg_aligns.txt'),
	#("MATIO_optimized_pos_aligns.txt", "MATIO_optimized_neg_aligns.txt"))
	#mat = algs.get_scoring_matrix("BLOSUM50")
	#algs.score_existing_align("AAAAAA", "AAAAA", mat)
	
	"""
	#optimize starting
	#read in and score new aligns for consistency in scoring
	pos_start = read_existing_aligns("BLOSUM50_pos_aligns.txt")
	neg_start = read_existing_aligns('BLOSUM50_neg_aligns.txt')

	cur_matrix = "MATIO"
	tp_b50 = [algs.score_existing_align(*(p), "MATIO") for p in pos_start]
	fp_b50 = [algs.score_existing_align(*(n), "MATIO") for n in neg_start]
	tp_b50o = read_alignment_scores("optimized_MATIO")
	fp_b50o = read_alignment_scores("optimized_MATO")

	total_diff = sum([x[0] - x[1] for x in zip(all_b50o, all_b50)]) / len(all_b50)
	pos_diff = sum([x[0] - x[1] for x in zip(tp_b50o, tp_b50)]) / len(tp_b50)
	neg_diff = sum([x[0] - x[1] for x in zip(fp_b50o, fp_b50)]) / len(fp_b50)

	print(total_diff, pos_diff, neg_diff)
	
	pos_aligns = read_existing_aligns(starting_pos)
	neg_aligns = read_existing_aligns(starting_neg)

	tp_b50 = read_alignment_scores("BLOSUM50_pos_aligns.txt")
	fp_b50 = read_alignment_scores("BLOSUM50_neg_aligns.txt")
	all_b50 = tp_b50 + fp_b50

	tp_b50o = read_alignment_scores("BLOSUM50_optimized_pos_aligns.txt")
	fp_b50o = read_alignment_scores("BLOSUM50_optimized_neg_aligns.txt")
	all_b50o = tp_b50o + fp_b50o

	total_diff = sum([x[0] - x[1] for x in zip(all_b50o, all_b50)]) / len(all_b50)
	pos_diff = sum([x[0] - x[1] for x in zip(tp_b50o, tp_b50)]) / len(tp_b50)
	neg_diff = sum([x[0] - x[1] for x in zip(fp_b50o, fp_b50)]) / len(fp_b50)

	print(total_diff, pos_diff, neg_diff)
	"""
	
