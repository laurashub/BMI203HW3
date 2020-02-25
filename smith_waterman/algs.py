import numpy as np
import sys
import matplotlib.pyplot as plt

def sw(seq1, seq2, score_matrix, gap_start, gap_extend, align = True):

	#pass in either a dict or a str
	if isinstance(score_matrix, str):
		#import scoring function
		scoring_matrix = get_scoring_matrix(score_matrix)
	else:
		scoring_matrix = score_matrix

	#append a space to the front, makes the lengths corrent for 
	#the matrix dimensions and avoids off by 1 errors later
	#Also fix certain amino acids (x) that were giving key errors
	seq1 = " " + sanitize_seq(seq1, scoring_matrix)
	seq2 = " " + sanitize_seq(seq2, scoring_matrix)

	#initialize 3 matrices
	#First matrix holds the scores for a local alignment ending 
	#in a match (or mismatch) as opposed to a gap
	M = np.zeros((len(seq1), len(seq2))) 
	Mp = np.zeros((len(seq1), len(seq2))) 

	#Second matrix holds scores for local alignment ending in a gap
	#in seq1
	X = np.zeros((len(seq1), len(seq2)))
	Xp = np.zeros((len(seq1), len(seq2))) 

	#Second matrix holds scores for local alignment ending in a gap
	#in seq2
	Y = np.zeros((len(seq1), len(seq2))) #ending w space in Y
	Yp = np.zeros((len(seq1), len(seq2))) #ending w space in Y

	matrices = [(M, Mp), (X, Xp), (Y, Yp)]

	for i, aa1 in enumerate(seq1):
		for j, aa2 in enumerate(seq2):

			#keep first row and cols to 0
			if i == 0 or j == 0:
				continue
			#update scoring matrices
			#Add score of current match to the best scoring local alignment
			#up to this point. No gaps, so no penalties applied
			M[i,j] = scoring_matrix[aa1][aa2] + max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1])

			#Pointer matrix - get the matrix that we came from by getting argmax
			Mp[i,j] = np.argmax([M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1]])

			#Gap in seq1
			#If coming from the other two matrices (M or Y), this is
			#the start of a gap, so apply gap_start penalty
			#If coming from the same matrix (X), this is an extension of
			#an already existing gap, so apply gap extend penalty
			X[i,j] = max(M[i, j-1] - gap_start,
						X[i, j-1] - gap_extend, 
						Y[i, j-1] - gap_start)
			#Pointer matrix - get the matrix that we came from by getting argmax
			Xp[i,j] = np.argmax([M[i, j-1] - gap_start,
						X[i, j-1] - gap_extend, 
						Y[i, j-1] - gap_start])

			#Gap in seq2
			#If coming from the other two matrices (M or X), this is
			#the start of a gap, so apply gap_start penalty
			#If coming from the same matrix (X), this is an extension of
			#an already existing gap, so apply gap extend penalty
			Y[i,j] = max(M[i - 1, j] - gap_start,
						X[i - 1, j] - gap_start,
						Y[i - 1, j] - gap_extend)
			#Pointer matrix - get the matrix that we came from by getting argmax
			Yp[i,j] = np.argmax([M[i - 1, j] - gap_start,
						X[i - 1, j] - gap_start,
						Y[i - 1, j] - gap_extend])
	
	#do full traceback
	if align:
		results = sw_traceback(seq1, seq2, matrices)

	#just want the score, max val across all 3 matrices
	else:
		results = max(np.amax(M), np.amax(X), np.amax(Y))
	#score for the alignment is the best score across all 3 matrices
	return results

def sw_traceback(seq1, seq2, matrices):
	#traceback to reconstruct the alignments

	#unpack matrices
	(M, Mp), (X, Xp), (Y, Yp) = matrices

	#get largest value across all matrices and indeces
	#this is the higherest scoring local alignment, and thus
	#the starting point for tracing back the local alignment
	max_mat = np.argmax([np.amax(M), np.amax(X), np.amax(Y)])
	cur_matrix, cur_pointer = matrices[max_mat]
	i, j = np.unravel_index(cur_matrix.argmax(), cur_matrix.shape)

	seq1_aligned = ""
	seq2_aligned = ""
	align = ""

	cur_score = cur_matrix[i, j]
	best_score = cur_score

	#get the next matrix to go to
	prev_matrix = matrices[int(cur_pointer[i,j])][0]

	#Go until 0 or until you reach the beginning of one of the
	#sequences
	while cur_matrix[i, j] > 0 and (i > 0 and j > 0):
		#adapted from explanation at
		# https://stackoverflow.com/questions/18101078/
		# traceback-in-smith-wateman-algorithm-with-affine-gap-penalty

		#save previous matrix for current position + next pointer
		#need to do this here because i and j are about to change
		prev_matrix, cur_pointer = matrices[int(cur_pointer[i,j])]

		#Update sequence based on current matrix
		#If we're currently in M, we matched the values,
		#add current ones to alignment
		#Also this means we were moving diagonally, so subtract
		#both i and j before switching matrices
		if cur_matrix is M:
			seq1_aligned = seq1[i] + seq1_aligned
			seq2_aligned = seq2[j] + seq2_aligned
			if seq1[i] == seq2[j]:
				align = "|" + align
			else:
				align = "*" + align
			i -= 1
			j -= 1
		#if we're in X, add a gap in seq1
		#This also means that we are moving horizontally, so 
		#subtract j before switching matrices
		elif cur_matrix is X:
			seq1_aligned = '-' + seq1_aligned
			seq2_aligned = seq2[j] + seq2_aligned
			align = " " + align
			j -= 1
		#Similarly if we're in Y, add a gap in seq2
		#Were moving vertically, subtract from i
		elif cur_matrix is Y:
			seq1_aligned = seq1[i] + seq1_aligned
			seq2_aligned = '-' + seq2_aligned
			align = " " + align
			i -= 1

		#Update our current matrix based on what our pointer matrix said
		#ie trace back to the previous position
		cur_matrix = prev_matrix

	return best_score, seq1_aligned, seq2_aligned, align

def sanitize_seq(seq, allowed):
	#Get rid of disallowed characters to avoid key errors, also set to uppercase
	new_seq = ""
	for char in seq.upper():
		if char not in allowed:
			new_seq += "*"
		else:
			new_seq += char
	return new_seq

def get_scoring_matrix(filename):
	#open file, skip comment lines, build dictionary of scores
	scores = {}
	with open(filename, 'r') as f:
		lines = [line for line in f.read().splitlines() if line[0] != "#"]
	aas = lines[0].split()
	lines = lines[1:]
	for i, a1 in enumerate(aas):
		for j, a2 in enumerate(aas):
			if a1 not in scores:
				scores[a1] = {}
			scores[a1][a2] = float(lines[i].split()[j])
	return scores

def align(seq1, seq2, score_matrix = "BLOSUM50", gap_start = 11, gap_extend = 1):
	"""
	Do full align, return aligned sequences
	"""
	return sw(seq1, seq2, score_matrix, gap_start, gap_extend, align = True)

def score(seq1, seq2, score_matrix = "BLOSUM50", gap_start = 11, gap_extend = 1):
	"""
	Align but no traceback, just get the alignment score
	"""
	return sw(seq1, seq2, score_matrix, gap_start, gap_extend, align = False)

def score_existing_align(seq1, seq2, matrix, gap_start = 11, gap_extend = 1):
	"""
	Scores an existing alignment between seq1 and seq2 by pairwise comparison
	"""

	#start at 0, no gap
	score = 0
	gap = False
	
	#go through pairs
	for char1, char2 in zip(seq1, seq2):
		#if this is a gp
		if char1 == "-"or char2 == "-":
			#is were already in a gap, extend
			if gap:
				score -= gap_extend
			#otherwise, this is a new gap, gap start
			else:
				score -= gap_start
				gap = True
		#otherside, add score, turn off gap
		else:
			score += matrix[char1][char2]
			if gap:
				gap = False
	return score


def roc(tss, fss, titles, save = False, show = False, filename = None):
	"""
	Calculates ROC curve for multiple matrices
	"""
	xs, ys = [], []

	#get x and y values for each true and false scores
	for ts, fs in zip(tss, fss):
		x, y = single_roc(ts, fs)
		xs.append(x)
		ys.append(y)

	#plot each line
	for x, y, title in zip(xs, ys, titles):
		label = "{0}: {1:.2f}".format(title, np.trapz(y[::-1],x[::-1]))
		plt.plot(x, y, label = label)

	#middle line for comparison
	plt.plot(np.linspace(0, 1, 100),np.linspace(0, 1, 100), linestyle='--')
	plt.legend()

	#determine save or show
	if save and filename:
		plt.savefig(filename, dpi=200)
	elif show:
		plt.show()
	else:
		return xs, ys


def single_roc(true_scores, false_scores):
	x = []
	y = []
	total_true = len(true_scores)
	total_false = len(false_scores)

	#calculate values for threshold, 100 total steps
	min_score = min(true_scores + false_scores)
	max_score = max(true_scores + false_scores)
	step = (max_score - min_score) / 100

	#loop through thresholds
	for threshold in np.arange(min_score-step, max_score+step, step):
		#start with true positive and false positive rates of 0
		tp = 0.0
		fp = 0.0

		#calculate TP rate for current threshold
		for ts in true_scores:
			if ts > threshold:
				tp += 1
		y.append(tp/total_true)

		#Calculate FP rate for current threshold
		for fs in false_scores:
			if fs > threshold:
				fp += 1
		x.append(fp/total_false)
	#return all calculated true positive and false positive rates for current scores
	return x,y
