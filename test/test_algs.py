import numpy as np

from smith_waterman import algs
import smith_waterman.__main__ as sw


def test_roc():
	pos = sw.read_pairs("Pospairs.txt")
	neg = sw.read_pairs("Negpairs.txt")

	true_scores = sw.calc_all_scores(pos, "BLOSUM50", 11, 1)
	false_scores = sw.calc_all_scores(neg, "BLOSUM50", 11, 1)

	xs, ys = algs.roc([[x[1] for x in true_scores]],
		[[x[1] for x in false_scores]], ["BLOSUM50"])

	assert all([all([x <= 1 and x >= 0 for x in xvals]) for xvals in xs])
	assert all([all([y <= 1 and y >= 0 for y in yvals]) for yvals in ys])

def test_read_matrix():
	blosum50 = algs.get_scoring_matrix("BLOSUM50")
	assert(blosum50['A']['A'] == 5)
	assert(blosum50['F']['A'] == -3)
	assert(blosum50['A']['F'] == -3)

	blosum62 = sw.algs.get_scoring_matrix("BLOSUM62")
	assert(blosum62['A']['A'] == 4)
	assert(blosum62['F']['A'] == -2)
	assert(blosum62['A']['F'] == -2)

	matio = sw.algs.get_scoring_matrix("MATIO")
	assert(matio['A']['A'] == 0)
	assert(matio['F']['A'] == 2)
	assert(matio['A']['F'] == 2)

def test_smithwaterman():
	score, seq1, seq2, align = sw.algs.align("AAAAAAA", "AAAAAAA", "BLOSUM50")
	assert(seq1 == seq2 == "AAAAAAA")

	score, seq1, seq2, align = sw.algs.align("CAGT", "CAAAGT", "BLOSUM50")
	assert(seq1 == "CA--GT" or seq1 == "C--AGT") #affine, should have double gap
	assert(seq2 == "CAAAGT")

def test_scoring():
	#7 * A:A (5)
    assert (algs.score("AAAAAAA", "AAAAAAA", "BLOSUM50") == 35)
    assert (algs.score("AAAAAAA", "GGGGGGG", "BLOSUM50") == 0)
