import numpy as np
from smith_waterman import algs

def test_roc():
    return None

def test_read_matrix():
	blosum50 = algs.get_scoring_matrix("BLOSUM50")
	assert(blosum50['A']['A'] == 5)
	assert(blosum50['F']['A'] == -3)
	assert(blosum50['A']['F'] == -3)

	blosum62 = algs.get_scoring_matrix("BLOSUM62")
	assert(blosum62['A']['A'] == 4)
	assert(blosum62['F']['A'] == -2)
	assert(blosum62['A']['F'] == -2)

	matio = algs.get_scoring_matrix("MATIO")
	assert(matio['A']['A'] == 0)
	assert(matio['F']['A'] == 2)
	assert(matio['A']['F'] == 2)

def test_smithwaterman():
	score, seq1, seq2, align = algs.align("AAAAAAA", "AAAAAAA", "BLOSUM50")
	assert(seq1 == seq2 == "AAAAAAA")

	score, seq1, seq2, align = algs.align("CAGT", "CAAAGT", "BLOSUM50")
	assert(seq1 == "CA--GT" or seq1 == "C--AGT") #affine, should have double gap
	assert(seq2 == "CAAAGT")

def test_scoring():
	#7 * A:A (5)
    assert (algs.score("AAAAAAA", "AAAAAAA", "BLOSUM50") == 35)
    assert (algs.score("AAAAAAA", "GGGGGGG", "BLOSUM50") == 0)
