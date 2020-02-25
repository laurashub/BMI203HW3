# HW3

[![Build
Status](https://travis-ci.org/laurashub/BMI203HW3.svg?branch=master)](https://travis-ci.org/laurashub/BMI203HW3)

## functions
* __main__.py
  * read__pairs
  * read_fasta
  * read_alignment_scores
  * read_existing_aligns
  * calc_all_scores
  * calc_all_aligns
  * calc_fp_rate
  * best_fp_rate
  * compare_matrices
  * compare_normalized
  * loss_function
  * check_symmetry
  * mutate_matrix
  * create_matrix_file
  * genetic_loop
  * optimize_score_matrix
  * full_optimization_run
  * get_score_diff
  * get_len_diff
* algs.py
  * sw
  * sw_traceback
  * sanitize_seq
  * get_scoring_matrix
  * align
  * score
  * apply
  * roc
  * single_roc

## main

Currently does nothing, was used to generate alignments and graphs

```
python -m smith_waterman
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.
