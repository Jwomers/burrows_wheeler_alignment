Burrows Wheeler Alignment in Python
=========================

This is a python implementation of the Burrows Wheeler Alignment for DNA sequence matching. I wrote this to better understand the BWA algorithm, and to practice my Python skills.
* It does NOT implement any of the heuristic methods included in the actual BWA algorithm, so implements the basic algorithm which returns 100% accurate results
* It uses the Burrows-Wheeler transform (BWT), Suffix Array (SA), and 2 other auxillary datastructures C and Occ (sometimes called O)
* It uses the D array to prune the tree search and speed up the inexact search algorithm.
* The search is case insensitive.

Differences between this code and the real BWA algorithm
=========================
* This is **not in any way optimised**. This has been coded for ease of understanding, ease of reading and with the goal of better understaning the basic BWA algorithm and its datastructures
* This code parses the suffix tree using a recursive depth-first search, while the real BWA uses breadth-first search (and a heap datastructure), and uses this to prioritise the best partial alignments first
* If BWA finds a result with a difference score of z, it only further considers matches with difference z+1, speeding up the process by ignoring worse results and pruning the search space more aggresively
* BWA sets a maximum allowed difference in the first few tens of bases (called the seed) resulting in 2.5x improvement, and very little loss in accuracy. This is not effective on shorter reads (<50bp)
* BWA reduces the required operating memory by storing small fractions of the Occ and SA arrays, and calculating the rest on the fly. This implementation does not do this.