def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    # for looooooops
    # rows = 
    # raise NotImplementedError()
    col = input("Seq A: ")
    row = input("Seq B: ")
    n = len(col) + 1
    m = len(row) + 1

    print(40 * '-')
    print("Scoring Scheme: ")
    gap = int(input("Gap penalty: "))
    penalty = int(input("Penalty: "))

    # intial values for scoring matrix
    score_matrix = [[0 for j in range(m)] for i in range(n)]
    for i in range(1, n):
        score_matrix[i][0] = score_matrix[i-1][0] + gap

    for j in range(1, m):
        score_matrix[0][j] = score_matrix[0][j-1] + gap

    for row in score_matrix:
        print(row)
    print(40 * '-')

    # fill in scoring matrix
    for i in range(1, n):
        for j in range(1, m):
            if col[i-1] == row[j-1]:
                diag = score_matrix[i-1][j-1]
            else:
                diag = score_matrix[i-1][j-1] + penalty
            up = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap
            score_matrix[i][j] = min(diag, up, left)
    for row in score_matrix:
        print(row)


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
