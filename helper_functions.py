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
    gap_score = -1
    rows = len(seq1) + 1
    columns = len(seq2) + 1
    score_matrix = [[0] * columns for _ in range(rows)]

    for i in range(1, rows):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap_score
    for j in range(1, columns):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap_score

    for i in range(1, rows):
        for j in range(1, columns):
            diagonal = score_matrix[i - 1][j - 1] + scoring_function(
                seq1[i - 1], seq2[j - 1]
            )
            up = score_matrix[i - 1][j] + gap_score
            left = score_matrix[i][j - 1] + gap_score
            score_matrix[i][j] = max(diagonal, up, left)

    aligned_seq1 = []
    aligned_seq2 = []
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and score_matrix[i][j]
            == score_matrix[i - 1][j - 1]
            + scoring_function(seq1[i - 1], seq2[j - 1])
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_score:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    return (
        "".join(reversed(aligned_seq1)),
        "".join(reversed(aligned_seq2)),
        float(score_matrix[-1][-1]),
    )


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
