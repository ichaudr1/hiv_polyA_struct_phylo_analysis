'''
@author ichaudr
Takes a list of base pairs in two RNA hairpin structres, and returns an alignment of the helical base pairs. 
'''

from itertools import product
from collections import deque
import csv


def rhelix_align(hp_one, hp_two):
    '''
    Returns an alignment between two RNA hairpin structures. The algorithm mimics the Needleman-Wunsch algorithm for sequence alignment; however, it has been adapted for aligning base pairs in a RNA helix. 
    Parameters
    ----------
    hp_one, hp_two: The two hairpins of interest. These must be a list of tuples representing the base pairs in the helix from "bottom to top".
    For example, the hairpin below would be passed in as: [ (A, U) , (G, C), (C, G), (U, A), (U, -), (U - A)]
    U-A
    U
    U-A
    C-G
    G-C
    A-U
    '''

    hp_one_len, hp_two_len = len(hp_one), len(hp_two)

    #Constants that hold values for the dx and dy during the traceback. 
    DIAG = (-1,-1)
    UP = (0, -1)
    LEFT = (-1, 0)

    #Set up two tables: one to hold the values of each posistion and another to hold the instructions for how to backtrace that movement (i.e. Diagonal, DOWN->UP, and RIGHT->LEFT)
    scores = {}
    tracing_pointers = {}

    scores[-1, -1] = 0

    for i in range(hp_one_len):
        scores[i, - 1] = -(i+1)
    for j in range(hp_two_len):
        scores[-1, j] = -(j+1)

    #The current table of scores looks like this now: 
    #  0  -1  -2  -3  -4 .... -hp_one_len
    # -1
    # -2
    # -3
    # .
    # .
    # .
    # -len_y

    #Calculate the scores and populate the scores and tracing_pointers list
    pointer_options = DIAG, LEFT, UP
    for i, j in product(range(hp_one_len), range(hp_two_len)):

        possible_scores = (
            scores[i - 1, j - 1] + get_base_pair_score(hp_one[i], hp_two[j]),
            scores[i - 1, j] - 1,
            scores[i, j - 1] - 1
        )

        scores[i, j] , tracing_pointers[i, j] = max(zip(possible_scores, pointer_options))

    #Find the optimal path to determine the alignment
    alignment = deque()

    i, j = hp_one_len - 1, hp_two_len - 1

    while i >= 0 and j >= 0:
        trace_dir = tracing_pointers[i, j]

        if trace_dir == DIAG:
            alignment_element = (i, j)
        elif trace_dir == LEFT:
            alignment_element = (i, None)
        elif trace_dir == UP:
            alignment_element = (None, j)
        
        alignment.appendleft(alignment_element)

        i = i + trace_dir[0]
        j = j + trace_dir[1]

    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1
    
    return (list(alignment), scores[hp_one_len - 1, hp_two_len - 1])


def rhelix_align_verbose(hp_one, hp_two):
    '''
    Returns an alignment between two RNA hairpin structures. The algorithm mimics the Needleman-Wunsch algorithm for sequence alignment; however, it has been adapted for aligning base pairs in a RNA helix.

    **Prints out the dynamic programming grid to follow intermediate scoring**

    Parameters
    ----------
    hp_one, hp_two: The two hairpins of interest. These must be a list of tuples representing the base pairs in the helix from "bottom to top".
    For example, the hairpin below would be passed in as: [ (A, U) , (G, C), (C, G), (U, A), (U, -), (U - A)]
    U-A
    U
    U-A
    C-G
    G-C
    A-U
    '''

    hp_one_len, hp_two_len = len(hp_one), len(hp_two)

    #Constants that hold values for the dx and dy during the traceback. 
    DIAG = (-1,-1)
    UP = (0, -1)
    LEFT = (-1, 0)

    #Set up two tables: one to hold the values of each posistion and another to hold the instructions for how to backtrace that movement (i.e. Diagonal, DOWN->UP, and RIGHT->LEFT)
    scores = {}
    tracing_pointers = {}

    scores[-1, -1] = 0

    for i in range(hp_one_len):
        scores[i, - 1] = -(i+1)
    for j in range(hp_two_len):
        scores[-1, j] = -(j+1)

    #The current table of scores looks like this now: 
    #  0  -1  -2  -3  -4 .... -hp_one_len
    # -1
    # -2
    # -3
    # .
    # .
    # .
    # -len_y

    #Calculate the scores and populate the scores and tracing_pointers list
    pointer_options = DIAG, LEFT, UP
    for i, j in product(range(hp_one_len), range(hp_two_len)):

        possible_scores = (
            scores[i - 1, j - 1] + get_base_pair_score(hp_one[i], hp_two[j]),
            scores[i - 1, j] - 1,
            scores[i, j - 1] - 1
        )

        scores[i, j] , tracing_pointers[i, j] = max(zip(possible_scores, pointer_options))

    #Find the optimal path to determine the alignment
    alignment = deque()

    i, j = hp_one_len - 1, hp_two_len - 1

    while i >= 0 and j >= 0:
        trace_dir = tracing_pointers[i, j]

        if trace_dir == DIAG:
            alignment_element = (i, j)
        elif trace_dir == LEFT:
            alignment_element = (i, None)
        elif trace_dir == UP:
            alignment_element = (None, j)
        
        alignment.appendleft(alignment_element)

        i = i + trace_dir[0]
        j = j + trace_dir[1]

    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1
    

    print(hp_one_len)
    print(hp_two_len)
    print(scores)

    rows = []
    for i in range(-1, hp_one_len):
        row = []
        for j in range(-1, hp_two_len):
            row.append(str(scores[i, j]))
        rows.append(row)
    with open('test.csv', 'w+') as f:
        fwriter = csv.writer(f)
        fwriter.writerows(rows)

    return (list(alignment), scores[hp_one_len - 1, hp_two_len - 1])


def get_base_pair_score(bp_one, bp_two):
    '''
    Calculates the similarity score between two base pairs. If they are identical, the score is 1. If they have the same nucleotides but are on opposite sides of the helix, the score is .75. If they share one nucleotide on the same side of the helix, the score is .5. If they share one nucelotide on opposite sides of the helix, the score is .25. Otherwise, the score is 0. 
    Parameters
    ---------
    bp_one, bp_two: The two base pairs to be compared. These are tuples representing the base pairs: i.e: (A, U)
    '''

    if bp_one[0] == bp_two[0] and bp_one[1] == bp_two[1]:
        return 1
    
    if bp_one[0] in bp_two and bp_one[1] in bp_two:
        return .75
    
    if bp_one[0] == bp_two[0] or bp_one[1] == bp_two[1]:
        return .5
    
    if bp_one[0] in bp_two or bp_one[1] in bp_two:
        return .25
    
    return 0




def get_pairs(sequence, vstring):
    '''
    Returns a list containing the pairs for this sequence and vienna string. 
    Parameters
    ----------
    sequence: string
        The nucleotide sequence (AUCG)
    vstring: string
        The vienna string describing the secondary structure for this sequence. 
    
    Returns
    -------
    pairs: list[(str, str, str, str)]
        A list of tuples for each of the pairs in the structure. The list is in sequential order from "botttom to top" of the hairpin. i.e. ('(', 'C', 'G', ')')
    '''

    #Checks if the entered data is valid
    if not len(sequence) == len(vstring):
        return []

    end_helix_pos = vstring.find(')')
    second_half_struct = vstring[end_helix_pos:]

    for c in second_half_struct:
        if c == '(':
            return []

    #List that holds the pairs
    pairs = []

    #Convert the sequence and vienna structure into a list
    temp_vstring = list(vstring)
    temp_seq = list(sequence)

    #Loop through the sequence and parse outt he pairs
    while len(temp_vstring) > 1:

        #Parse the 5' and 3' symbol and sequence
        five_end_char = temp_vstring[0]
        three_end_char = temp_vstring[-1]

        #If the current 5' and 3' positions belong in a pair, they are added
        if (five_end_char == '(' and three_end_char == ')') or (five_end_char == '.' and three_end_char == '.'):

            pairs.append((temp_seq.pop(0), temp_seq.pop(-1)))
            temp_vstring.pop(0)
            temp_vstring.pop(-1)

        #If either of the positions are unpaired, the pair that is added is adjusted.
        elif five_end_char == '.':

            pairs.append((temp_seq.pop(0), '-'))
            temp_vstring.pop(0)

        elif three_end_char == '.':
            
            pairs.append(('-', temp_seq.pop(-1)))
            temp_vstring.pop(-1)

    while len(temp_vstring) > 0:
        pairs.append(('-', temp_seq.pop(0)))
        temp_vstring.pop(0)
    
    return pairs





##Testing with normal RNA sequences
def test_seq_align(x, y):
    
    len_x = len(x)
    len_y = len(y)

    #Constants that hold values for the dx and dy during the traceback. 
    DIAG = (-1,-1)
    UP = (0, -1)
    LEFT = (-1, 0)

    #Set up two tables: one to hold the values of each posistion and another to hold the instructions for how to backtrace that movement (i.e. Diagonal, DOWN->UP, and RIGHT->LEFT)
    scores = {}
    tracing_pointers = {}

    scores[-1, -1] = 0

    for i in range(len_x):
        scores[i, - 1] = -i
    for j in range(len_y):
        scores[-1, j] = -j

    #The current table of scores looks like this now: 
    #  0  -1  -2  -3  -4 .... -len_x
    # -1
    # -2
    # -3
    # .
    # .
    # .
    # -len_y

    #Calculate the scores and populate the scores and tracing_pointers list
    pointer_options = DIAG, LEFT, UP
    for i, j in product(range(len_x), range(len_y)):

        possible_scores = (
            scores[i - 1, j - 1] + int(x[i] == y[j]),
            scores[i - 1, j] - 1,
            scores[i, j - 1] - 1
        )

        scores[i, j] , tracing_pointers[i, j] = max(zip(possible_scores, pointer_options))

    #Find the optimal path to determine the alignment
    alignment = deque()

    i, j = len_x - 1, len_y - 1

    while i >= 0 and j >= 0:
        trace_dir = tracing_pointers[i, j]

        if trace_dir == DIAG:
            alignment_element = (i, j)
        elif trace_dir == LEFT:
            alignment_element = (i, None)
        elif trace_dir == UP:
            alignment_element = (None, j)
        
        alignment.appendleft(alignment_element)

        i = i + trace_dir[0]
        j = j + trace_dir[1]

    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1
    
    return list(alignment)



def print_alignment(x, y, alignment):
    line = ''

    for a in alignment:
        if a[0] == None:
            line = line + '-' + '\t'
        else:
            line = line + x[a[0]][0] + ':' + x[a[0]][1] + '\t'
    
    print(line)

    line = ''
    for a in alignment:
        if a[1] == None:
            line = line + '-' + '\t'
        else:
            line = line + y[a[1]][0] + ':' + y[a[1]][1] + '\t'
    
    print(line)

def print_rhelix_aln(pairs1, pairs2, aln):
    pos = 0
    for a in aln:
        a1 = a[0]
        a2 = a[1] 
        to_print = str(pos) + ' '
        if a1 == None:
            to_print += ' - '
        else:
            to_print += str(pairs1[a1])
        
        to_print += '     ---------     '

        if a2 == None:
            to_print += ' - '
        else:
            to_print += str(pairs2[a2])

        print(to_print)
        pos += 1




hp_one = [
    ('A', 'U'),
    ('G', 'C'),
    ('U', 'A'),
    ('A', 'U'),
    ('-', 'U'),
    ('G', 'U'),
    ('-', 'U'),
    ('-', 'A'),
    ('-', 'A'),
    ('-', 'G'),
    ('-', 'C')
]
hp_two = [
    ('A', 'U'),
    ('U', 'A'),
    ('A', 'U'),
    ('-', 'U'),
    ('G', 'U'),
    ('-', 'U'),
    ('-', 'A'),
    ('-', 'A'),
    ('-', 'C')
]


hp_one = [
    ('U', 'A'),
    ('C', 'G'),
    ('C', 'G'),
    ('C', '-'),
    ('C', 'G')
]
hp_two = [
    ('U', 'A'),
    ('G', 'C'),
    ('-', 'U'),
    ('U', 'A'),
    ('G', 'C')
]
#print(rhelix_align_verbose(hp_one, hp_two))
print_rhelix_aln(hp_one, hp_two, rhelix_align_verbose(hp_one, hp_two)[0])