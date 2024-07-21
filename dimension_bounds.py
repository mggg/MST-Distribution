"""
This file contains code written by Jamie Tucker-Foltz,
as part of a collaboration with Eric Babson, Moon Duchin, Annina Iseli, Pietro Poggi-Corradini,
and Dylan Thurston, to compute the dimension of P_n (conjectured OEIS sequence A006231).
"""

import itertools
from fractions import Fraction
import random
import math
import sys


# From https://stackoverflow.com/questions/36683413/permutations-producing-cycle-notation
def to_cycles(perm):
    pi = {i: perm[i] for i in range(len(perm))}
    cycles = []
    while pi:
        elem0 = min(list(pi.keys())) # Alternative doesn't work: random.choice(list(pi.keys()))
        this_elem = pi[elem0]
        next_item = pi[this_elem]
        cycle = []
        while True:
            cycle.append(this_elem)
            del pi[this_elem]
            this_elem = next_item
            if next_item in pi:
                next_item = pi[next_item]
            else:
                break
        cycles.append(cycle)
    final_cycles = []
    for cycle in cycles:
        if len(cycle) > 1:
            final_cycles.append(tuple(cycle))
    return final_cycles


def get_all_multi_cycle_permutations(n):
    to_return = []
    for permutation in itertools.permutations(range(n)):
        cycles = to_cycles(permutation)
        if len(cycles) >= 2:
            to_return.append(cycles)
    return to_return


def is_subseq(x, y):
    it = iter(y)
    return all(c in it for c in x)


def rotate_cycle(cycle):
    min_index = cycle.index(min(cycle))
    rotated_cycle = cycle[min_index:] + cycle[:min_index]
    return rotated_cycle


def rotate_cycles(cycles):
    return [rotate_cycle(cycle) for cycle in cycles]


def commutator_sequence(first, rest):
    positive, negative = [(first,)], []
    for r in rest:
        new_char = (r,)
        new_positive = []
        new_negative = []
        for p in positive:
            new_positive.append(p + new_char)
            new_negative.append(new_char + p)
        for p in negative:
            new_negative.append(p + new_char)
            new_positive.append(new_char + p)
        positive = new_positive
        negative = new_negative
    return positive, negative


def generate_cotangent_basis(n):
    full_length_permutations = list(itertools.permutations(range(n)))
    M = [[1 for _ in full_length_permutations]]
    cycles_list = [rotate_cycles(cycles) for cycles in get_all_multi_cycle_permutations(n)]
    for cycles in cycles_list:
        positive_list = []
        negative_list = []
        for cycle in cycles:
            first = cycle[0]
            rest = cycle[1:]
            positive, negative = commutator_sequence(first, rest)
            positive_list.append(positive)
            negative_list.append(negative)
        row = []
        for p_2 in full_length_permutations:
            coefficients_multiplied = 1
            for positive, negative in zip(positive_list, negative_list):
                coefficient = 0
                for p_1 in positive:
                    if is_subseq(p_1, p_2):
                        coefficient += 1
                for p_1 in negative:
                    if is_subseq(p_1, p_2):
                        coefficient -= 1
                coefficients_multiplied *= coefficient
            row.append(coefficients_multiplied)
        M.append(row)
    return M


def word(s):
    return [ord(c) - 97 for c in s]


def draw_matrix(w):
    word_length = len(w)
    n = len(set(w))
    full_length_permutations = list(itertools.permutations(range(n)))
    to_return = []
    for index_must_include in range(word_length):
        letter_must_include = w[index_must_include]
        m = {(): 1}
        for i in range(word_length):
            next_letter = w[i]
            if i == index_must_include:
                m = {p + (letter_must_include,): m[p] for p in m}
            elif next_letter != letter_must_include:
                new_m = {}
                for p in m:
                    if p in new_m:
                        new_m[p] += m[p]
                    else:
                        new_m[p] = m[p]
                    if next_letter not in p:
                        new_p = p + (next_letter,)
                        if new_p in new_m:
                            new_m[new_p] += m[p]
                        else:
                            new_m[new_p] = m[p]
                m = new_m
        to_return.append([m[p] if p in m else 0 for p in full_length_permutations])
    return to_return


def rank(M, prime=None, clear_denoms=False):
    if clear_denoms:
        cf = 1
        for row in M:
            for entry in row:
                cf = math.lcm(cf, entry.denominator)
        print(cf)
        new_M = []
        for row in M:
            new_row = []
            for entry in row:
                new_entry = entry*cf
                if new_entry.denominator != 1:
                    raise Exception(f"Denominator is {new_entry.denominator}")
                new_row.append(new_entry.numerator)
            new_M.append(new_row)
        M = new_M
    M = [[x % prime for x in row] for row in M] if prime else \
        [[Fraction(x) for x in row] for row in M]
    num_rows = len(M)
    num_cols = len(M[0])
    if num_rows < num_cols:  # Transpose so that matrix is tall and skinny.
        M = [[M[j][i] for j in range(num_rows)] for i in range(num_cols)]
        temp = num_cols
        num_cols = num_rows
        num_rows = temp
    possible_pivot_rows = set(range(num_rows))
    rank = 0
    for pivot_col in range(num_cols):
        found = False
        for pivot_row in possible_pivot_rows:
            pivot = M[pivot_row][pivot_col]
            if pivot == 0:
                continue
            found = True
            rank += 1
            for row_to_knock_out in possible_pivot_rows:
                if row_to_knock_out != pivot_row:
                    if prime:
                        multiplier = M[row_to_knock_out][pivot_col]*\
                                     pow(M[pivot_row][pivot_col], -1, prime) % prime
                    else:
                        multiplier = M[row_to_knock_out][pivot_col]/M[pivot_row][pivot_col]
                    for col in range(pivot_col, num_cols):
                        M[row_to_knock_out][col] -= (multiplier*M[pivot_row][col])
                        if prime:
                            M[row_to_knock_out][col] = (M[row_to_knock_out][col] + prime) % prime
            break
        if found:
            possible_pivot_rows.remove(pivot_row)
    return rank


def print_nicely(a):
    for row in a:
        print(" ".join([str(x) for x in row]))


test_num = 0
if len(sys.argv) < 2:
    print("USAGE: 'python3 dimension_bounds.py test_num'")
else:
    test_num = int(sys.argv[1])


if test_num == 1:  # Verify rank computation works.
    M = [
        [1, 2, 3],
        [4, 5, 6]
    ]
    print(f"Rank: {rank(M)}")
    M = [
        [1, 2, 3],
        [4, 8, 12]
    ]
    print(f"Rank: {rank(M)}")
elif test_num == 2:  # Rank 2 as we expect.
    M = draw_matrix(word("aba"))
    print_nicely(M)
    print(f"Rank: {rank(M)}")
elif test_num == 3:  # This is longer than a universal word but the rank is still 5 - bad pattern.
    M = draw_matrix(word("abcabcabcabcabcabcabcabcabcabcabcabcabcabc"))
    print_nicely(M)
    print(f"Rank: {rank(M)}")
elif test_num == 4:  # Short word that gives draw matrix M of full rank, example appears in paper.
    M = draw_matrix(word("abcabcba"))
    print_nicely(M)
    print(f"Rank: {rank(M)}")
elif test_num == 5:  # Lower bounds for P3, P4, P5, P6, want to see 6, 21, 85, 410.
    for n in range(3, 7):
        print(f"n: {n}")
        word_length = 600
        print(f"Word length: {word_length}")
        w = [random.choice(range(n)) for _ in range(word_length)]
        M = draw_matrix(w)
        print(f"Rank: {rank(M, 997)}\n")
elif test_num == 6:  # Upper bounds for P3, P4, P5, P6, want to see 1, 4, 36, 311.
    for n in range(3, 7):
        print(f"n: {n}")
        M = generate_cotangent_basis(n)
        print(f"Rank: {rank(M)}")
        print("")
elif test_num == 7:  # Lower bound for P7, want to see 2366. DONE IN 6 HOURS ON HARVARD CLUSTER.
    n = 7
    print(f"n: {n}")
    word_length = 3000
    print(f"Word length: {word_length}")
    w = [random.choice(range(n)) for _ in range(word_length)]
    M = draw_matrix(w)
    print(f"Rank: {rank(M, 997)}\n")
elif test_num == 8:  # Upper bound for P7, want to see 2675 = 8! - 2365. DONE IN 1 HOUR.
    n = 7
    print(f"n: {n}")
    M = generate_cotangent_basis(n)
    print(f"Rank: {rank(M, 997)}")
    print("")

