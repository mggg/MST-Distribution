"""
This file contains sloppy, not cleaned up at all, code written by Jamie Tucker-Foltz,
as part of a collaboration with Eric Babson, Moon Duchin, Annina Iseli, Pietro Poggi-Corradini,
and Dylan Thurston, to compute the dimension of P_n.
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
    #print(f"Num constraints: {len(to_return)}")
    return to_return


def rank(J, prime=None, clear_denoms=False):
    if clear_denoms:
        cf = 1
        for row in J:
            for entry in row:
                cf = math.lcm(cf, entry.denominator)
        print(cf)
        new_J = []
        for row in J:
            new_row = []
            for entry in row:
                new_entry = entry*cf
                if new_entry.denominator != 1:
                    raise Exception(f"Denominator is {new_entry.denominator}")
                new_row.append(new_entry.numerator)
            new_J.append(new_row)
        J = new_J
    #print("A")
    M = [[x % prime for x in row] for row in J] if prime else \
        [[Fraction(x) for x in row] for row in J]
    #print(len(M))
    #print(len(M[0]))
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
                    #print(f"Multiplier: {multiplier}")
                    for col in range(pivot_col, num_cols):
                        M[row_to_knock_out][col] -= (multiplier*M[pivot_row][col])
                        if prime:
                            M[row_to_knock_out][col] = (M[row_to_knock_out][col] + prime) % prime
            break
        if found:
            possible_pivot_rows.remove(pivot_row)
    #print("B")
    return rank#, M


def word(s):
    return [ord(c) - 97 for c in s]


def print_nicely(a):
    for row in a:
        print(" ".join([str(x) for x in row]))


def jacobian(w):
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
        #print(m)
        to_return.append([m[p] if p in m else 0 for p in full_length_permutations])
    return to_return


def ordered_correctly(cycle, permutation):
    raise Exception("No longer works - use is_subseq instead.")
    for i in range(len(cycle) - 1):
        if permutation[cycle[i]] > permutation[cycle[i + 1]]:
            return False
    return True


def all_ordered_correctly(cycles, permutation):
    for cycle in cycles:
        # These next three lines demonstrate difference between the two functions - T/F's agree.
        #oc = ordered_correctly(cycle, permutation)
        #iss = is_subseq(cycle, inverse(permutation))
        #print(f"{cycle} {permutation} {oc} {iss}")
        if not is_subseq(cycle, permutation):  # Used to call ordered_correctly instead.
            return False
    return True


def inverse(permutation):
    n = len(permutation)
    to_return = [None for _ in range(n)]
    for i in range(n):
        to_return[permutation[i]] = i
    return to_return


def compute_upper_bound_entry(cycles, permutation):
    n = len(permutation)
    common_factor = Fraction(1)
    sum_of_derivative_terms = Fraction(0)
    for cycle in cycles:
        #print(f"Cycle: {cycle}")
        cycle_length_factorial = Fraction(math.factorial(len(cycle)))
        common_factor /= cycle_length_factorial
        if is_subseq(cycle, permutation):
            #print("Ordered correctly")
            sum_of_derivative_terms += cycle_length_factorial
        #print(f"Sum: {sum_of_derivative_terms}")
        #print(f"CF: {common_factor}")
    single_term_contribution = Fraction(1) if all_ordered_correctly(cycles, permutation) else 0
    return sum_of_derivative_terms*common_factor - single_term_contribution


def upper_bound_matrix(constraints, permutations):
    return [[Fraction(1) for _ in range(len(permutations))]] + [[compute_upper_bound_entry(cycles, permutation) for permutation in permutations] for cycles in constraints]


def is_subseq(x, y):
    it = iter(y)
    return all(c in it for c in x)


def generate_dylan_tangent_basis(n, cap=None, top=None):
    if cap is None:
        cap = n
    if top is None:
        top = 2
    full_length_permutations = list(itertools.permutations(range(n)))
    M = []
    for subset_size in range(top, cap + 1):
        for subset in itertools.combinations(range(n), subset_size):
            first = subset[0]
            for rest in itertools.permutations(subset[1:]):
                positive, negative = commutator_sequence(first, rest)
                row = [0 for _ in full_length_permutations]
                for j in range(len(full_length_permutations)):
                    p_2 = full_length_permutations[j]
                    for p_1 in positive:
                        if is_subseq(p_1, p_2):
                            row[j] += 1
                    for p_1 in negative:
                        if is_subseq(p_1, p_2):
                            row[j] -= 1
                M.append(row)
    return M


# From ChatGPT
def rotate_cycle(cycle):
    min_index = cycle.index(min(cycle))
    rotated_cycle = cycle[min_index:] + cycle[:min_index]
    return rotated_cycle


def rotate_cycles(cycles):
    return [rotate_cycle(cycle) for cycle in cycles]


def normalize_tangent_vector(v):
    average = Fraction(sum(v), len(v))
    return [x - average for x in v]


def normalize_tangent_vectors(J):
    return [normalize_tangent_vector(v) for v in J]


def generate_dylan_cotangent_basis(n):
    full_length_permutations = list(itertools.permutations(range(n)))
    M = [[1 for _ in full_length_permutations]]
    cycles_list = [rotate_cycles(cycles) for cycles in get_all_multi_cycle_permutations(n)]
    #print(cycles_list)
    for cycles in cycles_list:
        print(cycles)
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


def generate_jamie_cotangent_basis(n):
    permutations = list(itertools.permutations(range(n)))
    constraints = get_all_multi_cycle_permutations(n)
    return upper_bound_matrix(constraints, permutations)


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


def find_dylans_cotangent(n, v):
    full_length_permutations = list(itertools.permutations(range(n)))
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
        if row == v:
            return cycles
    return None


def find_dylans_tangent(n, v, cap=None):
    if cap is None:
        cap = n
    full_length_permutations = list(itertools.permutations(range(n)))
    for subset_size in range(2, cap + 1):
        for subset in itertools.combinations(range(n), subset_size):
            first = subset[0]
            for rest in itertools.permutations(subset[1:]):
                positive, negative = commutator_sequence(first, rest)
                row = [0 for _ in full_length_permutations]
                for j in range(len(full_length_permutations)):
                    p_2 = full_length_permutations[j]
                    for p_1 in positive:
                        if is_subseq(p_1, p_2):
                            row[j] += 1
                    for p_1 in negative:
                        if is_subseq(p_1, p_2):
                            row[j] -= 1
                if row == v:
                    return first, rest
    return None


def dot(v_1, v_2):
    return sum(v_1_i * v_2_i for v_1_i, v_2_i in zip(v_1, v_2))


test_num = 0
if len(sys.argv) < 2:
    print("USAGE: 'python3 dimension_bounds.py test_num'")
else:
    test_num = int(sys.argv[1])


if test_num == 1:
    J = [
        [1, 2, 3],
        [4, 5, 6]
    ]
    print(rank(J))
    J = [
        [1, 2, 3],
        [4, 8, 12]
    ]
    print(rank(J))
elif test_num == 2:
    J = [
        [1, 2, 3],
        [2, 3, 4],
        [3, 4, 5]
    ]
    print(rank(J))
elif test_num == 3:  # Rank 2 as we expect.
    J = jacobian(word("aba"))
    print_nicely(J)
    print(rank(J))
elif test_num == 4:  # Rank 5, we expected 6. Prints out matrix from pic (rows in diff order).
    J = jacobian(word("abcabcab"))
    print_nicely(J)
    print(rank(J))
elif test_num == 5:  # This is longer than a universal word and the rank is still 5.
    J = jacobian(word("abcabcabcabcabcabcabcabcabcabcabcabcabcabc"))
    print(rank(J))
elif test_num == 6:  # Still 5.
    n = 3
    for word_length in range(n, 100):
        w = []
        letter = 0
        for _ in range(word_length):
            w.append(letter)
            letter = (letter + 1) % n
        J = jacobian(w)
        print(rank(J))
elif test_num == 7:  # Try different n. Dimensions max out at 2, 5, 12, 20, 30, 42, 56.
    n = 8
    for word_length in range(n, 100):
        w = []
        letter = 0
        for _ in range(word_length):
            w.append(letter)
            letter = (letter + 1) % n
        J = jacobian(w)
        print(rank(J))
elif test_num == 8:  # Exhaustive search over better words: 
    for word_length in range(8, 12):
        ranks = [0 for _ in range(7)]
        print(f"Word length: {word_length}")
        for w in itertools.product(range(3), repeat=word_length):
            J = jacobian(w)
            if len(set(w)) < 3:
                continue
            nontrivial = True
            for j in range(6):
                all_zero = True
                for i in range(word_length):
                    if J[i][j] > 0:
                        all_zero = False
                        break
                if all_zero:
                    nontrivial = False
                    break
            if not nontrivial:
                continue
            for i in range(word_length - 1):
                if w[i] == w[i + 1]:
                    nontrivial = False
                    break
            if nontrivial:
                ranks[rank(J)] += 1
        print(ranks)
elif test_num == 9:  # Switching last two letters from blackboard example gets full rank.
    J = jacobian(word("abcabcba"))
    print_nicely(J)
    print(rank(J))
elif test_num == 10:  # Randomized search for P4.
    for word_length in range(1, 31):
        ranks = {i: 0 for i in range(25)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(4)) for _ in range(word_length)] for _ in range(100)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J)] += 1
        print(ranks)
elif test_num == 11:  # Randomized search for P4 just with long words.
    for word_length in range(30, 32):
        ranks = {i: 0 for i in range(25)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(4)) for _ in range(word_length)] for _ in range(10000)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J)] += 1
        print(ranks)
elif test_num == 12:  # Randomized search for P5 just with long words.
    for word_length in range(120, 141):
        ranks = {i: 0 for i in range(121)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(5)) for _ in range(word_length)] for _ in range(20)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J)] += 1
        for i in range(121):
            r = ranks[i]
            if r > 0:
                print(f"{i}: {r}")
elif test_num == 13:  # Randomized search for P6 just with long words, best known lower bound.
    for word_length in range(600, 700):
        ranks = {i: 0 for i in range(700)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(6)) for _ in range(word_length)] for _ in range(1)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J, 997)] += 1
        for i in range(700):
            r = ranks[i]
            if r > 0:
                print(f"{i}: {r}")
elif test_num == 14:  # Like test 12 but modulo 997.
    for word_length in range(120, 141):
        ranks = {i: 0 for i in range(121)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(5)) for _ in range(word_length)] for _ in range(20)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J, 997)] += 1
        for i in range(121):
            r = ranks[i]
            if r > 0:
                print(f"{i}: {r}")
elif test_num == 15:  # Like test 14 but with larger words.
    for word_length in range(150, 300):
        ranks = {i: 0 for i in range(121)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(5)) for _ in range(word_length)] for _ in range(1)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J, 997)] += 1
        for i in range(121):
            r = ranks[i]
            if r > 0:
                print(f"{i}: {r}")
elif test_num == 16:  # Randomized search for P6 just with long words.
    for word_length in range(10, 2501, 10):
        ranks = {i: 0 for i in range(2500)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(7)) for _ in range(word_length)] for _ in range(1)]:
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J, 997)] += 1
        for i in range(2500):
            r = ranks[i]
            if r > 0:
                print(f"{i}: {r}")
elif test_num == 17:
    print(compute_upper_bound_entry([(0, 1), (2, 3)], (0, 1, 2, 3)))
elif test_num == 18:
    J = upper_bound_matrix([[(0, 1), (2, 3)], [(1, 0), (2, 3)]], list(itertools.permutations(range(4))))
    for row in J:
        print(list(map(str, row)))
elif test_num == 19:
    J = upper_bound_matrix([[(0, 1), (2, 3)], [(1, 0), (2, 3)]], list(itertools.permutations(range(4))))
    print(rank(J))
elif test_num == 20:
    J = upper_bound_matrix([[(0, 1), (2, 3)], [(1, 0), (2, 3)], [(0, 1), (3, 2)], [(1, 0), (3, 2)]], list(itertools.permutations(range(4))))
    print(rank(J))
elif test_num == 21:
    J = upper_bound_matrix([[(0, 1), (2, 3)], [(1, 0), (2, 3)], [(0, 1), (3, 2)], [(1, 0), (3, 2)], [(0, 2), (1, 3)]], list(itertools.permutations(range(4))))
    print(rank(J))
elif test_num == 22:
    for permutation in itertools.permutations(range(4)):
        print(to_cycles(permutation))
elif test_num == 23:  # This is the main test that gave the upper bound.
    for n in range(3, 7):
        print(f"n: {n}")
        J = upper_bound_matrix(get_all_multi_cycle_permutations(n), list(itertools.permutations(range(n))))
        #print(J)
        print(f"Rank: {rank(J)}")
        print("")
elif test_num == 24:
    n = 5
    cycles = [[(a, b, c), (d, e)] for a, b, c, d, e in itertools.permutations(range(n))]
    J = upper_bound_matrix(cycles, list(itertools.permutations(range(n))))
    print(f"Rank: {rank(J)}")
elif test_num == 25:
    n = 5
    J = upper_bound_matrix(get_all_multi_cycle_permutations(n), list(itertools.permutations(range(n))))
    for row in J:
        print("".join([(" " if v == 0 else "X") for v in row]))
elif test_num == 26:  # Upper bound attempt for P7, want to see 2675. DONE IN 1 HOUR.
    for n in range(3, 8):
        print(f"n: {n}")
        J = upper_bound_matrix(get_all_multi_cycle_permutations(n), list(itertools.permutations(range(n))))
        #print(J)
        print(f"Rank: {rank(J, 997, clear_denoms=True)}")
        print("")
elif test_num == 27:  # Lower bound attempt for P7, want to see 2366. DONE IN 6 HOURS ON FASRC.
    for word_length in range(3000, 3001):
        ranks = {i: 0 for i in range(10000)}
        print(f"Word length: {word_length}")
        for w in [[random.choice(range(7)) for _ in range(word_length)] for _ in range(1)]:
            if len(set(w)) < 4:
                continue
            J = jacobian(w)
            nontrivial = True
            ranks[rank(J, 997)] += 1
        for i in range(10000):
            r = ranks[i]
            if r > 0:
                print(f"{i}: {r}")
elif test_num == 28:  # Test Dylan's April 10 conjecture.
    targets = [None, None, 2, 6, 21, 85, 410, 2366]  # One more than dimension of P_n
    for n in range(3, 5):#range(2, 7):
        target_dimension = targets[n]
        for word_length in range(1, 999):
            w = [random.choice(range(n)) for _ in range(word_length)]
            J_1 = jacobian(w)
            rank_1 = rank(J_1)
            if rank_1 == target_dimension:
                print(f"For n = {n}, first universal word found was of length {word_length}.")
                print(f"Universal word: {''.join(map(str, w))}")
                print(f"Column indices:")
                print_nicely([list(itertools.permutations(range(n)))])
                print("Combined matrix:")
                J_2 = generate_dylan_tangent_basis(n)
                rank_2 = rank(J_2)
                J_3 = J_1 + J_2
                rank_3 = rank(J_3)
                print_nicely(J_3)
                print(f"Rank of universal word (upper part of matrix): {rank_1}")
                print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
                print(f"Combined rank (whole matrix): {rank_3}\n")
                break
elif test_num == 29:
    print_nicely(generate_dylan_tangent_basis(3))
    print("")
    print_nicely(generate_dylan_tangent_basis(4))
elif test_num == 30:
    print(commutator_sequence(0, [1, 2]))
elif test_num == 31:  # Test Dylan's April 10 conjecture, actually using uniform distribution.
    targets = [None, None, 2, 6, 21, 85, 410, 2366]  # One more than dimension of P_n
    words = [None, None, None, "cabacabac", "dcabaccccabacdddcabaccccabacdddcabaccccabacd"]
    for n in range(3, 5):#range(2, 7):
        target_dimension = targets[n]
        w = word(words[n])
        J_1 = jacobian(w)
        rank_1 = rank(J_1)
        if True:#rank_1 == target_dimension:
            print(f"For n = {n}:.")
            print(f"Universal word: {''.join(map(str, w))}")
            print(f"Column indices:")
            print_nicely([list(itertools.permutations(range(n)))])
            print("Combined matrix:")
            J_2 = generate_dylan_tangent_basis(n)
            rank_2 = rank(J_2)
            J_3 = J_1 + J_2
            rank_3 = rank(J_3)
            print_nicely(J_3)
            print(f"Rank of universal word (upper part of matrix): {rank_1}")
            print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
            print(f"Combined rank (whole matrix): {rank_3}\n")
elif test_num == 32:  # Like test 31 but with equal num occurrences of each letter.
    targets = [None, None, 2, 6, 21, 85, 410, 2366]  # One more than dimension of P_n
    words = [None, None, None, "ccccaaabbbbbbaaaccccaaabbbbbbaaacccc", "dddddddddccccaaaaaabbbbbbbbbbbbaaaaaaccccccccccccccccaaaaaabbbbbbbbbbbbaaaaaaccccdddddddddddddddddddddddddddccccaaaaaabbbbbbbbbbbbaaaaaaccccccccccccccccaaaaaabbbbbbbbbbbbaaaaaaccccdddddddddddddddddddddddddddccccaaaaaabbbbbbbbbbbbaaaaaaccccccccccccccccaaaaaabbbbbbbbbbbbaaaaaaccccddddddddd"]
    for n in range(3, 5):#range(2, 7):
        target_dimension = targets[n]
        w = word(words[n])
        J_1 = jacobian(w)
        rank_1 = rank(J_1)
        if True:#rank_1 == target_dimension:
            print(f"For n = {n}:.")
            print(f"Universal word: {''.join(map(str, w))}")
            print(f"Column indices:")
            print_nicely([list(itertools.permutations(range(n)))])
            print("Combined matrix:")
            J_2 = generate_dylan_tangent_basis(n)
            rank_2 = rank(J_2)
            J_3 = J_1 + J_2
            rank_3 = rank(J_3)
            print_nicely(J_3)
            print(f"Rank of universal word (upper part of matrix): {rank_1}")
            print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
            print(f"Combined rank (whole matrix): {rank_3}\n")
elif test_num == 33:  # Like test 31 but just with non-nested commutators.
    targets = [None, None, 2, 6, 21, 85, 410, 2366]  # One more than dimension of P_n
    words = [None, None, None, "cabacabac", "dcabaccccabacdddcabaccccabacdddcabaccccabacd"]
    for n in range(3, 5):#range(2, 7):
        target_dimension = targets[n]
        w = word(words[n])
        J_1 = jacobian(w)
        rank_1 = rank(J_1)
        if True:#rank_1 == target_dimension:
            print(f"For n = {n}:.")
            print(f"Universal word: {''.join(map(str, w))}")
            print(f"Column indices:")
            print_nicely([list(itertools.permutations(range(n)))])
            print("Combined matrix:")
            J_2 = generate_dylan_tangent_basis(n, cap=2)
            rank_2 = rank(J_2)
            J_3 = J_1 + J_2
            rank_3 = rank(J_3)
            print_nicely(J_3)
            print(f"Rank of universal word (upper part of matrix): {rank_1}")
            print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
            print(f"Combined rank (whole matrix): {rank_3}\n")
elif test_num == 34:  # Test Dylan's April 19 conjecture.
    for n in range(4, 6):
        J_1 = upper_bound_matrix(get_all_multi_cycle_permutations(n), list(itertools.permutations(range(n))))
        print(f"For n = {n}:")
        print(f"Column indices:")
        print_nicely([list(itertools.permutations(range(n)))])
        J_2 = generate_dylan_cotangent_basis(n)
        #print_nicely(J_1)
        #print(J_1)
        #print_nicely(J_2)
        #print(J_2)
        print("Combined matrix:")
        J_3 = J_1 + J_2
        print_nicely(J_3)
        rank_1 = rank(J_1)
        rank_2 = rank(J_2)
        rank_3 = rank(J_3)
        print(f"Rank of original construction (upper part of matrix): {rank_1}")
        print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
        print(f"Combined rank (whole matrix): {rank_3}\n")
elif test_num == 35:
    print(get_all_multi_cycle_permutations(5))
elif test_num == 36:
    print(rotate_cycle((3,5,1,2,4)))
elif test_num == 37:  # My tangent and cotangent are NOT orthogonal.
    n = 4
    J = jacobian(word("dcabaccccabacdddcabaccccabacdddcabaccccabacd"))
    print_nicely(J)
    M = upper_bound_matrix(get_all_multi_cycle_permutations(n), list(itertools.permutations(range(n))))
    print_nicely(M)
    for v1 in M:
        for v2 in J:
            p = 0
            for v1i, v2i in zip(v1, v2):
                p += v1i*v2i
            print(p)
elif test_num == 38:  # Dylans tangent and cotangent ARE orthogonal for n = 4, but...
    n = 4
    J = generate_dylan_tangent_basis(n, cap=2)
    print_nicely(J)
    print("")
    M = generate_dylan_cotangent_basis(n)
    print_nicely(M)
    for v1 in J:
        for v2 in M:
            p = 0
            for v1i, v2i in zip(v1, v2):
                p += v1i*v2i
            print(p)
elif test_num == 39:  # ... not for n = 5.
    n = 5
    J = generate_dylan_tangent_basis(n, cap=2)
    print_nicely(J)
    print("")
    M = generate_dylan_cotangent_basis(n)
    print_nicely(M)
    for v1 in J:
        for v2 in M:
            p = 0
            for v1i, v2i in zip(v1, v2):
                p += v1i*v2i
            print(p)
elif test_num == 40:
    n = 5
    J = generate_dylan_tangent_basis(n, cap=2)
    M = generate_dylan_cotangent_basis(n)
    for v1 in J:
        for v2 in M:
            p = 0
            for v1i, v2i in zip(v1, v2):
                p += v1i*v2i
            if p != 0:
                print(p)
                print(v1)
                print(v2)
                raise Exception("STOP")
elif test_num == 41:  # Track down one of the non-orthogonal pairs from previous test.
    n = 5
    print(find_dylans_tangent(n, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1], cap=2))
    print(find_dylans_cotangent(n, [1, 0, 1, 1, 0, 0, 1, 0, -1, 0, -1, 1, -1, -1, -1, 0, 0, 0, -1, -1, -1, 1, 1, 1, 1, 0, 1, 1, 0, 0, -1, 0, -1, -1, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 1, -1, -1, 0, -1, -1, 0, 0, 1, 0, 1, 1, 0, 0, 1, -1, 1, 1, -1, -1, -1, -1, -1, 0, 0, 0, -1, -1, 1, 1, -1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1]))
elif test_num == 42:  # My tangent vectors are orthogonal to Dylan's cotangent vectors.
    n = 4
    J = jacobian(word("dcabaccccabacdddcabaccccabacdddcabaccccabacd"))
    print_nicely(J)
    print("")
    M = generate_dylan_cotangent_basis(n)
    print_nicely(M)
    for v1 in J:
        for v2 in M:
            p = 0
            for v1i, v2i in zip(v1, v2):
                p += v1i*v2i
            print(p)
elif test_num == 43:
    cycles = [(0, 1), (2, 3)]
    permutation = (0, 3, 1, 2)
    print(compute_upper_bound_entry(cycles, permutation))
elif test_num == 44:
    print(is_subseq((2, 3), (0, 3, 1, 2)))
elif test_num == 45:  # First test after resolving bug! Similar to test 18.
    M = generate_jamie_cotangent_basis(4)
    for row in M:
        print(list(map(str, row)))
elif test_num == 46:  # My cotangent space is the same as Dylan's.
    for n in range(3, 7):
        print(f"n = {n}:")
        M_1 = generate_jamie_cotangent_basis(n)
        M_2 = generate_dylan_cotangent_basis(n)
        #print_nicely(M_1 + M_2)
        print(rank(M_1))
        print(rank(M_2))
        print(rank(M_1 + M_2))
elif test_num == 47:  # Verify spaces are orthogonal (they're not, but are in next test).
    n = 4
    J = jacobian(word("dcabaccccabacdddcabaccccabacdddcabaccccabacd"))
    for M in [generate_jamie_cotangent_basis(n), generate_dylan_cotangent_basis(n)]:
        print_nicely(M)
        for J_i in J:
            for M_i in M:
                print(dot(J_i, M_i))
elif test_num == 48:  # Same as test 47 but normalizing so that it's in P_n and not P_n tilde.
    n = 4
    J = normalize_tangent_vectors(jacobian(word("dcabaccccabacdddcabaccccabacdddcabaccccabacd")))
    for M in [generate_jamie_cotangent_basis(n), generate_dylan_cotangent_basis(n)]:
        print_nicely(M)
        for J_i in J:
            for M_i in M:
                print(dot(J_i, M_i))
elif test_num == 49:  # Like test 31, but normalizing.
    targets = [None, None, 2, 6, 21, 85, 410, 2366]  # One more than dimension of P_n
    words = [None, None, None, "cabacabac", "dcabaccccabacdddcabaccccabacdddcabaccccabacd"]
    for n in range(3, 5):#range(2, 7):
        target_dimension = targets[n]
        w = word(words[n])
        J_1 = normalize_tangent_vectors(jacobian(w))
        rank_1 = rank(J_1)
        if True:#rank_1 == target_dimension:
            print(f"For n = {n}:.")
            print(f"Universal word: {''.join(map(str, w))}")
            print(f"Column indices:")
            print_nicely([list(itertools.permutations(range(n)))])
            print("Combined matrix:")
            J_2 = generate_dylan_tangent_basis(n)
            rank_2 = rank(J_2)
            J_3 = J_1 + J_2
            rank_3 = rank(J_3)
            print_nicely(J_3)
            print(f"Rank of universal word (upper part of matrix): {rank_1}")
            print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
            print(f"Combined rank (whole matrix): {rank_3}\n")
elif test_num == 50:  # Figure out that it's only the middle ones that are the troublemakers.
    targets = [None, None, 2, 6, 21, 85, 410, 2366]  # One more than dimension of P_n
    words = [None, None, None, "cabacabac", "dcabaccccabacdddcabaccccabacdddcabaccccabacd"]
    for n in range(4, 5):#range(2, 7):
        target_dimension = targets[n]
        w = word(words[n])
        J_1 = normalize_tangent_vectors(jacobian(w))
        rank_1 = rank(J_1)
        if True:#rank_1 == target_dimension:
            print(f"For n = {n}:.")
            print(f"Universal word: {''.join(map(str, w))}")
            print(f"Column indices:")
            print_nicely([list(itertools.permutations(range(n)))])
            print("Combined matrix:")
            J_2 = generate_dylan_tangent_basis(n, cap=3, top=3)
            rank_2 = rank(J_2)
            J_3 = J_1 + J_2
            rank_3 = rank(J_3)
            print_nicely(J_3)
            print(f"Rank of universal word (upper part of matrix): {rank_1}")
            print(f"Rank of Dylan's construction (lower part of matrix): {rank_2}")
            print(f"Combined rank (whole matrix): {rank_3}\n")
elif test_num == 51:
    n = 5
    for size in range(2, n + 1):
        print(f"size: {size}")
        for v in generate_dylan_tangent_basis(n, cap=size, top=size):
            orth = True
            for dual in generate_dylan_cotangent_basis(n):
                if sum(x*y for x, y in zip(v, dual)) != 0:
                    orth = False
                    break
            print(f"[{orth}]")
elif test_num == 52:  # Following May 13 meeting, see if [a, [b, cde]] is in tangent space.
    n = 5
    v = []
    for p in itertools.permutations(range(n)):
        if p[0] == 0 and p[1] == 1:
            v.append(1)
        elif p[3] == 1 and p[4] == 0:
            v.append(1)
        elif p[0] == 0 and p[4] == 1:
            v.append(-1)
        elif p[0] == 1 and p[4] == 0:
            v.append(-1)
        else:
            v.append(0)
    print(f"v: {v}")
    for dual in generate_dylan_cotangent_basis(n):
        print(f"\nCotangent vector: {dual}")
        print(dot(dual, v))



