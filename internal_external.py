"""
This file contains semi-cleaned-up versions of various code written by Jamie Tucker-Foltz,
as part of a collaboration with Eric Babson, Moon Duchin, Annina Iseli, Pietro Poggi-Corradini,
and Dylan Thurston, to compute and compare MST probabilities of certain trees.
"""

import itertools
import networkx as nx
from fractions import Fraction
import math
import sys


def can_add_edge(tree_edges, island_1, island_2):
    for v_1 in island_1:
        for v_2 in island_2:
            if (v_1, v_2) in tree_edges:
                return True
    return False
        

# Computes the probability of sampling a given tree from a complete graph using the internal formula.
def internal(tree):
    tree_edges_1 = list(tree.edges())
    tree_edges_2 = list(map(lambda e: (e[1], e[0]), tree_edges_1))
    tree_edges = set(tree_edges_1 + tree_edges_2)
    initial_islands = tuple(sorted(frozenset([x]) for x in tree.nodes()))
    num_vertices = len(tree.nodes())
    if len(tree_edges_1) + 1 != num_vertices:
        raise Exception("Not a spanning tree.")
    return internal_helper(tree_edges, initial_islands, {}, int(num_vertices*(num_vertices - 1)/2))


def internal_helper(tree_edges, islands, cache, num_edges_remaining):
    if num_edges_remaining == 0:
        print(f"Hit first base case (shouldn't happen): NER = {num_edges_remaining}.")
        return 1
    elif len(islands) == 2:
        return Fraction(1, num_edges_remaining)
    elif islands in cache:
        return cache[islands]
    else:
        p = 0
        num_islands = len(islands)
        for i in range(num_islands):
            for j in range(i):
                if can_add_edge(tree_edges, islands[i], islands[j]):
                    new_islands_list = []
                    for k in range(num_islands):
                        if k != i and k != j:
                            new_islands_list.append(islands[k])
                    new_islands_list.append(frozenset(list(islands[i]) + list(islands[j])))
                    p += internal_helper(tree_edges, tuple(sorted(new_islands_list)), cache, num_edges_remaining - len(islands[i])*len(islands[j]))
        p /= num_edges_remaining
        cache[islands] = p
        return p


# Computes the probability of sampling a given tree from a graph g using the external formula, where the edges in g missing from the tree are specified.
def external(g, missing_edges):
    test_g = g.copy()
    for e in missing_edges:
        test_g.remove_edge(*e)
    if nx.is_connected(test_g):
        return external_helper(g, missing_edges, {})
    else:
        return 0


def external_helper(g, missing_edges, cache):
    key = tuple(sorted(map(lambda e: tuple(sorted(e)), g.edges())))
    if key in cache:
        return cache[key]
    bridges = set(nx.bridges(g))
    num_bridges = len(bridges)
    num_removable = g.size() - num_bridges
    if len(missing_edges) == 1:
        probability = Fraction(1, num_removable)
        cache[key] = probability
        return probability
    else:
        numerator = 0
        for e in missing_edges:
            new_g = g.copy()
            new_g.remove_edge(*e)
            new_missing_edges = missing_edges - {e}
            numerator += external_helper(new_g, new_missing_edges, cache)
        probability = numerator/num_removable
        cache[key] = probability
        return probability


# Quick way of drawing graphs using wasd as arrow keys to move around a grid and space to draw an edge.
def wasd_code_graph(s):
    x = 0
    y = 0
    last = (x, y)
    m = {last: 0}
    g = nx.Graph()
    g.add_node(0)
    next_id = 1
    for c in s:
        if c == "a":
            x -= 1
        elif c == "d":
            x += 1
        elif c == "s":
            y -= 1
        elif c == "w":
            y += 1
        elif c == " ":
            if (x, y) not in m:
                m[(x, y)] = next_id
                g.add_node(next_id)
                next_id += 1
            g.add_edge(m[last], m[(x, y)])
            last = (x, y)
    return g


def get_path(T, source, target):
    p = nx.shortest_path(T, source=source, target=target)
    to_return = set()
    for i in range(len(p) - 1):
        to_return.add(frozenset((p[i], p[i + 1])))
    return to_return


def find_bijection(T1, T2):
    nodes = list(T1.nodes())
    edges_T1 = list(map(frozenset, T1.edges()))
    edges_T2 = list(map(frozenset, T2.edges()))

    for p in itertools.permutations(edges_T2):
        bijection = {edges_T1[i]: p[i] for i in range(len(p))}
        G = nx.Graph()
        left_nodes = [(False, frozenset((a, b))) for a, b in itertools.product(nodes, nodes) \
                      if a < b and not T1.has_edge(a, b)]
        right_nodes = [(True, frozenset((a, b))) for a, b in itertools.product(nodes, nodes) \
                      if a < b and not T2.has_edge(a, b)]
        G.add_nodes_from(left_nodes)
        G.add_nodes_from(right_nodes)

        for _, edge1 in left_nodes:
            for _, edge2 in right_nodes:
                path1 = get_path(T1, source=min(edge1), target=max(edge1))
                path2 = get_path(T2, source=min(edge2), target=max(edge2))

                if {bijection[edge] for edge in path1}.issubset(path2):
                    G.add_edge((False, edge1), (True, edge2))
        matching = nx.max_weight_matching(G)
        if len(matching) == len(G.nodes())/2:
            return bijection, matching

    return None


test_num = 0
if len(sys.argv) < 2:
    print("USAGE: 'python3 jamie_mst_probabilities.py test_num'")
else:
    test_num = int(sys.argv[1])


if test_num == 1:  # Computes probabilities of sampling paths in K_n; numerators are OEIS sequence A374293.
    for i in range(2, 17):
        tree = nx.path_graph(i)
        p = internal(tree)
        print(f"\nn = {i}")
        print(f"Probability of sampling a specific path from K_n:")
        print(p)
        print(f" = {float(p)}")
        print(f"Numerator (num orderings of edges leading to the path):")
        print(p * math.factorial(int(i*(i - 1)/2)))
elif test_num == 2:  # Shows how to use external formula function.
    g = nx.Graph([(1, 2), (0, 1), (0, 2)])
    missing_edges = {(0, 1)}
    print(external(g, missing_edges))
elif test_num == 3:  # Shows no cycle-increasing bijection exists between either pair of graphs:
    # T1: O
    #     |
    #     O
    #     |
    # O-O-O-O-O
    #
    #
    # T2:
    #
    # O-O-O-O-O-O-O
    #
    T1 = wasd_code_graph("d d d d a a w w ")
    T2 = wasd_code_graph("d d d d d d ")
    print(find_bijection(T1, T2))
    print(find_bijection(T2, T1))

"""
Output of test 1:

n = 2
Probability of sampling a specific path from K_n:
1
 = 1.0
Numerator (num orderings of edges leading to the path):
1

n = 3
Probability of sampling a specific path from K_n:
1/3
 = 0.3333333333333333
Numerator (num orderings of edges leading to the path):
2

n = 4
Probability of sampling a specific path from K_n:
11/180
 = 0.06111111111111111
Numerator (num orderings of edges leading to the path):
44

n = 5
Probability of sampling a specific path from K_n:
113/15120
 = 0.007473544973544973
Numerator (num orderings of edges leading to the path):
27120

n = 6
Probability of sampling a specific path from K_n:
21881/32432400
 = 0.000674664841331508
Numerator (num orderings of edges leading to the path):
882241920

n = 7
Probability of sampling a specific path from K_n:
46253/966984480
 = 4.7832205124946785e-05
Numerator (num orderings of edges leading to the path):
2443792425984000

n = 8
Probability of sampling a specific path from K_n:
131494614599/47359225289376000
 = 2.7765364360488794e-06
Numerator (num orderings of edges leading to the path):
846533597741050576896000

n = 9
Probability of sampling a specific path from K_n:
401720718947/2954952703422720000
 = 1.3594827371743958e-07
Numerator (num orderings of edges leading to the path):
50571850611494440562578575851520000

n = 10
Probability of sampling a specific path from K_n:
1958264695513589873/341074932465086139640320000
 = 5.7414500718666745e-09
Numerator (num orderings of edges leading to the path):
686805008584962439650318114385825747697664000000

n = 11
Probability of sampling a specific path from K_n:
3918430726447198567960531/18414082814237990859865088915712000
 = 2.1279532442514164e-10
Numerator (num orderings of edges leading to the path):
2701735270674169239689693528384644314472371275610193920000000000

n = 12
Probability of sampling a specific path from K_n:
1004211915275585293951070659343/143100424990201202655916245395129395200000
 = 7.017532724618733e-12
Numerator (num orderings of edges leading to the path):
3819958423456547324072333722421751679308286064300212197312630212725309440000000000

n = 13
Probability of sampling a specific path from K_n:
4412801128905638844377098968509009257/21190741193182066471630651785018406937726976000000
 = 2.082419434354386e-13
Numerator (num orderings of edges leading to the path):
2358190320559038013253038734002134501056785955033525100306183708271864557601556171595972608000000000000

n = 14
Probability of sampling a specific path from K_n:
34724813087641197329987933150401131669402059/6186975914245301867997740030004288258847465273436160000000
 = 5.612566392522799e-15
Numerator (num orderings of edges leading to the path):
758819833688728668682269168314076661888759143276671968466480158793861501876004861227804730198555357876721601740800000000000000

n = 15
Probability of sampling a specific path from K_n:
108183610375360277148317376866292834550831527746957417/781245180337980228317513278385849316415941866170148961941898240000000
 = 1.3847587556131633e-16
Numerator (num orderings of edges leading to the path):
149747362926493391018694406416033310179041677974350548939852245430365028195570327376952660380841611867271255283255613993669149247668224000000000000000000

n = 16
Probability of sampling a specific path from K_n:
183128253489049988956227768685942522498182431613727572162524233806103/58162052290933150592235732138903087476380671551040597371493148484696625723801600000000
 = 3.148586514331059e-18
Numerator (num orderings of edges leading to the path):
21062478660864250324062173458528475756409338492756619268477477451431572215501936895383635684968030696015138892447582728576466051633673409957851455315237202087116800000000000000000000
"""
