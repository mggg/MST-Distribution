"""
A small file containing all of the functions necessary to compute the probability
of a particular minimum spanning tree.

Authors:
    Pietro Poggi-Corradini (pietro@math.ksu.edu)
    Peter Rock (peter@mggg.org)
"""

from itertools import permutations
from fractions import Fraction
import networkx as nx
from multiprocessing import Pool, Manager, Queue
from tqdm import tqdm
from math import factorial


def broken_cycle(graph, spanning_forest, edge_u):
    """
    Given a spanning forest and an edge not in that forest, return the broken cycle
    corresponding to that edge.

    Args:
        spanning_forest (nx.Graph): A spanning forest (trees are forests of 1).
        edge_u (tuple): An edge not in the spanning forest.

    Returns:
        set: The broken cycle corresponding to the edge.
    """
    try:
        gamma = nx.shortest_path(spanning_forest, edge_u[0], edge_u[1])
        broken_cycle_u = set()
        broken_cycle_u.add(edge_u)
        for i in range(len(gamma) - 1):
            if (gamma[i], gamma[i + 1]) in list(graph.edges()):
                broken_cycle_u.add((gamma[i], gamma[i + 1]))
            else:
                broken_cycle_u.add((gamma[i + 1], gamma[i]))
        return broken_cycle_u
    except nx.NetworkXNoPath:
        return set()


def product_cycles(graph, spanning_tree, perm_U):
    """
    Given a graph, a spanning tree, and a permutation of the edges not in the spanning tree,
    compute the product of the sizes of the broken cycles corresponding to the edges in the
    permutation.

    Args:
        graph (nx.Graph): The graph.
        spanning_tree (nx.Graph): The spanning tree.
        perm_U (list): The permutation of the edges not in the spanning tree.

    Returns:
        int: The product of the sizes of the broken cycles.
    """
    broken_cycle_edges = set()
    prod_broken_cycle_union_len = 1
    seen = set()  # To avoid recomputation of broken cycles

    for j in range(len(perm_U)):
        new_elements = set()
        for i in range(j + 1):
            if perm_U[i] not in seen:
                seen.add(perm_U[i])
                new_elements.update(broken_cycle(graph, spanning_tree, perm_U[i]))
        if new_elements:
            broken_cycle_edges.update(new_elements)
        prod_broken_cycle_union_len *= len(broken_cycle_edges)

    return prod_broken_cycle_union_len


def process_external_permutation(graph, spanning_tree, perm_U):
    r"""
    For a particular $\sigma \in S_{m-n+1}$, computes the product

    $$
    \prod_{j=1}^{m-n+1} \frac{1}{d_j(\sigma)}
    $$

    where $d_j(\sigma)$ is the size of the union of the broken cycles corresponding
    to the first $j$ edges in $U$ given the order determined by $\sigma$.

    Args:
        graph (nx.Graph): The graph.
        spanning_tree (nx.Graph): The spanning tree.
        perm_U (list): The permutation of the edges not in the spanning tree.

    Returns:
        Fraction: The value of the product as an exact fraction.
    """
    prod = product_cycles(graph, spanning_tree, perm_U)
    return Fraction(1, prod) if prod != 0 else 0


def update_result(result, x, queue):
    """
    Used as a callback function to update the result value and put a value in the queue
    during the multithreading process.

    Args:
        result (Value): The shared value to update.
        x (Fraction): The value to add to the result.
        queue (Queue): The queue to put a value in.
    """
    result.value += x
    queue.put(1)


def mst_prob_external(G, T, show_progress=False):
    r"""
    Computes the probability of a minimum spanning tree given a graph and a spanning tree
    using the external formula. The external formula is given as follows:

    Let $U = E(G) \setminus E(T)$. Given an edge $u = \{x, y\} \in U$, there is a unique path
    $\gamma_u$ in $T$ connecting $x$ to $y$. Consider the cycle

    $$
    C_u := \gamma_u \cup \{u\}.
    $$

    We say that $C_u$ is the _broken cycle_ corresponding to $u$.

    Enumerate $U = \{u_i\}_{i=1}^{m-n+1}$.
    Given a permutation $\sigma \in S_{m-n+1}$, for $j = 1, \ldots, m-n+1$, write

    $$
    D_j(\sigma) := \bigcup_{i=1}^j C_{u_{\sigma(i)}}.
    $$

    In words, $D_j(\sigma)$ contains all the edges in the broken cycles created by the
    first $j$ edges in $U$, given the order determined by $\sigma$. Set
    $d_j(\sigma) := |D_j(\sigma)|$. Then, we get the following _external formula_:

    $$
    \mathbb{P}_{\rm MST}(T) = \sum_{\sigma \in S_{m-n+1}} \prod_{j=1}^{m-n+1} \frac{1}{d_j(\sigma)}.
    $$

    Args:
        G (nx.Graph): The graph.
        T (nx.Graph): A spanning tree of the graph $G$.
        show_progress (bool): Whether to show progress bars. Defaults to False.

    Returns:
        Fraction: The probability of the minimum spanning tree as an exact fraction.
    """

    U = [e for e in G.edges() if e not in T.edges()]
    result = Manager().Value(Fraction, 0)
    queue = Manager().Queue()
    total_permutations = factorial(len(U))

    with Pool() as pool:
        for perm in tqdm(
            permutations(U),
            total=total_permutations,
            desc="Generating permutations...",
            disable=not show_progress,
        ):
            pool.apply_async(
                process_external_permutation,
                args=(G, T, perm),
                callback=lambda x: update_result(result, x, queue),
            )

        pool.close()

        with tqdm(
            total=total_permutations,
            desc="Processing permutations...",
            disable=not show_progress,
        ) as pbar:
            for _ in range(total_permutations):
                queue.get()  # Wait for a task to complete
                pbar.update()

        pool.join()

    return result.value


def cutset_size(graph, spanning_forest):
    """
    Given a graph and a spanning forest, compute the number of cutsets in the graph.
    That is to say, that an edge is a member of the cut set if there is no path from
    one endpoint to the other contained completely in the spanning forest.

    Args:
        graph (nx.Graph): The graph.
        spanning_forest (nx.Graph): The spanning forest.

    Returns:
        int: The number of cutsets in the graph with respect to the spanning forest.
    """

    H = nx.Graph()
    H.add_edges_from(spanning_forest)
    H.add_nodes_from(graph.nodes())

    c = 0
    for e in graph.edges():
        if e not in spanning_forest:
            if len(broken_cycle(graph, H, e)) == 0:
                c += 1
    return c


def grow_forests(order):
    """
    Given an ordering on the edges for a particular spanning tree,
    grow a set of forests by adding edges in the order given by the ordering.
    The returned forests from this function contain only the non-trivial
    trees within the forest, so trees containing a single vertex are not
    included here even though they are a part of the spanning forest.

    Args:
        order (list): The ordering of the edges in a given spanning tree.

    Returns:
        list: A list of forests contained within the given spanning tree.
    """

    forests = [[e for e in list(order)[: (j + 1)]] for j in range(len(order) - 1)]
    return forests


def product_cutsets(graph, perm_T, cutset_cache):
    """
    Given a graph and a permutation of the edges in a spanning tree, compute the product
    of the number of cutsets in the graph with respect to the forests generated by the
    permutation and the function `grow_forests`.

    Args:
        graph (nx.Graph): The graph.
        perm_T (list): The permutation of the edges in the given spanning tree.
        cutset_cache (dict): A cache for storing the number of cutsets for a particular forest.

    Returns:
        int: The product of the number of cutsets in the graph with respect to the forests.
    """
    forests = grow_forests(perm_T)
    c = len(graph.edges())
    for forest in forests:
        forest_tuple = tuple(sorted(forest))
        if forest_tuple not in cutset_cache:
            cutset_cache[forest_tuple] = cutset_size(graph, forest)
        c *= cutset_cache[forest_tuple]
    return c


def process_internal_permutation(graph, perm_T, cutset_cache):
    r"""
    For a particular $\sigma \in S_{n-1}$, computes the product

    $$
    \prod_{j=1}^{n-1} \frac{1}{c_j(\sigma)}
    $$

    where $c_j(\sigma)$ is the the number of edges that do not form a cycle within
    the first $j$ edges in the given permutation `perm_T` = $\sigma(E(T))$ of
    the edges of a particular spanning tree $T$.

    Args:
        graph (nx.Graph): The graph.
        perm_T (list): The permutation of the edges in the given spanning tree.
        cutset_cache (dict): A cache for storing the number of cutsets for a particular forest.

    Returns:
        Fraction: The reciprocal of the product of the cutsets as an exact fraction.
    """

    prod = product_cutsets(graph, perm_T, cutset_cache)
    return Fraction(1, prod) if prod != 0 else 0


def mst_prob_internal(G, T, show_progress=False):
    r"""
    Computes the probability of a minimum spanning tree given a graph and a spanning tree
    using the internal formula. The internal formula is given as follows:

    Let $T$ be a spanning tree of $G$ and let $\omega \in S_{n-1}$. Let
    $E(T) = \{e_1, \ldots, e_{n-1}\}$ and
    $\{\tilde{e}_1, \ldots, \tilde{e}_{n-1}\} = \{e_{\omega(1)}, \ldots, e_{\omega(n-1)}\}$.
    Then we may construct a sequence of forests
    $F_0, \ldots, F_{n-1}$ by letting $V(F_0) = V(G)$ with $E(F_0) = \emptyset$, and

    $$
    V(F_j) = V(G)\qquad \text{with} \qquad E(F_j) = \{\tilde{e}_1, \ldots, \tilde{e}_{j}\}.
    $$

    where $1 \leq j \leq n-2$. We may then define the cut set $C_j(\omega)$ to be the set of
    edges between distinct connected components of $F_j$, and let $c_j(\omega) = |C_j(\omega)|$.

    Then, the _internal formula_ for the probability of the minimum spanning tree is given by

    $$
    \mathbb{P}_{\rm MST}(T)=\sum_{\omega\in S_{n-1}} \prod_{k=0}^{n-2} \frac{1}{c_k(\omega)}.
    $$


    Args:
        G (nx.Graph): The graph.
        T (nx.Graph): A spanning tree of the graph $G$.
        show_progress (bool): Whether to show progress bars. Defaults to False.

    Returns:
        Fraction: The probability of the minimum spanning tree as an exact fraction.
    """
    # form the set T_edges
    total_permutations = factorial(len(T.edges()))
    result = Manager().Value(Fraction, 0)
    queue = Queue()
    cutset_cache = Manager().dict()

    with Pool() as pool:
        for perm in tqdm(
            permutations(T.edges),
            total=total_permutations,
            desc="Generating permutations...",
            disable=not show_progress,
        ):
            pool.apply_async(
                process_internal_permutation,
                args=(G, perm, cutset_cache),
                callback=lambda x: update_result(result, x, queue),
            )

        pool.close()

        with tqdm(
            total=total_permutations,
            desc="Processing permutations...",
            disable=not show_progress,
        ) as pbar:
            for _ in range(total_permutations):
                queue.get()  # Wait for a task to complete
                pbar.update()

        pool.join()

    return result.value


def compute_mst_prob(graph, tree, show_progress=False, force=False, print_method=False):
    """
    Compute the probability of a minimum spanning tree given a graph and a spanning tree.
    The function will automatically choose the most efficient method (i.e. internal or
    external) to compute the probability based on the number of permutations required.
    If the number of permutations is below 1e8, the function will compute the probability
    using the appropriate formula. Otherwise, the function will raise an error unless the
    `force` parameter is set to `True`.

    Args:
        graph (nx.Graph): The graph.
        tree (nx.Graph): A spanning tree of the graph $G$.
        show_progress (bool): Whether to show progress bars. Defaults to False.
        force (bool): Whether to force the computation of the probability. Defaults to False.
        print_method (bool): Whether to print the method used to compute the probability.
            Defaults to False.

    Returns:
        Fraction: The probability of the minimum spanning tree as an exact fraction.

    Raises:
        ValueError: If the number of permutations is too large to compute.
    """
    external_permutations = factorial(len(graph.edges()) - len(tree.edges()))
    internal_permutations = factorial(len(tree.edges()))

    if min(external_permutations, internal_permutations) > 1e8 and not force:
        raise ValueError(
            "Too many permutations to compute. To continue with the computation anyway, "
            "re-run the 'compute_mst_prob' function with force=True."
        )

    if external_permutations < internal_permutations:
        if print_method:
            print("Using external formula.")
        return mst_prob_external(graph, tree, show_progress=show_progress)
    else:
        if print_method:
            print("Using internal formula.")
        return mst_prob_internal(graph, tree, show_progress=show_progress)
