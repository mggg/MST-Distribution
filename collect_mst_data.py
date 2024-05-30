from mst_lib import *
import networkx as nx
import numpy as np


def lattice_hexagon_points(radius):
    points = set([(0, 0)])

    current_rad = radius

    while current_rad > 0:
        # Diagonal sides
        a, b = current_rad, current_rad
        while b >= 0:
            points.add((a, b))
            points.add((a, -b))
            points.add((-a, b))
            points.add((-a, -b))
            a += 1
            b -= 1

        # Top and bottom
        a, b = current_rad, current_rad
        while a >= 0:
            points.add((a, b))
            points.add((a, -b))
            points.add((-a, b))
            points.add((-a, -b))
            a -= 2
        current_rad -= 1

    return points


def make_hex_grid(radius):
    G = nx.Graph()

    vertices = lattice_hexagon_points(radius)
    for v in vertices:
        G.add_node(v)

    for v in vertices:
        a, b = v
        possible_neighbors = [
            (a + 1, b + 1),
            (a + 1, b - 1),
            (a - 1, b + 1),
            (a - 1, b - 1),
            (a + 2, b),
            (a - 2, b),
        ]
        for n in possible_neighbors:
            if n in vertices:
                G.add_edge(v, n)

    return G


if __name__ == "__main__":

    import json
    import jsonlines

    print("Computing probabilities for hexagonal grid of radius 1")
    radius = 1
    G = make_hex_grid(radius)

    n_spanning_trees = int(np.linalg.det(nx.laplacian_matrix(G).toarray()[:-1, :-1]))

    with jsonlines.open("./enumerations/hex_grid_radius_1.jsonl", "w") as writer:
        for T in tqdm(nx.SpanningTreeIterator(G), total=n_spanning_trees):
            diameter = nx.diameter(T)
            probability = compute_mst_prob(G, T)
            writer.write(
                {
                    "diameter": diameter,
                    "probability": str(probability),
                    "tree_edges": list(T.edges),
                }
            )

    for i in range(2, 4):
        print(f"Computing probabilities for {i}x{i} grid")
        with jsonlines.open(f"./enumerations/mst_{i}x{i}_grid.jsonl", "w") as writer:
            G = nx.grid_2d_graph(i, i)
            n_spanning_trees = int(
                np.linalg.det(nx.laplacian_matrix(G).toarray()[:-1, :-1])
            )

            for T in tqdm(nx.SpanningTreeIterator(G), total=n_spanning_trees):
                diameter = nx.diameter(T)
                probability = compute_mst_prob(G, T)

                writer.write(
                    {
                        "diameter": diameter,
                        "probability": str(probability),
                        "tree_edges": list(T.edges),
                    }
                )

    for i in range(2, 5):
        if i == 3:
            continue
        with jsonlines.open(f"./enumerations/mst_{3}x{i}_grid.jsonl", "w") as writer:
            G = nx.grid_2d_graph(3, i)
            n_spanning_trees = int(
                np.linalg.det(nx.laplacian_matrix(G).toarray()[:-1, :-1])
            )

            for T in tqdm(nx.SpanningTreeIterator(G), total=n_spanning_trees):
                diameter = nx.diameter(T)
                probability = compute_mst_prob(G, T)

                writer.write(
                    {
                        "diameter": diameter,
                        "probability": str(probability),
                        "tree_edges": list(T.edges),
                    }
                )

    for i in range(2, 8):
        print(f"Computing probabilities for complete graph with {i} nodes")
        with jsonlines.open(f"./enumerations/mst_complete_{i}.jsonl", "w") as writer:
            G = nx.complete_graph(i)
            n_spanning_trees = int(
                np.linalg.det(nx.laplacian_matrix(G).toarray()[:-1, :-1])
            )

            for T in tqdm(nx.SpanningTreeIterator(G), total=n_spanning_trees):
                diameter = nx.diameter(T)
                probability = compute_mst_prob(G, T)

                writer.write(
                    {
                        "diameter": diameter,
                        "probability": str(probability),
                        "tree_edges": list(T.edges),
                    }
                )
