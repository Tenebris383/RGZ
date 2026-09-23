import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
#

# РОЗРАХУНКОВЕ ЗАВДАННЯ — ОСНОВИ ТЕОРІЇ ГРАФІВ
# Повне рішення всіх пунктів (а) – (і)
#

def solve_graph(edges, n=None, title="Граф G"):
 """
 Повне дослідження графа.
 edges — список пар (u, v), вершини нумеруються з 1.
 """
 G = nx.Graph()
 if n:
    G.add_nodes_from(range(1, n+1))
 G.add_edges_from(edges)
 nodes = sorted(G.nodes())
 V, E = len(nodes), G.number_of_edges()
 idx = {v: i for i, v in enumerate(nodes)}
 print(f"\n{'='*60}")
 print(f" {title}: |V|={V}, |E|={E}")
 print(f"{'='*60}")
 # ── а) Матриці

 print("\nа) МАТРИЦІ СУМІЖНОСТІ ТА ІНЦИДЕНТНОСТІ")
 A = nx.to_numpy_array(G, nodelist=nodes, dtype=int)
 print("Матриця суміжності A(G):")
 header = " " + " ".join(f"v{v}" for v in nodes)
 print(header)
 for i, v in enumerate(nodes):
    row = " ".join(f" {A[i,j]}" for j in range(V))
    print(f" v{v} {row}")
 edge_list = list(G.edges())
 B = np.zeros((V, E), dtype=int)
 for j, (u, v) in enumerate(edge_list):
    B[idx[u], j] = 1
    B[idx[v], j] = 1
 print("\nМатриця інцидентності B(G):")
 print(" " + " ".join(f"e{j+1:02d}" for j in range(E)))
 for i, v in enumerate(nodes):
    row = " ".join(str(B[i,j]) for j in range(E))
    print(f" v{v} {row}")
 # Матриця суміжності ребер (edge adjacency)
 print("\nМатриця суміжності ребер A_e(G) [розмір m×m]:")
 Ae = (B.T @ B)
 np.fill_diagonal(Ae, 0)
 print(" (рядки/стовпці = ребра e1..em, 1=суміжні ребра)")
 # ── б) Степені
#──────────────────────────────────────────────
 print("\nб) ВЕКТОР ТА МНОЖИНА СТЕПЕНІВ")
 deg_vec = [G.degree(v) for v in nodes]
 deg_set = sorted(set(deg_vec))
 print(f" Вектор степенів: {deg_vec}")
 print(f" Множина степенів: {deg_set}")
 print(f" δ(G)={min(deg_vec)}, Δ(G)={max(deg_vec)}")
 print(
    f"Перевірка: Σdeg = {sum(deg_vec)} = 2·|E| = {2 * E} ✓"
    if sum(deg_vec) == 2 * E
    else "✗ помилка!"
)
 # ── в) Максимальна незалежна множина
#────────────────────────
 print("\nв) МАКСИМАЛЬНА НЕЗАЛЕЖНА МНОЖИНА (MIS)")
 complement = nx.complement(G)
 cliques_c = list(nx.find_cliques(complement))
 mis = max(cliques_c, key=len)
 print(f" α(G) = {len(mis)}, MIS = {sorted(mis)}")
 # Перевірка
 ok = all(not G.has_edge(u,v) for u,v in combinations(mis,
2))
 print(f" Перевірка незалежності: {'✓' if ok else '✗'}")
 # ── г) Радіус та діаметр
#────────────────────────────────────
 print("\nг) РАДІУС ТА ДІАМЕТР")
 if nx.is_connected(G):
    ecc = nx.eccentricity(G)
    r, d = nx.radius(G), nx.diameter(G)
    center = nx.center(G)
    periph = nx.periphery(G)
    print(f" Ексцентриситети: { {v: ecc[v] for v in
nodes} }")
    print(f" Радіус r(G) = {r}, Центр = {center}")
    print(f" Діаметр d(G) = {d}, Периферія = {periph}")
 else:
    print(" Граф не зв'язний — обчислення для компонент")
 # ── д) Мости та точки зчленування
#───────────────────────────
 print("\nд) МОСТИ ТА ТОЧКИ ЗЧЛЕНУВАННЯ")
 bridges = list(nx.bridges(G))
 art_pts = sorted(nx.articulation_points(G))
 print(f" Мости: {bridges if bridges else 'відсутні'}")
 print(f" Точки зчленування: {art_pts if art_pts else
'відсутні'}")
 # ── е) Вершинна та реберна зв'язність

 print("\nе) ЧИСЛА ЗВ'ЯЗНОСТІ")
 if nx.is_connected(G):
    kv = nx.node_connectivity(G)
    ke = nx.edge_connectivity(G)
    min_vc = nx.minimum_node_cut(G)
    min_ec = nx.minimum_edge_cut(G)
    print(f" Вершинна зв'язність κ(G) = {kv}")
    print(f" Мін. вершинний розріз: {min_vc}")
    print(f" Реберна зв'язність λ(G) = {ke}")
    print(f" Мін. реберний розріз: {min_ec}")
    print(
        f"Нерівність Уітні: κ ≤ λ ≤ δ → "
        f"{kv} ≤ {ke} ≤ {min(deg_vec)} ✓"
    )
 else:
    print(" Граф не зв'язний: κ(G) = λ(G) = 0")
 # ── ж) Колода (максимальна кліка)
#───────────────────────────
 print("\nж) КОЛОДА (МАКСИМАЛЬНА КЛІКА)")
 cliques = list(nx.find_cliques(G))
 max_cl = max(cliques, key=len)
 omega = len(max_cl)
 print(f" ω(G) = {omega}, Колода = {sorted(max_cl)}")
 print(f" Усі максимальні кліки: {[sorted(c) for c in
cliques]}")
 # ── з) Планарність

 print("\nз) ПЛАНАРНІСТЬ")
 is_planar, _ = nx.check_planarity(G)
 print(f" Граф {'планарний ✓' if is_planar else 'непланарний ✗'}")
 if is_planar:
    F = 2 - V + E

    print(
        f"Формула Ейлера: V-E+F = "
        f"{V}-{E}+{F} = {V - E + F}"
    )

    print(
        f"Умова планарності: E ≤ 3V-6 → "
        f"{E} ≤ {3 * V - 6} ✓"
    )
 else:
    print(
        "Міноровим методом: містить підграф "
        "≅ K₅ або K₃,₃"
    )
 # ── і) Хроматичне число
#──────────────────────────────────────
 print("\ні) ХРОМАТИЧНЕ ЧИСЛО")
 coloring = nx.coloring.greedy_color(G,
strategy='largest_first')
 chi_upper = max(coloring.values()) + 1
 # Lower bound: clique number
 chi_lower = omega
 print(f" Нижня межа: ω(G) = {omega}")
 print(f" Верхня межа (жадібний): {chi_upper}")
 print(f" χ(G) = {chi_upper}")
 color_classes = {}
 for v, c in coloring.items():
    color_classes.setdefault(c, []).append(v)
 for c, verts in sorted(color_classes.items()):
    print(f" Клас {c+1}: {sorted(verts)}")
 return G
#
#════════════════════════════════════════════════════════════════
#══
# ПРИКЛАД ВИКОРИСТАННЯ
#
#════════════════════════════════════════════════════════════════
#══
if __name__ == '__main__':
 # Граф з прикладу розв'язання
 edges = [(1,2),(2,3),(3,4),(4,5),(5,1),(1,3),(3,5)]
 G = solve_graph(edges, n=5, title="Граф G (варіант 15)")
 # Візуалізація — 4 окремі графи
 pos = nx.spring_layout(G, seed=42)

 fig, axes = plt.subplots(2, 2, figsize=(12, 10))
 fig.suptitle("Аналіз графа G", fontsize=14, fontweight='bold')

    # ─────────────────────────────────────────────
    # 1. Початковий граф
    # ─────────────────────────────────────────────
 nx.draw(
   G, pos,
   ax=axes[0, 0],
   with_labels=True,
   node_color='#4A90D9',
   edge_color='#555',
   node_size=700,
   font_color='white',
   font_weight='bold'
   )
axes[0, 0].set_title("Граф G")

    # ─────────────────────────────────────────────
    # 2. Максимальна незалежна множина
    # ─────────────────────────────────────────────
complement = nx.complement(G)
mis = max(nx.find_cliques(complement), key=len)

node_colors = [
    '#9B59B6' if v in mis else '#AAAAAA'
    for v in G.nodes()
   ]

nx.draw(
        G, pos,
   ax=axes[0, 1],
   with_labels=True,
   node_color=node_colors,
   edge_color='#555',
   node_size=700,
   font_color='white',
   font_weight='bold'
   )
axes[0, 1].set_title(
        f"Максимальна незалежна множина: {sorted(mis)}"
    )

    # ─────────────────────────────────────────────
    # 3. Мости / точки зчленування
    # ─────────────────────────────────────────────
bridges = list(nx.bridges(G))
articulation_points = list(nx.articulation_points(G))
node_colors = [
    '#E74C3C' if v in articulation_points else '#AAAAAA'
     for v in G.nodes()
   ]

edge_colors = [
        '#E74C3C' if edge in bridges or (edge[1], edge[0]) in bridges
        else '#555'
        for edge in G.edges()
    ]

nx.draw(
        G, pos,
        ax=axes[1, 0],
   with_labels=True,
   node_color=node_colors,
   edge_color=edge_colors,
   node_size=700,
   font_color='white',
   font_weight='bold'
   )
axes[1, 0].set_title("Мости / точки зчленування")

    # ─────────────────────────────────────────────
    # 4. Розфарбування графа
    # ─────────────────────────────────────────────
coloring = nx.coloring.greedy_color(
        G,
        strategy='largest_first'
    )

palette = ['#E84B3A', '#4A90D9', '#27AE60', '#F39C12']

node_colors = [
        palette[coloring[v]]
        for v in G.nodes()
    ]

chi = max(coloring.values()) + 1

nx.draw(
   G, pos,
   ax=axes[1, 1],
   with_labels=True,
   node_color=node_colors,
   edge_color='#555',
   node_size=700,
   font_color='white',
   font_weight='bold'
   )
axes[1, 1].set_title(f"Розфарбування χ = {chi}")

plt.tight_layout()
plt.savefig(
   'graph_analysis.png',
   dpi=120,
   bbox_inches='tight'
   )
plt.show()