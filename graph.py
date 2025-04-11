import networkx as nx
import matplotlib.pyplot as plt  

# To create an empty undirected graph
G = nx.DiGraph(directed=True)


node_map = {
    "Glucose": 1,
    "G6P": 2,
    "F6P": 3,
    "F16P": 4,
    "2P": 5,
    "3P": 6,
    "13B": 7,
    "G3P": 8,
    "DP": 9,
    "PP": 10,
    "Pyruvate": 11,
    "aCoa": 12,


    "mCoa": 13,
    "aACP": 14,
    "mACP": 15,
    "Beta-K-ACP": 16,
    "Beta-H-ACP": 17,
    "e-ACP": 18,
    "b-ACP": 19,
    "h-ACP": 20,
    "o-ACP": 21,
    "dec-ACP": 22,
    "dod-ACP": 23,
    "tetra-ACP": 24,
    "pan-ACP": 25,
    "palm-ACP": 26,
    "ACP": 27,
}

reverse_node_map = {value: key for key, value in node_map.items()}

def add_edge(key1, key2):
    G.add_edge(reverse_node_map[key1], reverse_node_map[key2])


# ------------------------------------------------------------------


for name in node_map.keys():
    G.add_node(name)
  

add_edge(1, 2)
add_edge(2, 3)
add_edge(3, 4)
add_edge(4, 9)
add_edge(4, 8)
add_edge(8, 7)
add_edge(7, 6)
add_edge(6, 5)
add_edge(5, 10)
add_edge(10, 11)


add_edge(11, 12)


add_edge(12, 14)
add_edge(12, 13)
add_edge(13, 15)
add_edge(15, 16)
add_edge(16, 17)
add_edge(17, 18)
add_edge(18, 19)
add_edge(19, 16)

add_edge(19, 20)
add_edge(20, 21)
add_edge(21, 22)
add_edge(22, 23)
add_edge(23, 24)
add_edge(24, 25)
add_edge(25, 26)
add_edge(25, 27)


nx.draw(G, with_labels=True)
plt.savefig('plotgraph.png', dpi=600)
plt.show()
