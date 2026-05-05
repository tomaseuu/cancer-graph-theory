import networkx as nx
import random
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np


random.seed(42)

# 3. The Graph

# initialize an undirected graph
G = nx.Graph()

# load interaction data 
 # read this file as UTF-8, but if there is a BOM at the start, ignore it (encoding=utf-8-sig)
with open("data/interacting_proteins.txt", "r", encoding="utf-8-sig") as file:   
    for line in file:
        protein1, protein2 = line.strip().split()  # strip removes extra spaces and split splits line by spaces
        G.add_edge(protein1.strip().upper(), protein2.strip().upper())  # creates an edge between two nodes - upper makes it all caps


# 3.1 The Shortest Path
# testing biological hypothesis: Tumor-causing (oncogenic) genes are more closely
# connected to each other in the protein interaction network than random genes are

# load onco_genes data
 # read this file as UTF-8, but if there is a BOM at the start, ignore it (encoding = utf-8-sig)
oncogenes = []
with open("data/onco_genes.txt", "r", encoding="utf-8-sig") as file:
    for line in file:
        if line.strip():
            oncogenes.append(line.strip().upper())

# give me all unique pairs of two genes from the oncogenes list
# nested loop
onco_pairs = []
for i in range(len(oncogenes)): 
    for j in range(i + 1, len(oncogenes)): # don't want to repeat or reverse any pair, so start at i + 1
        onco_pairs.append((oncogenes[i], oncogenes[j]))

distances = [] # storing all shortest path distance

for gene1, gene2 in onco_pairs:
    if nx.has_path(G, gene1, gene2):  # checks if there is a path
        length = nx.shortest_path_length(G, source=gene1, target=gene2) #CALCULATE SHORTEST PATH LENGTH between two nodes
        distances.append(length)


# checks to see if it collected any distances
if distances:
    avg_onco_distance = sum(distances) / len(distances) # calculate average
    print(f"Average shortest path between tumorigenic genes: {avg_onco_distance:.3f}") # keeps 3 digits after the decimal point
else:
    print("There was no connected oncogene pairs found.")

# set -> converts it into a set so you can do set operations
# create a set of nodes from G
all_nodes = set(G.nodes())
non_oncogenes = sorted(all_nodes - set(oncogenes)) 

# RANDOM SETS
samples = 1000 # random 1000 sets
sample_size = len(oncogenes)
background_scores = [] # use to store 1000 average shortest path scores

for i in range(samples):
    sample_genes = random.sample(non_oncogenes, sample_size)
    pairs =[]
    for i in range(len(sample_genes)):
        for j in range(i + 1, len(sample_genes)):
            pairs.append((sample_genes[i], sample_genes[j]))

    distances = []
    for gene1, gene2 in pairs:
        if nx.has_path(G, gene1, gene2):
            length = nx.shortest_path_length(G, source=gene1, target=gene2)
            distances.append(length)
    
    if distances: 
        avg_dist = sum(distances) / len(distances)
        background_scores.append(avg_dist)

# plotting the distribution of average shortest path

plt.figure(figsize= (10, 6))
plt.hist(background_scores, bins=30, color= "lightblue", edgecolor="black")
plt.axvline(avg_onco_distance, color="red", linestyle="dashed",linewidth=2, label= "Oncogene Score")
plt.title("Distribution of Average Shortest Path (Random Sets)")
plt.xlabel("Average Shortest Path")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True)
plt.show()


# P-VALUE
# asking: how many random values are as extreme or more extreme than my observed value?
num_extreme = sum(score <= avg_onco_distance for score in background_scores)
p_value = num_extreme / len(background_scores)
print(f"p-value: {p_value:.3f}")





# 4. Graph Diffusion and Random Walks


# 4.1.1
#add edges based on the image.
G = nx.Graph()

edges = [
    ('A', 'B'), ('A', 'C'), ('A', 'D'), ('A', 'E'),
    ('B', 'F'), ('C', 'G'), ('D', 'H'), ('E', 'I'),
    ('F', 'J'), ('G', 'K'), ('H', 'L'), ('I', 'M'),
    ('J', 'K'), ('K', 'L'), ('L', 'M')
]

G.add_edges_from(edges)

for u, v in G.edges():
    G[u][v]['weight'] = 1  # default weight

#G['A']['C']['weight'] = 10     # Test Case # 1 (very heavy edge)
#G['A']['D']['weight'] = 0.1    # Test Case # 2 (reduce weight to near zero)
#G['E']['I']['weight'] = 15     # Test Case # 3 (very heavy edge to a leaf node)
#G['K']['L']['weight'] = 8       # Test Case # 4 ( midpoint path heavier )
 
restart_prob= 0.3
start_node = 'A'
num_steps = 100000  # simulate 100,000 steps

visit_counts = {node: 0 for node in G.nodes()} # give every node a visit count of 0
current = start_node

# random walk with restarts:  
#    restart at the start node or move to a neighboring node based on edge weights
for i in range(num_steps):
    visit_counts[current] += 1
    if random.random() < restart_prob:
        current = start_node  # restart
    else:
        neighbors = list(G.neighbors(current))
        if neighbors:
            weights = [G[current][n].get('weight', 1) for n in neighbors]
            current = random.choices(neighbors, weights=weights)[0]

# stationary distribution
# for each node - divide its visit count by total visits to get frequency
total_visits = sum(visit_counts.values())
stationary = {node: count / total_visits for node, count in visit_counts.items()}

# get layout positions for nodes from nx library
pos = nx.shell_layout(G)
node_colors = [stationary[node] for node in G.nodes()]
labels_freq = {n: f"{stationary[n]:.2f}" for n in G.nodes()}
labels_node = {n: n for n in G.nodes()}


# normalizing color to make sure it looks consistent and nice
norm = mcolors.Normalize(vmin=min(node_colors), vmax=max(node_colors))
cmap = plt.cm.YlOrRd
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 

# edge widths by weight
edge_weights = [G[u][v].get('weight', 1) for u, v in G.edges()]
edge_widths = [w for w in edge_weights]  # scale if needed

# plot
fig, ax = plt.subplots(figsize=(10, 7))
nx.draw(
    G,
    pos,
    node_color=node_colors,
    cmap=cmap,
    node_size=1600,
    width=edge_widths,  # widths changed based off edge (heavy vs light)
    edge_color='blue',
    ax=ax
)
# draw frequency labels nodes
nx.draw_networkx_labels(G, pos, labels=labels_freq, font_size=12, font_color='black', ax=ax)

# draw node names (A, B, C...) nodes so it can show
offset_pos = {k: (v[0], v[1] + 0.07) for k, v in pos.items()}
nx.draw_networkx_labels(G, offset_pos, labels=labels_node, font_size=10, font_color='black', ax=ax)

# add colorbar
cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.03)
cbar.set_label("Stationary Frequency")

ax.set_title("Random Walk with Restart from Node A", pad=20)
ax.axis('off')
plt.show()



# 4.1.2 (CHANGE y = 1 or y = 0 to test) 
degrees = dict(G.degree()) 
stationary_scores = stationary # results from the RWR walk (stationary frequencies)

nodes_sorted = list(G.nodes())
degree_list = [degrees[node] for node in nodes_sorted]
stationary_list = [stationary_scores[node] for node in nodes_sorted]

correlation, _ = spearmanr(degree_list, stationary_list) # from scipy.stats (easier)
print(f"Spearman rank correlation (degree vs. stationary): {correlation:.4f}")


# 4.2 REAL DATA

#create new graph for this section
G = nx.Graph()

with open("data/interacting_proteins.txt", "r", encoding="utf-8-sig") as file:   
    for line in file:
        protein1, protein2 = line.strip().split()
        G.add_edge(protein1.strip().upper(), protein2.strip().upper())

# run_rwr function - simulates random walk with restart from a given start node
def run_rwr(graph, start_node, restart_prob=0.3, num_steps=100000):
    visit_counts = {node: 0 for node in graph.nodes()}
    current = start_node

    # random walk with restarts:  
    #    restart at the start node or move to a neighboring node based on edge weights
    for i in range(num_steps):
        visit_counts[current] += 1
        if random.random() < restart_prob:
            current = start_node  # restart
        else:
            neighbors = list(graph.neighbors(current))
            if neighbors:
                weights = [graph[current][n].get('weight', 1) for n in neighbors]
                current = random.choices(neighbors, weights=weights)[0]

    # stationary distribution
    # for each node - divide its visit count by total visits to get frequency
    total_visits = sum(visit_counts.values())
    stationary = {node: count / total_visits for node, count in visit_counts.items()}
    return stationary



# keep only genes that exist in the graph
filtered_oncogenes = []
for gene in oncogenes:
    if gene in G.nodes():
        filtered_oncogenes.append(gene)
oncogenes = filtered_oncogenes

# calculating average proximity for oncogenes using RWR
restart_prob = 0.3
num_steps = 100000
proximity_scores = []

for gene in oncogenes:
    stationary = run_rwr(G, gene, restart_prob, num_steps)
    other_oncos = []    # make a list of all other oncogenes except current one
    for g in oncogenes:
        if g != gene:
            other_oncos.append(g)
    total = 0    # sum of stationary values of all other oncogenes
    for g in other_oncos:
        total += stationary.get(g,0)
    avg_freq = total / len (other_oncos) # average visit frequency
    proximity_scores.append(avg_freq)

oncogene_rwr_score = sum(proximity_scores) / len(proximity_scores)
print(f"RWR proximity score for oncogenes: {oncogene_rwr_score:.6f}")



# background distribution using random non-oncogenes
samples = 1000 
sample_size = len(oncogenes)
background_rwr_scores = []
non_oncogenes = sorted(set(G.nodes()) - set(oncogenes))

for i in range(samples):
    sample_genes = random.sample(non_oncogenes, sample_size)
    proximity_scores = []

    for gene in sample_genes:
        stationary = run_rwr(G, gene, restart_prob, num_steps)
        others = []
        for g in sample_genes:
            if g != gene:
                others.append(g)

        total = 0
        for g in others:
            total += stationary.get(g, 0)
        avg_freq = total / len(others)
        proximity_scores.append(avg_freq)

    background_rwr_scores.append(sum(proximity_scores) / len(proximity_scores))

x_max = max(np.percentile(background_rwr_scores, 95), oncogene_rwr_score * 1.05)
plt.figure(figsize=(10, 6))
plt.hist(background_rwr_scores, bins=30, color="lightblue", edgecolor="black")
plt.axvline(oncogene_rwr_score, color="red", linestyle="dashed", linewidth=2, label="Oncogene RWR Score")
plt.xlim(0, x_max)
plt.title("RWR Proximity Scores (Background vs Oncogenes)")
plt.xlabel("Average RWR Proximity Score")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True)
plt.show()


# p-value
num_extreme = sum(score >= oncogene_rwr_score for score in background_rwr_scores)
p_value = num_extreme / len(background_rwr_scores)
print(f"RWR p-value: {p_value:.4f}")


# 4.2.3 predict novel tumorigenic genes

# first initialize frequency counter for all the nodes
combine_freqs = {node: 0 for node in G.nodes()}

# count how many times each gene would get visit in RWRs 
# starting from each known oncogenes
for gene in oncogenes:
    stationary = run_rwr(G, gene, restart_prob=0.3, num_steps=100000)
    for node, freq in stationary.items():
        combine_freqs[node] += freq
# average the frequency across all walks
average_freqs = {}
for node, freq in combine_freqs.items():
    average_freqs[node] = freq / len(oncogenes)

# remove og oncogenes from the list
og_genes = {}
for gene, score in average_freqs.items():
    if gene not in oncogenes:
        og_genes[gene] = score

# sort by highest average frequency
ranked_genes = sorted(og_genes.items(), key=lambda item: item[1], reverse=True)

# extract top 30 predicted genes
top_30_predicted = []
for gene, i in ranked_genes[:30]:
    top_30_predicted.append(gene)

# Save to file
with open("Thomas_Le_onco_predictions.txt", "w") as f:
    for gene in top_30_predicted:
        f.write(f"{gene}\n")

print("Top 30 novel tumorigenic genes written to Thomas_Le_onco_predictions.txt")



# 4.3 Graph Randomization   

# copy graph and randomize it
G_randomized = G.copy()
nx.double_edge_swap(G_randomized, nswap=10 * G_randomized.number_of_edges(), max_tries=100 * G_randomized.number_of_edges())

# using the RWR function from before

# filter oncogenes in randomized graph


random_oncogenes = []
for gene in oncogenes:
    if gene in G_randomized.nodes:
        random_oncogenes.append(gene)

restart_prob = 0.3
num_steps = 100000

random_proximity_scores = []

for gene in random_oncogenes:
    stationary = run_rwr(G_randomized, gene, restart_prob, num_steps)
    others = []
    for g in random_oncogenes:
        if g != gene:
            others.append(g)
    total = 0 
    for g in others:
        if g in stationary:
            total += stationary[g]
        else:
            total += 0
    average_freq = total / len(others)
    random_proximity_scores.append(average_freq)

random_onco_rwr_score = sum(random_proximity_scores) / len(random_proximity_scores)
print(f"Randomized graph RWR proximity score for oncogenes: {random_onco_rwr_score: .6f}")


# backgruond distribution in randomized graph
samples = 1000
sample_size = len(random_oncogenes)
random_background_rwr_scores = []
random_non_oncogenes = sorted(set(G_randomized.nodes()) - set (random_oncogenes))

for i in range(samples):
    sample_genes = random.sample(random_non_oncogenes, sample_size)
    proximity_scores = []

    for gene in sample_genes:
        stationary = run_rwr(G_randomized, gene, restart_prob, num_steps)
        others = []
        for g in random_oncogenes:
            if g != gene:
                others.append(g)
        total = 0 
        for g in others:
            if g in stationary:
                total += stationary[g]
            else:
                total += 0
        average_freq = total / len(others)
        proximity_scores.append(average_freq)
    random_background_rwr_scores.append(sum(proximity_scores) / len(proximity_scores))


#plot

x_max = max(np.percentile(random_background_rwr_scores, 95), random_onco_rwr_score * 1.05)
plt.figure(figsize=(10, 6))
plt.hist(random_background_rwr_scores, bins=30, color="lightblue", edgecolor="black")
plt.axvline(random_onco_rwr_score, color="red", linestyle="dashed", linewidth=2, label="Oncogene RWR Score (Randomized)")
plt.xlim(0, x_max)
plt.title("RWR Proximity Scores (Randomized Graph)")
plt.xlabel("Average RWR Proximity Score")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True)
plt.show()

#p-value
num_extreme = sum(score >= random_onco_rwr_score for score in random_background_rwr_scores)
rand_p_value = num_extreme / len(random_background_rwr_scores)
print(f"RWR p-value (randomized graph): {rand_p_value:.4f}")