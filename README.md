# Cancer and Graph Theory

This project explores how rare gene mutations in cancer can still matter, even if they don’t show up often across patients. Instead of only looking at how often a gene is mutated, I used **graph theory** to see how these genes are connected in a **protein–protein interaction network**. The main idea is that even “rare” mutations might still be close to important cancer genes, meaning they could still play a big role in tumor development.

This project was inspired by a 2015 study published in *Nature Genetics* that used a graph diffusion method to find hidden cancer driver genes. I wanted to recreate a smaller version of that idea — to see if I could find meaningful connections between known oncogenes using Python and graph algorithms.

---

## Dataset
- Known kidney cancer oncogenes: `data/onco_genes.txt`  
- Protein–protein interactions: `data/interacting_proteins.txt`

Each protein is treated as a **node**, and every interaction between two proteins is an **edge** in the graph.

---

## Objectives
- Build a **graph** from the protein interaction data using NetworkX.  
- Measure **average shortest paths** between oncogenes and compare them to random gene sets.  
- Run **Random Walk with Restart (RWR)** to see how mutation effects spread across the network.  
- Test how **edge weights** and **restart probabilities (γ)** change diffusion behavior.  
- Predict **new tumor-related genes** based on their proximity to known oncogenes.  
- Check significance using **randomized graphs** and **p-values**.

---

## Tools & Libraries
- **Python 3**
- **NetworkX** – graph modeling and random walks  
- **NumPy** – matrix operations and random sampling  
- **SciPy** – stats and correlation tests  
- **Matplotlib** – histograms and graph visualizations  
- **Random** – reproducible simulations

---

## Key Findings
- The average shortest path between oncogenes was **5.667**, not significantly closer than random sets (p = 0.401).  
- The Random Walk with Restart (RWR) analysis showed **oncogenes are much more connected** than random genes (p = 0.0030).  
- Randomizing the network removed that pattern (p = 0.3890), showing the **real network structure matters**.  
- Predicted **30 possible new tumor-related genes** based on proximity scores.  
- Edge weights and restart probability clearly affected how the walk spread — higher weights pulled the walk toward specific nodes, while lower weights pushed it away.

---

## Project Report
For full explanations, figures, and math steps, see:  
📄 [Cancer and Graph Theory.pdf](https://github.com/user-attachments/files/22854158/Cancer.and.Graph.Theory.pdf)



---

## Author
**Thomas Le**  
Student in Computer Science & Bioinformatics
