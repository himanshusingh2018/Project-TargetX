"""
Module 3: Disease Modeling & Target Discovery
Script: disease_knowledge_graph.py
Constructs a biological knowledge graph and evaluates network centrality
to identify and rank high-priority therapeutic targets in disease pathways.
"""

import os
import networkx as nx
import pandas as pd


def build_oncology_kg() -> nx.Graph:
    """Builds a representative protein-protein and disease-pathway interaction graph.
    Nodes: Genes/Proteins and Biological Pathways.
    Edges: Physical interaction, phosphorylation, or pathway membership.
    """
    G = nx.Graph()

    # Define biological interactions in the EGFR/MAPK/PI3K oncogenic network
    interactions = [
        # Upstream Receptors
        ("EGFR", "GRB2", "binds"),
        ("EGFR", "PIK3CA", "phosphorylates"),
        ("ERBB2", "EGFR", "heterodimerizes"),
        ("ERBB2", "GRB2", "binds"),
        
        # MAPK/ERK Signaling Cascade
        ("GRB2", "SOS1", "recruits"),
        ("SOS1", "KRAS", "activates"),
        ("KRAS", "BRAF", "activates"),
        ("BRAF", "MAP2K1", "phosphorylates"), # MEK1
        ("MAP2K1", "MAPK1", "phosphorylates"),# ERK2
        ("MAPK1", "MYC", "transcription_activation"),
        
        # PI3K/AKT/MTOR Survival Axis
        ("PIK3CA", "PTEN", "antagonized_by"),
        ("PIK3CA", "AKT1", "activates"),
        ("AKT1", "MTOR", "activates"),
        ("AKT1", "BAD", "inhibits_apoptosis"),
        ("MTOR", "MYC", "translation_activation"),
        
        # Disease Phenotype Anchor Nodes
        ("MYC", "Proliferation_Phenotype", "drives"),
        ("BAD", "Apoptosis_Evasion", "drives"),
        ("AKT1", "Apoptosis_Evasion", "drives")
    ]

    for u, v, relation in interactions:
        G.add_edge(u, v, relationship=relation)

    return G


def prioritize_targets(G: nx.Graph) -> pd.DataFrame:
    """Calculates network centrality metrics to prioritize therapeutic targets.
    - Degree Centrality: Measures how many direct biological connections a node has.
    - Betweenness Centrality: Measures how often a node falls on the shortest communication
      pathway between other proteins (information bottlenecks).
    """
    degree_cent = nx.degree_centrality(G)
    between_cent = nx.betweenness_centrality(G, normalized=True)

    # Exclude phenotypic endpoint sink nodes from drug target rankings
    endpoints = {"Proliferation_Phenotype", "Apoptosis_Evasion"}
    protein_nodes = [node for node in G.nodes() if node not in endpoints]

    records = []
    for node in protein_nodes:
        records.append({
            "Gene_Target": node,
            "Degree": G.degree[node],
            "Degree_Centrality": round(degree_cent[node], 3),
            "Betweenness_Centrality": round(between_cent[node], 3),
        })

    df = pd.DataFrame(records)
    # Composite Target Priority Score: balanced combination of connectivity and bottlenecking
    df["Priority_Score"] = round((df["Degree_Centrality"] * 0.5) + (df["Betweenness_Centrality"] * 0.5), 3)
    df = df.sort_values(by="Priority_Score", ascending=False).reset_index(drop=True)
    return df


if __name__ == "__main__":
    os.makedirs("data/processed", exist_ok=True)
    kg = build_oncology_kg()
    print(f"[+] Knowledge Graph built successfully:")
    print(f"    - Total Nodes (Proteins + Phenotypes): {kg.number_of_nodes()}")
    print(f"    - Total Interactions (Edges)         : {kg.number_of_edges()}")

    target_rankings = prioritize_targets(kg)
    output_path = "data/processed/target_prioritization_rankings.csv"
    target_rankings.to_csv(output_path, index=False)
    print(f"\n[+] Top 5 Prioritized Drug Targets saved -> {output_path}")
    print(target_rankings.head(8).to_string(index=False))