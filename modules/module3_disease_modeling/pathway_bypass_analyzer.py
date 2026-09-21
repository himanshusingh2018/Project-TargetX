"""
Module 3: Disease Modeling
Script: pathway_bypass_analyzer.py
Simulates target knockout and identifies oncogenic bypass escape pathways
to recommend rational combination therapeutic targets.
"""

import networkx as nx
from disease_knowledge_graph import build_oncology_kg


def find_bypass_routes(G: nx.Graph, source: str, target: str, blocked_node: str) -> list[list[str]]:
    """Identifies alternative communication paths between source and target
    when the primary target (blocked_node) is inhibited.
    """
    # Create a copy of the graph simulating the drug-inhibited state
    G_inhibited = G.copy()
    if blocked_node in G_inhibited:
        G_inhibited.remove_node(blocked_node)

    # Find all simple paths (length <= 5) that bypass the blocked target
    try:
        bypass_paths = list(nx.all_simple_paths(G_inhibited, source=source, target=target, cutoff=5))
        return bypass_paths
    except nx.NetworkXNoPath:
        return []


def analyze_combination_strategy():
    G = build_oncology_kg()
    primary_target = "EGFR"
    driver_source = "ERBB2"
    disease_sink = "MYC"

    print(f"[+] Simulating monotherapy target blockade: {primary_target}")
    
    # Check baseline paths from receptor upstream to cancer sink
    baseline_paths = list(nx.all_simple_paths(G, source=driver_source, target=disease_sink, cutoff=5))
    print(f"[+] Total baseline signaling routes ({driver_source} -> {disease_sink}): {len(baseline_paths)}")

    # Check remaining routes when EGFR is knocked out
    bypass_routes = find_bypass_routes(G, source=driver_source, target=disease_sink, blocked_node=primary_target)
    
    print(f"\n--- Bypass Resistance Analysis ({primary_target} Inhibited) ---")
    if bypass_routes:
        print(f"[!] Warning: {len(bypass_routes)} alternative escape route(s) detected:")
        for idx, path in enumerate(bypass_routes, 1):
            print(f"    Route {idx}: {' -> '.join(path)}")
        
        # Identify non-source, non-sink nodes in the escape routes as combo targets
        escape_nodes = set()
        for path in bypass_routes:
            escape_nodes.update(path[1:-1])
        
        print(f"\n[+] Recommended Rational Combination Targets: {sorted(list(escape_nodes))}")
    else:
        print("[+] Pathway fully shut down. No bypass routes detected within cutoff.")


if __name__ == "__main__":
    analyze_combination_strategy()