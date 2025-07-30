#!/usr/bin/env python3
"""
Draw phylogenetic trees from STRdust log files.

This script parses a STRdust log file, extracts tree information,
and generates a visualization of the hierarchical structure.

Usage:
  python draw_tree.py <log_file>
"""

import sys
import os
import argparse
import pygraphviz as pgv

def parse_log_and_draw_trees(logfile):
    """Parse STRdust log file and draw hierarchical trees."""
    trees = {}
    seq_weights = {}

    print(f"Parsing log file: {logfile}", file=sys.stderr)
    
    with open(logfile) as f:
        for line in f:
            if "DEBUG STRdust::phase_insertions] " in line and "with dissimilarity" in line:
                # Extract the name of the tree
                name = line.split(' ')[3].rstrip(':')
                
                # Initialize tree if not already created
                if name not in trees:
                    trees[name] = pgv.AGraph(directed=True)
                    trees[name].node_attr['style'] = 'filled'
                    seq_weights[name] = {}
                
                # Extract node and children information
                node = line.split(' ')[5]
                children = (line.split(' ')[11], line.split(' ')[14])
                seqs = (line.split(' ')[12], line.split(' ')[15].rstrip('\n'))
                
                # Add edges to the tree
                trees[name].add_edge(node, children[0])
                trees[name].add_edge(node, children[1])
                
                # Record sequence lengths for color weighting
                for child, seq in zip(children, seqs):
                    if len(seq) > 2:
                        seq_weights[name][child] = len(seq)
            
            if "Roots for this tree:" in line:
                name = line.split(' ')[3].rstrip(':')
                # Highlight root nodes
                for root in line.split(':')[-1].replace(' ', '').replace('\n', '').replace('[', '').replace(']', '').split(','):
                    if root:  # Skip empty strings
                        n = trees[name].get_node(root)
                        n.attr['fillcolor'] = "#CCCCFF"
    
    # Apply colors based on sequence lengths and generate output
    if not trees:
        print("No trees found in the log file.", file=sys.stderr)
        return
    
    print(f"Found {len(trees)} trees in the log file.", file=sys.stderr)
    
    for label, G in trees.items():
        # Skip empty trees or trees with no weights
        if not seq_weights.get(label):
            print(f"Skipping tree {label}: no sequence weights", file=sys.stderr)
            continue
            
        maxlength = max(seq_weights[label].values())
        for node, length in seq_weights[label].items():
            color = int((length / maxlength) * 255)
            n = G.get_node(node)
            n.attr['fillcolor'] = f"#{255:02x}{255 - color:02x}{255 - color:02x}"
        
        G.layout(prog='dot')
        
        # Output DOT format to stdout
        print(f"# Tree: {label}")
        print(G.string())

def main():
    parser = argparse.ArgumentParser(description="Generate tree visualizations from STRdust log files")
    parser.add_argument("logfile", help="Path to the STRdust log file")
    parser.add_argument("--format", choices=["dot", "png"], default="dot", 
                        help="Output format (default: dot format to stdout)")
    parser.add_argument("--output-dir", help="Directory to save PNG files (for --format=png)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.logfile):
        print(f"Error: Log file '{args.logfile}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if args.format == "png" and not args.output_dir:
        print("Error: --output-dir must be specified when using --format=png", file=sys.stderr)
        sys.exit(1)
    
    if args.format == "png" and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    if args.format == "dot":
        # For DOT format, we print to stdout
        parse_log_and_draw_trees(args.logfile)
    else:
        # For PNG format, we save to files
        trees = {}
        seq_weights = {}
        
        with open(args.logfile) as f:
            for line in f:
                if "DEBUG STRdust::phase_insertions] " in line and "with dissimilarity" in line:
                    name = line.split(' ')[3].rstrip(':')
                    if name not in trees:
                        trees[name] = pgv.AGraph(directed=True)
                        trees[name].node_attr['style'] = 'filled'
                        seq_weights[name] = {}
                    node = line.split(' ')[5]
                    children = (line.split(' ')[11], line.split(' ')[14])
                    seqs = (line.split(' ')[12], line.split(' ')[15].rstrip('\n'))
                    trees[name].add_edge(node, children[0])
                    trees[name].add_edge(node, children[1])
                    for child, seq in zip(children, seqs):
                        if len(seq) > 2:
                            seq_weights[name][child] = len(seq)
                if "Roots for this tree:" in line:
                    name = line.split(' ')[3].rstrip(':')
                    for root in line.split(':')[-1].replace(' ', '').replace('\n', '').replace('[', '').replace(']', '').split(','):
                        if root:
                            n = trees[name].get_node(root)
                            n.attr['fillcolor'] = "#CCCCFF"
                            
        for label, G in trees.items():
            if not seq_weights.get(label):
                print(f"Skipping tree {label}: no sequence weights", file=sys.stderr)
                continue
                
            maxlength = max(seq_weights[label].values())
            for node, length in seq_weights[label].items():
                color = int((length / maxlength) * 255)
                n = G.get_node(node)
                n.attr['fillcolor'] = f"#{255:02x}{255 - color:02x}{255 - color:02x}"
            
            G.layout(prog='dot')
            output_file = os.path.join(args.output_dir, f"graph_{label}.png")
            G.draw(output_file)
            print(f"Generated: {output_file}", file=sys.stderr)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Operation cancelled by user.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)