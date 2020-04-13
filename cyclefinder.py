#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 11:21:55 2020

@author: nelson
"""

import graph_tool as gt
from graph_tool.all import *
import random
import numpy as np


class CycleFinder:
    
    def __init__(self, k=21, verbose=False, reads_mode=False):
        self.k = k
        self.graphs = []
        self.verbose = verbose
        self.reads_mode = reads_mode
        
    
    def generate_kmer_transitions(self, sequence, k=21):
    
        kmer_transitions = {}
        for i in range(0,len(sequence)-k-1):
            current_kmer = sequence[i:i+k]
            next_kmer = sequence[i+1:i+k+1]
            if current_kmer+">"+next_kmer in kmer_transitions.keys():
                kmer_transitions[current_kmer+">"+next_kmer] += 1
            else:
                kmer_transitions[current_kmer+">"+next_kmer] = 1  
        print("Dictionary contains " + str(len(kmer_transitions.keys())) + " unique kmer transitions!")
        return kmer_transitions

    
    def generate_kmer_dict(self, sequence, k=21):
    
        kmer_dict = {}
        for i in range(0,len(sequence)-k):
            current_kmer = sequence[i:i+k]
            if current_kmer in kmer_dict.keys():
                kmer_dict[current_kmer] += 1
            else:
               kmer_dict[current_kmer] = 1  
        if self.verbose: print("Dictionary contains " + str(len(kmer_dict.keys())) + " unique kmers!")
        return kmer_dict
        
    def generate_graph_from_sequence(self, sequence, k=21, filter_singles=False):
        
        graph = gt.Graph(directed=True)
        kmer_dict = self.generate_kmer_dict(sequence, k)
        kmer_transitions = self.generate_kmer_transitions(sequence, k)
        
        kmer_transitions_copy = kmer_transitions.copy()
        if filter_singles:
            for kmer_transition in kmer_transitions_copy.keys():
                if kmer_transitions[kmer_transition] < 2:
                    del kmer_transitions[kmer_transition]

        vert_prop_kmer = graph.new_vertex_property("string")
        vert_prop_counts = graph.new_vertex_property("int")
        edge_prop_counts = graph.new_edge_property("float")

        for kmer in list(kmer_dict.keys()):
            v = graph.add_vertex()
            vert_prop_kmer[v] = kmer
            vert_prop_counts[v] = kmer_dict[kmer]

        for transition in list(kmer_transitions.keys()):
            trans = transition.split(">")
            current_kmer = trans[0]
            next_kmer = trans[1]
            v1 = gt.util.find_vertex(graph, vert_prop_kmer, current_kmer)[0]
            v2 = gt.util.find_vertex(graph, vert_prop_kmer, next_kmer)[0]
            e1 = graph.add_edge(v1, v2)
            edge_prop_counts[e1] = kmer_transitions[transition]
    
        graph.vertex_properties["kmer"] = vert_prop_kmer 
        graph.vertex_properties["counts"] = vert_prop_counts
        graph.edge_properties["counts"] = edge_prop_counts
        
        if self.verbose: print("Graph generated!")
        return graph
    
    def negate_edge_weights(self, graph):
        
        for vertex in graph.vertices():
            if graph.get_out_degrees([vertex]) > 1:
                total = 0
                for e in vertex.out_edges():
                    total += graph.edge_properties["counts"][e]
                for e in vertex.out_edges():
                    graph.edge_properties["counts"][e] /= -np.log(graph.edge_properties["counts"][e]/total)
            else:
                for e in vertex.out_edges():
                    graph.edge_properties["counts"][e] = 0
                    
        if self.verbose: print("Edge weights negated!")
        return graph
    
    def prune_graph(self, graph):
        
        in_degrees = np.array(graph.get_in_degrees(graph.get_vertices()))
        out_degrees = np.array(graph.get_out_degrees(graph.get_vertices()))
        total_degrees = in_degrees + out_degrees
        while 0 in total_degrees or 1 in total_degrees:
            graph.remove_vertex(np.where(total_degrees < 2))
            in_degrees = np.array(graph.get_in_degrees(graph.get_vertices()))
            out_degrees = np.array(graph.get_out_degrees(graph.get_vertices()))
            total_degrees = in_degrees + out_degrees 
        
        if self.verbose: print("Graph pruned!")
        return graph
    
    def split_into_subgraphs(self, graph):
        
        sub_graph_property_map, num_sub_graphs = gt.topology.label_components(graph)
        sub_graphs = []
        for i in range(0, len(num_sub_graphs)):
            sub_graph_view = gt.GraphView(graph, vfilt=sub_graph_property_map.a == i)
            sub_graph = gt.Graph(sub_graph_view, prune=True)
            if sub_graph.num_vertices() > 1:
                sub_graphs.append(sub_graph)
        
        if self.verbose: print("Graph split into subgraphs!")
        return sub_graphs
    
    # def extract_consensus(self, graph):
        
    #     #if self.verbose: print(graph.num_vertices())
    #     sequences = []
    #     for path in gt.topology.all_shortest_paths(graph, 1, 0, weights=graph.edge_properties["counts"]):
    #         sequence = ""
    #         for vertex in path:
    #             sequence += graph.vertex_properties["kmer"][vertex][-1]
    #         sequences.append(sequence)
            
    #     if self.verbose: print("Consensus sequences extracted!")
    #     #return sequences / Should normally only return one consensus sequence per graph
    #     return sequences

    def extract_consensus(self, graph, num_samples=20):
        
        paths = []
        num_vertices = []
        for i in range(num_samples):
            source = random.randint(1, graph.num_vertices()-1)
            path = gt.topology.shortest_path(graph, source, source-1, weights=graph.edge_properties["counts"])
            paths.append(path[0])
            num_vertices.append(len(path[0]))
                
        sequence = ""
        for vertex in paths[num_vertices.index(max(num_vertices))]:
            sequence += graph.vertex_properties["kmer"][vertex][-1]
            
        return sequence
        
    
    def plot_graph(self, graphs):
        
        for graph in graphs:
            gt.draw.graph_draw(graph, vertex_size=graph.vertex_properties["counts"])
            
    def get_graph_and_consensus(self, graphs):
        
        graph_and_consensus = []
        for graph in graphs:
            consensus = self.extract_consensus(graph)
            graph_and_consensus.append((graph, consensus))
        return graph_and_consensus
    
    def get_intervals_from_graph(self, k, sequence, consensus, graph, offset = 0):
        
        #Determine where graph-kmers are found in the reference sequence
        kmers = [graph.vertex_properties["kmer"][i] for i in range(0, graph.num_vertices())]
        positions = []
        indices = []
        meta_data = sequence.id.split(",")
        contig = meta_data[0]
        seq_start = int(meta_data[1])
        sequence = sequence.seq
        
        for i in range(0, len(sequence)):
            current_kmer = sequence[i:i+k]
            if current_kmer in kmers:
                positions.append(1)
                indices.append(i)
            else:
                positions.append(0)
        
        #Merge annotations
        cutoff = len(consensus)*2
        for i in range(0, len(indices)-1):
            current_count = indices[i]
            next_count = indices[i+1]
            if next_count - current_count <= cutoff:
                for i in range(current_count, next_count):
                    positions[i] = 1
        
        #Extracts intervals where merged kmers mapped above
        positions = np.array(positions)
        indices = np.where(positions == 1)[0]
        intervals = []
        start = indices[0]
        for i in range(1, len(indices)-1):
            if indices[i+1] - indices[i] > 1:
                intervals.append([start + offset, indices[i] + offset])
                start = indices[i+1]
        intervals.append([start + offset, indices[i] + offset])
        
        if self.verbose: print("Intervals generated!")       
        return intervals
                
        
        
            
