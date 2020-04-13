#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 13:02:51 2020

@author: nelson
"""

from complexitymasker import ComplexityMasker
from cyclefinder import CycleFinder
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonWarning
import csv
import ray
import argparse
import warnings

class Tandem:
    
    def __init__(self, reads_mode=False):
        self.reads_mode = reads_mode
        self.complexity_filter_params = []
    
    @ray.remote
    def cyclefinder_graph_generation(self, kmer_len, sequence):
        g = CycleFinder(verbose=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore') 
            graph = g.generate_graph_from_sequence(sequence.seq, k=kmer_len, filter_singles=True)
            graph = g.negate_edge_weights(graph)
            graph = g.prune_graph(graph)
            subs = g.split_into_subgraphs(graph)
            graphs_and_consensus = g.get_graph_and_consensus(subs)
        return graphs_and_consensus
    
    @ray.remote
    def cycle_finder_interval_curation(self, kmer_len, sequence, consensus, graph, offset = 0):
        g = CycleFinder(verbose=True)
        intervals = g.get_intervals_from_graph(kmer_len, sequence, consensus, graph, offset)
        return intervals
    
    def find_tandem_repeats_in_genome(self, genome_filepath, size_cutoff = 50, length_cutoff = 100, kmer_len = 21,
                                      window_size = 10000, step_size = 2500, cutoff = 0.8, padding = 20000,  
                                      complexity_filter = True, verbose = True,
                                      cm_fasta_filename = "Low_Complexity_Intervals.fasta", cm_gtf_filename = "Low_Complexity_Intervals.gtf",
                                      cf_fasta_filename = "Consensus_Repeats.fasta", cf_gtf_filename = "Tandem_Repeat_Intervals.gtf"):
        
        if complexity_filter:
            #Determines regions of low-complexity and saves them as sequences and
            #as intervals
            cm = ComplexityMasker(verbose=verbose)
            cm.window_size = window_size
            cm.step_size = step_size
            cm.cutoff = cutoff
            cm.kmer_length = kmer_len
            cm.padding = padding
            cm.load_genome(genome_filepath)
            cm.analyze_genome(cm_fasta_filename, cm_gtf_filename)
        
        sequences = []
        if complexity_filter:
            sequences = list(SeqIO.parse(cm_fasta_filename, "fasta"))
        else:
            sequences = list(SeqIO.parse(genome_filepath, "fasta"))
        
        graphs_and_consensus = []
        for sequence in sequences:
            print("Analyzing " + sequence.id.split(",")[0])
            graphs_and_consensus.append(self.cyclefinder_graph_generation.remote(self, kmer_len, sequence))
        
        intervals_obj = []
        consensus_sequences = []
        meta_data = []
        for i, gandc in enumerate(graphs_and_consensus, 0):
            graphs = ray.get(gandc)
            for graph in graphs:
                consensus = graph[1]
                if graph[0].num_vertices() > size_cutoff:
                    consensus_sequences.append(graph[1])
                    start = int(sequences[i].id.split(",")[1])
                    intervals_obj.append(self.cycle_finder_interval_curation.remote(self, kmer_len, sequences[i], graph[1], graph[0], start))
                    meta_data.append([sequences[i].id.split(",")[0], consensus])
        
            
        gtf_intervals = []
        seq_records = []
        for x, intervals in enumerate(intervals_obj, 0):
            intervals = ray.get(intervals)
            for i, interval in enumerate(intervals, 0):
                start = int(intervals[i][0])
                end = int(intervals[i][1])
                if end - start > length_cutoff:
                    gtf_intervals.append([meta_data[x][0], "Low Complexity", "Repeat", start, end, 0, "+", ".", meta_data[x][1]])
                    record = SeqRecord(Seq(meta_data[x][1]), id=meta_data[x][0] + ":" + str(start) + "," + str(end), description = str(end - start))
                    seq_records.append(record)
                    
        with open(cf_gtf_filename, 'w', newline='') as f_output:
              tsv_output = csv.writer(f_output, delimiter='\t', escapechar=' ', quoting=csv.QUOTE_NONE)
              for entry in gtf_intervals:
                  tsv_output.writerow(entry)

        SeqIO.write(seq_records, cf_fasta_filename, "fasta") 
        
        unique_consensus = []
        lengths = []
        ids = []
        for record in seq_records:
            if record.seq not in unique_consensus:
                unique_consensus.append(record.seq)
                ids.append(record.id)
                lengths.append(int(record.description))
            else:
                index = unique_consensus.index(record.seq)
                ids[index] += "; " + record.id
                lengths[index] += int(record.description)
        
        unique_consensus_records = []
        for i, consensus in enumerate(unique_consensus, 0):
            record = SeqRecord(consensus, id=ids[i], description=str(lengths[i]))
            unique_consensus_records.append(record)
            
        SeqIO.write(unique_consensus_records, "Unique_Consensus_Sequences.fasta", "fasta") 
        
        print("Done! Thanks for using TANDEM :)")



def __main__():
    parser = argparse.ArgumentParser(description='Finds tandem repeats in genomes using a graph-based repeat assember.')
    parser.add_argument('-genome', metavar='--genome', type=str,
                    help='the path to the genome you want to analyze', required=True)
    parser.add_argument('-size_cutoff', metavar='--size_cutoff', type=int, default=50,
                    help='the minimum length of tandem repeat monomer you wish to annotate', required=False)
    parser.add_argument('-length_cutoff', metavar='--length_cutoff', type=int, default=200,
                    help='the minimum length of tandem repeat locus you wish to annotate', required=False)
    parser.add_argument('-k', metavar='--kmer_length', type=int, default=21,
                    help='the kmer size used for all data processing', required=False)
    parser.add_argument('-window_size', metavar='--window_size', type=int, default=10000,
                    help='the batch size of nuclotides analyzed in complexity filter', required=False)
    parser.add_argument('-step_size', metavar='--step_size', type=int, default=2500,
                    help='the increment by which the window is moved along the genome for each iteration of complexiy filter', required=False)
    parser.add_argument('-cutoff', metavar='--cutoff', type=int, default=0.8,
                    help='the complexity cutoff value below which genome regions are further analyzed', required=False)
    parser.add_argument('-padding', metavar='--padding', type=int, default=20000,
                    help='extra nucleotides from the left and right of a region designated for further anlaysis', required=False)
    parser.add_argument('-complexity_filter', metavar='--complexity_filter', type=bool, default=True,
                    help='determines whether the sequences are complexity filtered before being analyzed. Caution, only set to false if the sequences are small in length!', required=False)
    parser.add_argument('-verbose', metavar='--verbose', type=int, default=True,
                    help='output extra debugging statements', required=False)
    parser.add_argument('-cm_fasta_filename', metavar='--cm_fasta_filename', type=str, default="Low_Complexity_Intervals.fasta",
                    help='filename for the complexity masked genome regions (as a .fasta file)', required=False)
    parser.add_argument('-cm_gtf_filename', metavar='--cm_gtf_filename', type=str, default="Low_Complexity_Intervals.gtf",
                    help='filename for the complexity masked genome intervals (as a .gtf file)', required=False)
    parser.add_argument('-cf_fasta_filename', metavar='--cf_fasta_filename', type=str, default="Consensus_Repeats.fasta",
                    help='filename of annotated tandem repeats (as a .fasta file)', required=False)
    parser.add_argument('-cf_gtf_filename', metavar='--cf_gtf_filename', type=str, default="Tandem_Repeat_Intervals.gtf",
                    help='filename of annotated tandem repeat intervals (as a .gtf file)', required=False)
    parser.add_argument('-v', action='version', version='TANDEM v0.1')

    args = vars(parser.parse_args())
    
    genome_filepath = args['genome']
    size_cutoff = args['size_cutoff']
    length_cutoff = args['length_cutoff']
    kmer_length = args['k']
    window_size = args['window_size']
    step_size = args['step_size']
    cutoff = args['cutoff']
    padding = args['padding']
    complexity_filter = args['complexity_filter']
    verbose = args['verbose']
    cm_fasta_filename = args['cm_fasta_filename']
    cm_gtf_filename = args['cm_gtf_filename']
    cf_fasta_filename = args['cf_fasta_filename']
    cf_gtf_filename = args['cf_gtf_filename']
    
    ray.init()
    try:
        tan = Tandem()
        tan.find_tandem_repeats_in_genome(genome_filepath, size_cutoff, length_cutoff, kmer_length, window_size, step_size,
                                          cutoff, padding, complexity_filter, verbose, cm_fasta_filename, cm_gtf_filename,
                                          cf_fasta_filename, cf_gtf_filename)
    except:
        ray.shutdown()
    ray.shutdown()
    
if __name__ == '__main__':
    __main__()
        
