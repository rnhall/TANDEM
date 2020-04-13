#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:40:24 2020

@author: nelson
"""
import csv
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import ray
import numpy as np

class ComplexityMasker:
    """This class takes a genome (multi-sequence fasta file) or a single sequence
    as a biopython SeqRecord object and returns a .GTF interval file containing
    annotated regions of low-complexity as well as a fasta file containing the
    sequences contained within these intervals.

    Attributes:
        verbose: Used to print various helpful debuging statements.
        
        window_size: The size of the window of sequence which is
        kmerized. The window size rougly defines the resolution of your
        analysis. Smaller window sizes can identify smaller repeats but 
        consumes more resources. 
        
        step_size: The amount of nucleotides that the window slides
        over each iteration. This parameter also defines the resolution of
        your analysis. A smaller step size increases the number of samples
        taken from the genome, but it consumes more resources.
        
        cutoff: The percentage drop in complexity which is considered
        significant enough to be annotated as a region of interest.
        
        kmer_length: The kmer length that is used for all analysis.
        
        genome: Holds the dictionary object containing the genome to be
        analyzed.
        
        sequence: Holds the SeqRecord object which contains the single sequence
        to be analyzed.
    """
    def __init__(self, num_cores=8, verbose=True, window_size=10000, step_size=5000, 
                 cutoff=0.8, kmer_length=21, padding=10000):
        
        """Constructs a ComplexityFilter object"""
        self.verbose = verbose
        self.window_size = window_size
        self.step_size = step_size
        self.cutoff = cutoff
        self.kmer_length = kmer_length
        self.padding = padding
        self.genome = None
        self.sequence = None
        self.num_cores = num_cores
    
    
    def load_genome(self, filepath):
        """Loads a genome using the Biopython SeqIO module. It is stored as
        a dictionary with the .fasta sequence ID as the key and a SeqRecord
        object as the value.
        
        Args:
            filepath: A string which provides the filepath to the genome file
            to be loaded.
            
        Returns:
            None
            
        Raises:
            FileNotFoundError: The specific filepath is invalid.
        """
        print(filepath)
        try:
            self.genome = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))
        except FileNotFoundError:
            print("FileNotFoundError: the filepath you specified does not exist!")
            return
        self.genome_or_sequence = "Genome"
        if self.verbose: print("Genome loaded!")
    

    def load_sequence(self, sequence):
        """Loads a single sequence as a Biopython SeqRecord object.
        
        Args:
            sequence: The sequence to be loaded. Can be a generic string, 
            a biopython Seq object, or a SeqRecord object.
            
        Returns:
            None
            
        Raises:
            TypeError: The passed sequences was not a valid datatype.
            """
        
        if type(sequence) is SeqIO.SeqRecord:
            self.sequence = sequence
        elif type(sequence) is Seq:
            self.sequence = SeqIO.SeqRecord(sequence)
        elif type(sequence) is str:
            sequence = Seq(sequence, alphabet=generic_dna)
            self.sequence = SeqIO.SeqRecord(sequence)
        else:
            print("Sequence must be a biopython SeqRecord object!")
            raise TypeError


    def non_canonical_kmerize(self, sequence):
        """PRIVATE: Generates a (non-canonical) kmer dictionary of a specified sequence.
        
        Args:
            sequence: The sequence to be kmerized
            
        Returns:
            kmer_dict: A kmer dictionary containing each kmer as a key and the
            number of times it appears in the sequence as the value.
            
        Raises:
            None
        """
        kmer_dict = {}
        for i in range(0,len(sequence)-self.kmer_length):
            kmer = sequence[i:i+self.kmer_length]
            if kmer in kmer_dict:
                kmer_dict[str(kmer)] += 1
            else:
                kmer_dict[str(kmer)] = 1
        return kmer_dict
    
    
    def generate_kmer_trace(self, sequence):
        """Scans the provided sequence with a resolution of the provided window
        size and determines the number of unique kmers within that region.
        
        Args:
            sequence: The sequence to be analyzed
            
        Returns:
            indices, kmer_counts
            
        Raises:
            None
        """
        index = 0
        sub_sequence = sequence[index:index+self.window_size]
        kmer_counts = []
        while(index+self.window_size < len(sequence)):
            kmer_dict = self.non_canonical_kmerize(sub_sequence)
            kmer_counts.append(len(kmer_dict.keys())/self.window_size)
            index += self.step_size
            sub_sequence = sequence[index:index+self.window_size]
        return range(0, len(sequence), self.step_size)[0:len(kmer_counts)], kmer_counts
    
    @ray.remote
    def helper_generate_kmer_trace(self, sequence):
        """Scans the provided sequence with a resolution of the provided window
        size and determines the number of unique kmers within that region.
        
        Args:
            sequence: The sequence to be analyzed
            
        Returns:
            indices, kmer_counts
            
        Raises:
            None
        """
        index = 0
        sub_sequence = sequence[index:index+self.window_size]
        kmer_trace = []
        while(index+self.window_size < len(sequence)):
            kmer_dict = self.non_canonical_kmerize(sub_sequence)
            kmer_trace.append(len(kmer_dict.keys())/self.window_size)
            index += self.step_size
            sub_sequence = sequence[index:index+self.window_size]
        return range(0, len(sequence), self.step_size)[0:len(kmer_trace)], kmer_trace
    
    
    def multi_generate_kmer_trace(self, sequence):
        """Scans the provided sequence with a resolution of the provided window
        size and determines the number of unique kmers within that region.
        
        Args:
            sequence: The sequence to be analyzed
            
        Returns:
            indices, kmer_counts
            
        Raises:
            None
        """
        #ray.shutdown()
        #ray.init(num_cpus=self.num_cores)
        anchors = range(0, len(sequence), self.step_size)
        chunks = np.array_split(anchors, self.num_cores)
        kmer_traces_obj = []
        for chunk in chunks:
            try:
                kmer_traces_obj.append(self.helper_generate_kmer_trace.remote(self, sequence[chunk[0]:chunk[-1] + self.window_size + self.step_size]))
            except IndexError:
                kmer_traces_obj.append(self.helper_generate_kmer_trace.remote(self, sequence[chunk[0]:len(sequence)]))
        kmer_trace = []
        for trace in kmer_traces_obj:
            kmer_trace += ray.get(trace)[1]
        return range(0, len(sequence), self.step_size)[0:len(kmer_trace)], kmer_trace
        #ray.shutdown()
    
    def analyze_sequence(self, contig_name, sequence, export_gtf = False, export_sequences = False):
        
        if self.verbose: print("Analyzing " + contig_name)
        if len(sequence) >= self.window_size * self.num_cores * 10:
            coordinates, kmer_trace = self.multi_generate_kmer_trace(sequence)
        else:
            coordinates, kmer_trace = self.generate_kmer_trace(sequence)
        intervals = self.find_ROIs(len(sequence), coordinates, kmer_trace, self.padding)
        gtf_intervals = []
        for interval in intervals:
                gtf_intervals.append([contig_name, "Low-Complexity", "exon", interval[0], interval[1], 0, "+", ".", "gene_id\"Low_Complexity\""])
        if self.verbose: print("Completed analyzing " + contig_name)
        return gtf_intervals
    
    def analyze_genome(self, cm_fasta_filename = "Low_Complexity_Intervals.fasta", cm_gtf_filename = "Low_Complexity_Intervals.gtf"):
        
        if self.genome == None: 
            print("No genome has been loaded!") 
            return None
        #ray.shutdown()
        #ray.init(num_cpus=self.num_cores)
        gtf_intervals = []
        for contig in self.genome.keys():
            gtf_intervals.append(self.analyze_sequence(contig, self.genome[contig].seq))
        #ray.shutdown()
        self.export_as_GTF(gtf_intervals, filename = cm_gtf_filename)
        self.export_sequences(gtf_intervals, filename = cm_fasta_filename)
            
    #Code from Uvar (https://stackoverflow.com/questions/49071081/merging-overlapping-intervals-in-python)
    #Merges overlapping intervals into a single super-interval.
    def recursive_merge(self, intervals, start_index = 0):
        for i in range(start_index, len(intervals) - 1):
            if intervals[i][1] >= intervals[i+1][0]:
                new_start = intervals[i][0]
                new_end = intervals[i+1][1]
                intervals[i] = [new_start, new_end]
                del intervals[i+1]
                return self.recursive_merge(intervals.copy(), start_index=i)
        return intervals
    
    #Finds regions of interest (ROIs), which are 4 element tuples consisting of (Genomic Contig, Start, End, Sequence)
    def find_ROIs(self, length, coordinates, kmer_trace, padding):
        ROIs = []
        for index in range(0, len(coordinates)):
            if kmer_trace[index] < self.cutoff:
                if coordinates[index]-padding <= 0: 
                    ROIs.append([0, padding])
                elif coordinates[index]+padding >= length: 
                    ROIs.append([coordinates[index], length])
                else:
                    ROIs.append([coordinates[index] - padding, coordinates[index] + padding])
        return self.recursive_merge(ROIs)
    
    def search_sequences(self, ids):
        intervals = []
        for contig in ids:
            if self.verbose:
                print("Analyzing " + contig)
            coordinates, kmer_trace = self.generate_kmer_trace(self.genome[contig].seq)
            interval_list = self.find_ROIs(len(self.genome[contig].seq), coordinates, kmer_trace, 2 * self.window_size)
            for interval in interval_list:
                intervals.append([contig, "Low-Complexity", "exon", interval[0], interval[1], 0, "+", ".", "gene_id\"Low_Complexity\""])
        return intervals
    
    def search_sequence(self, contig_name):
        intervals = []
        if self.verbose:
            print("Analyzing sequence!")
        coordinates, kmer_trace = self.generate_kmer_trace(self.sequence)
        interval_list = self.find_ROIs(len(self.sequence), coordinates, kmer_trace, 2 * self.window_size)
        for interval in interval_list:
            intervals.append([contig_name, "Low-Complexity", "exon", interval[0], interval[1], 0, "+", ".", "gene_id\"Low_Complexity\""])
        return intervals
    
    def extract_sequences(self, intervals):
        sequences = []
        for contig in intervals:
            for region in  contig:
                sequence = self.genome[region[0]].seq[region[3]:region[4]]
                record = SeqRecord.SeqRecord(sequence, id=region[0] + ","+str(region[3])+","+str(region[4]))
                sequences.append(record)
        return sequences
    
    def export_sequences(self, intervals, filename):
        sequences = self.extract_sequences(intervals)
        SeqIO.write(sequences, filename, "fasta")
    
    def export_as_GTF(self, intervals, filename):
        with open(filename, 'w', newline='') as f_output:
            tsv_output = csv.writer(f_output, delimiter='\t', escapechar=' ', quoting=csv.QUOTE_NONE)
            for contig in intervals:
                for region in  contig:
                    tsv_output.writerow(region)
    