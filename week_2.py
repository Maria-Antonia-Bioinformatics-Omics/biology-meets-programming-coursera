# File Replication

import sys

genome_test = "CATGGGTATCGATACTACCACCAACATCATACACCGACCACACACACACCGACCACACCCACACACACACTACCACCACACACACACCACACCACACACACCACACCAC"

# FasterSymbolArray: Tracks symbol frequencies across circular genomes using an optimized sliding window algorithm.
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i] + 1
    return array


def PatternCount(Pattern, Text):
    count = 0 
    pattern_length = len(Pattern)
    text_length = len(Text)
    for i in range(text_length - pattern_length + 1):
        if Text[i:i+pattern_length] == Pattern:
            count += 1          
    return count

# MinimumSkew: Identifies genomic positions where the GC skew reaches its minimum value, locating candidate oriC sites.
def MinimumSkew(Genome):
    positions = []
    skew_raw = SkewArray(Genome)
    min_value = min(skew_raw)
    
    for i in range(len(skew_raw)):
        if skew_raw[i] == min_value:
            positions.append(i)
            
    return positions

# SkewArray: Calculates the running difference between Guanine and Cytosine counts (G - C) across the genome.
def SkewArray(Genome):
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == "G":
            skew.append(skew[-1] + 1)
        elif Genome[i] == "C":
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

# ApproximatePatternCount: Computes the total occurrences of a motif in a genome within a d-mismatch threshold.
def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1
    return count

# HammingDistance: Counts the number of point mutations (mismatches) between two equal-length DNA sequences.
def HammingDistance(p, q):
    distance=0
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

# ApproximatePatternMatching: Maps starting positions of a motif within text, allowing up to d mismatches.
def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

# Prints Results

if __name__ == "__main__":
    genome_sample = "CATGGGTATCGATACTACCACCAACATCATACACCGACCACACACACACCGACCACACCCACACACACACTACCACCACACACACACCACACCACACACACCACACCAC"
    
    print("1. Skew Array:", SkewArray(genome_sample[:20]))
    print("2. Minimum Skew:", MinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGTATTCCCAATCTCCACGCATAGGCACTCTACACTACACCTCAACTCCACACTCTA"))
    print("3. Hamming Distance:", HammingDistance("GGGCGACACTCC", "GGGACTACCCGT"))
    print("4. Approximate Pattern Matching:", ApproximatePatternMatching("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", "ATTCTGGA", 3))
    print("5. Approximate Pattern Count:", ApproximatePatternCount("AAAAA", "AACAAGCTGATAAACACTCACTACTACCAACTACCAACTACCAACACCAACTACCAACCCACACACACACACCACACACACCAC", 2))
    print("6. Faster Symbol Array:", FasterSymbolArray("AAAAGGGG", "A"))