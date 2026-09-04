# week_1_hidden_messages/pattern_matching.py
# Here, we will set Text equal to the oriC and the trinucleotide to the Pattern

Text = "ACGTTGCATGCGGTGCGCGCATGATGCGGCGCGATGAGCGGAGCT"
Pattern = "GCG"
Genome = "AACTCGCGTATACCTCGCGCTTTGCGTTGCGGTCGAATTGCGTGTGTGATTT"
k = 3

# PatternCount: Counts how many times a specific sequence appears in a DNA strand to measure motif density.
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

# FrequencyMap: Maps every k-mer to its occurrence count, building a full frequency profile of the sequence.
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freq:
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1
    return freq

# FrequentWords: Extracts only the most frequent k-mers ($k$-mers with peak occurrences) to highlight potential biological signals over noise.
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

# Reverse: Inverts the DNA sequence orientation (3' -> 5' or 5' -> 3').
def Reverse(Pattern):
    rev = ""
    for char in Pattern:
        rev = char + rev  # Coloca o caractere atual no início da string rev
    return rev

# Complement: Replaces each nucleotide with its Watson-Crick base pair (A<->T, C<->G)
def Complement(Pattern):
    comp = ""
    for char in Pattern:
        if char == "A":
            comp += "T"
        elif char == "T":
            comp += "A"
        elif char == "C":
            comp += "G"
        elif char == "G":
            comp += "C"
    return comp


# ReverseComplement: Generates the matching opposite strand of DNA (5 -> 3'), essential because algorithms must search both strands for biological signals.
def ReverseComplement(Pattern):
    rev_comp = ""
    for char in Pattern:
        if char == "A":
            rev_comp = "T" + rev_comp
        elif char == "T":
            rev_comp = "A" + rev_comp
        elif char == "C":
            rev_comp = "G" + rev_comp
        elif char == "G":
            rev_comp = "C" + rev_comp
    return rev_comp

# PatternMatching: Returns the exact starting positions (0-based indices) of a motif to locate where specific elements occur in the genome.
def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions

# Printing the results:
if __name__ == "__main__":
    Pattern = "ACG"
    Genome = "AAACCCGGT"
    print("--- Week 1 ---")
    print(f"PatternCount: {PatternCount('GCGCG', 'GCG')}")
    print(f"FrequencyMap: {FrequencyMap('ACGTTGCATGTCACGTTTTGCGTCACTACACAGGACGTT', 3)}")
    print(f"FrequentWords: {FrequentWords('ACGTTGCATGTCACGTTTTGCGTCACTACACAGGACGTT', 3)}")
    print(f"Reverse: {Reverse(Pattern)}")
    print(f"Complement: {Complement(Pattern)}")
    print(f"ReverseComplement: {ReverseComplement('AAAACCCGGT')}")
    print(f"PatternMatching: {PatternMatching(Pattern, Genome)}")