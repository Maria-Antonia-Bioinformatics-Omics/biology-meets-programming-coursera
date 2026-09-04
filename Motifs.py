# Count: Builds a frequency matrix counting occurrences of each nucleotide per position in motifs.
def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Profile: Calculates the probability matrix by dividing nucleotide counts by the total number of motifs.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)
    return profile

# Consensus: Forms a representative DNA string using the most frequent nucleotide per column.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Score: Measures total motif variation by counting mismatches against the consensus sequence.
def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    t = len(Motifs)
    k = len(consensus)
    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score

# Pr: Computes the overall probability of a k-mer based on a given profile matrix.
def Pr(text, profile):
    p = 1
    t = len(text)
    for i in range(t):
        p = p * profile[text[i]][i]
    return p

# ProfileMostProbableKmer: Finds the most likely k-mer within a text sequence using the profile matrix.
def ProfileMostProbableKmer(text, k, profile):
    max_prob = -1.0
    most_probable = text[0:k]
    
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = Pr(kmer, profile)
        if prob > max_prob:
            max_prob = prob
            most_probable = kmer
            
    return most_probable

# GreedyMotifSearch: Implements a greedy heuristic search algorithm to locate candidate regulatory motifs across DNA sequences.
def GreedyMotifSearch(Dna, k, t):
    best_motifs = [Dna[i][0:k] for i in range(t)]
    n = len(Dna[0])
    
    for i in range(n - k + 1):
        motifs = [Dna[0][i:i + k]]
        for j in range(1, t):
            profile = Profile(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
            
    return best_motifs