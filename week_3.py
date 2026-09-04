from Motifs import *


if __name__ == "__main__":
    Dna_test = [
        "GGCGTACACGTAGACCAAACAAACACGGGACACACACACACACACACACACACACAC",
        "CAACAGGCAACTACCAACACCAACTACCAACCCACACACACACACCACACACACCAC",
        "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACC",
        "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACC"
    ]
    sample_motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGGGGG", "CACCAC"]
    sample_profile = {
        "A": [0.8, 0.0, 0.0, 0.2],
        "C": [0.0, 0.6, 0.2, 0.0],
        "G": [0.2, 0.2, 0.8, 0.0],
        "T": [0.0, 0.2, 0.0, 0.8]
    }

    print("1. Count Matrix:", Count(sample_motifs))
    print("2. Profile Matrix:", Profile(sample_motifs))
    print("3. Consensus Sequence:", Consensus(sample_motifs))
    print("4. Score:", Score(sample_motifs))
    print("5. Probability (Pr):", Pr("ACCT", sample_profile))
    print("6. Profile Most Probable K-mer:", ProfileMostProbableKmer("ACGCGTCACG", 4, sample_profile))
    print("7. Greedy Motif Search:", GreedyMotifSearch(Dna_test, 3, len(Dna_test)))