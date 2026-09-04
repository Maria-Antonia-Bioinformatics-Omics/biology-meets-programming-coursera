# Biology Meets Programming — UC San Diego (Coursera)

Code, exercises, and notes from UC San Diego's algorithmic bioinformatics introductory course. The main goal here was to understand the computational logic behind genomic sequence processing and implement the algorithms from scratch.

### Key Concepts Covered
* **Pattern Matching & k-mers:** Algorithms for finding frequent motifs and counting subsequences in DNA sequences.
* **Replication Origin (*oriC*):** Locating bacterial replication origins using skew diagrams and Hamming distance (handling mutations/mismatches).
* **Regulatory Motif Finding:** Implementing profile matrix searches to identify transcription factor binding sites.

### Languages & Skills
* **Language:** Python
* **Focus:** String manipulation, basic data structures, and algorithmic logic applied to biological problems.

### 📂 Repository Structure
* `week_1_hidden_messages/` — Pattern counting and k-mer identification scripts.
* `week_2_replication_origin/` — Skew array calculations and mismatch search algorithms.
* `week_3_4_motif_finding/` — Algorithms for regulatory motif discovery.

--------------------------------------------------------------------------------------------------------------------------------------

## 🧬 Week 1: Finding Hidden Messages in DNA (oriC)

### Biological Context
Cellular replication depends on the precise identification of the **origin of replication (*oriC*)**, where initiator proteins (such as DnaA) bind to short, repeated DNA sequences called **motifs** (or *k-mers*). Computationally identifying these frequent repetitions allows us to map functional genomic regions without relying solely on expensive wet-lab experimentation.

### Implemented Algorithms (*week_1.py*)

* **`PatternCount(Text, Pattern)`**: Calculates the exact occurrence count of a specific motif within a genomic sequence.
* **`FrequencyMap(Text, k)`**: Generates a complete frequency distribution profile for all *k-mers* of length *k*.
* **`FrequentWords(Text, k)`**: Filters and returns the most frequent *k-mers* (candidate protein-binding sites).
* **`Reverse(Pattern)`**: Inverts the DNA sequence orientation (3' -> 5' or 5' -> 3').
* **`Complement(Pattern)`**: Applies Watson-Crick base-pairing rules (A <-> T, C <-> G).
* **`ReverseComplement(Pattern)`**: Models the opposing strand of the DNA double helix for bidirectional pattern searches.
* **`PatternMatching(Pattern, Genome)`**: Returns the exact starting positions (0-indexed) of a motif across a genome.

--------------------------------------------------------------------------------------------------------------------------------------


## 🧬 Week 2: Finding Replication Origins in Bacterial Genomes (Skew & Mismatches)

### Biological Context
DNA replication is asymmetrical: one strand (leading) is synthesized continuously, while the other (lagging) is synthesized in fragments. This leads to a deamination process where Cytosine (C) mutates into Thymine (T), causing a shortage of C on the single-stranded leading strand. Measuring the imbalance between Guanine and Cytosine (**GC Skew**) helps locate the exact origin of replication (*oriC*) where the skew reaches its minimum value. Additionally, real binding sites contain evolutionary mutations, requiring algorithms tolerant to **mismatches** (Hamming Distance).

### Implemented Algorithms (Replication.py / week_2.py)

* **`SkewArray(Genome)`**: Computes the running difference between Guanine (G) and Cytosine (C) counts across the genome.
* **`MinimumSkew(Genome)`**: Finds the genomic positions where the GC skew reaches its minimum value (identifying candidate *oriC* locations).
* **`HammingDistance(p, q)`**: Calculates the number of point mutations (mismatches) between two equal-length DNA sequences.
* **`ApproximatePatternMatching(Text, Pattern, d)`**: Locates starting positions of a motif allowing up to d mismatches.
* **`ApproximatePatternCount(Pattern, Text, d)`**: Counts total occurrences of a motif within a genome given a mismatch threshold d.
* **`FasterSymbolArray(Genome, symbol)`**: Uses a sliding window mechanism to track symbol frequencies across circular genomes efficiently.


--------------------------------------------------------------------------------------------------------------------------------------

## 🧬 Week 3: Motif Finding (Regulatory Motifs & Greedy Search)

### Biological Context
Transcription factors bind to specific DNA patterns called **motifs** to regulate gene expression. Unlike exact sequence matches, motifs often exhibit evolutionary variations across different genes or organisms. To identify these unknown binding sites, we construct frequency and probability matrices (**Count** and **Profile**) to determine the most representative sequence (**Consensus**) and evaluate overall alignment quality (**Score**).

### Implemented Algorithms (Motifs.py / week_3.py)

* **`Count(Motifs)`**: Constructs a matrix counting nucleotide occurrences across motif alignments.
* **`Profile(Motifs)`**: Generates a position probability matrix by normalizing nucleotide counts.
* **`Consensus(Motifs)`**: Derives the most frequent nucleotide string across aligned motifs.
* **`Score(Motifs)`**: Calculates total mismatches between individual motifs and the consensus sequence.
* **`Pr(text, profile)`**: Computes the probability of a specific k-mer based on the profile matrix.
* **`ProfileMostProbableKmer(text, k, profile)`**: Scans a sequence to identify the most probable k-mer for a given profile.
* **`GreedyMotifSearch(Dna, k, t)`**: Uses a greedy heuristic to search for optimal regulatory motifs across multiple DNA sequences.