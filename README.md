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
