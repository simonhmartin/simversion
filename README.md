# simversion: Simulate genome assemblies with variants

* Evolve two assemblies from an ancestral sequence, with SNVs and SVs (insertions, deletions and inversions).
* Takes either an input genome in fasta format or a sequence length for a randomly generated genome.
* Variants are poisson distributed in the genome and will never overlap.
* simversion is useful for testing variant detection tools, but it is not designed as an accurae simulator or molecular evoluton.
* SV lengths can be either drawn from a random exponential distribution with predefined mean, or from a set of discrete length categories.
* Run `python simversion.py -h` for help.
