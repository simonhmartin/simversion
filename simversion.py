#!/usr/bin/env python

# A script to generate simulated genome assemblies.
# Evolves two assemblies from an ancestral sequence, with SNVs, insertions, deletions and inversions.
# Takes either an input genome in fasta format or a sequence length for a randomly generated genome.
# run python simversion.py -h for help.


import gzip, sys, argparse
import numpy as np

def break_into_contigs(sequence, n_contigs, features=[], pruning=[0,0], padding=[0,0]):
    l = len(sequence)
    max_contig = int(l/n_contigs)
    #assert max_edge >= min_edge, "Max edge cannot be shorter than mid edge."
    #assert 2*max_edge < max_contig, "Contigs are too short relative to the edges."
    contigs = []
    contig_intervals = []
    contig_pruned_lens = []
    contig_padded_lens = []
    feature_indices_by_contig = []
    features_by_contig = []
    j = 0
    for k in range(n_contigs):
        left_prune = abs(int(np.random.normal(*pruning)))
        right_prune = abs(int(np.random.normal(*pruning)))
        contig_pruned_lens.append((left_prune, right_prune))
        contig_len = max_contig - left_prune - right_prune
        j += left_prune #slide along to the new start after pruning
        c = sequence[j:(j+contig_len)]
        
        feature_indices = [i for i in range(len(features)) if features[i][0] >= j and features[i][1] <= j + contig_len]
        
        feature_indices_by_contig.append(feature_indices)
        
        _features_ = []
        
        left_padding_len = abs(int(np.random.normal(*padding)))
        right_padding_len = abs(int(np.random.normal(*padding)))
        
        for i in feature_indices:
            new_feature = features[i][:]
            new_feature[0] = new_feature[0] - j + left_padding_len
            new_feature[1] = new_feature[1] - j + right_padding_len
            _features_.append(new_feature)
        
        features_by_contig.append(_features_)
        
        left_padding_seq = "".join(np.random.choice(["a","c","g","t"], left_padding_len))
        right_padding_seq = "".join(np.random.choice(["a","c","g","t"], right_padding_len))
        contig_padded_lens.append((left_padding_len, right_padding_len))
        
        contigs.append(left_padding_seq + c + right_padding_seq)
        
        contig_intervals.append((j, j+contig_len,))
        
        j += contig_len + right_prune
    
    return (contigs, contig_intervals, features_by_contig, feature_indices_by_contig, contig_pruned_lens, contig_padded_lens)


def evolve_sequences(anc, SNV_rate, ins_rate, del_rate, inv_rate,
                     ins_mean=10, del_mean=10, inv_mean=1000,
                     ins_categories=None, del_categories = None, inv_categories=None):
    l = len(anc)
    
    rates = np.array([SNV_rate, ins_rate, del_rate, inv_rate])
    rel_rates = rates/rates.sum()
    
    #average gap between variants
    mean_gap = int(np.ceil(1/rates.sum()))
    
    var_anc_pos = []
    var_seq1_pos = []
    var_seq2_pos = []
    var_genome = []
    var_type = []
    var_len = []
    
    seq1 = ""
    seq2 = ""
    i_anc = 0
    i_seq1 = 0
    i_seq2 = 0
    while i_anc <= l:
        gap = int(np.ceil(np.random.exponential(mean_gap)))
        seq1 += anc[i_anc:(i_anc+gap)]
        seq2 += anc[i_anc:(i_anc+gap)]
        if i_anc + gap >= l: break
        i_anc += gap
        i_seq1 += gap
        i_seq2 += gap
        vartype = np.random.choice(["snv", "ins", "del", "inv"], p=rel_rates)
        
        if vartype == "snv":
            if np.random.random() > 0.5:
                #change is in seq1
                genome = 1
                seq1 += np.random.choice([l for l in "acgt" if l != anc[i_anc]])
                seq2 += anc[i_anc]
            else:
                genome = 2
                seq2 += np.random.choice([l for l in "acgt" if l != anc[i_anc]])
                seq1 += anc[i_anc]
            
            var_anc_pos.append([i_anc, i_anc+1])
            var_seq1_pos.append([i_seq1, i_seq1+1])
            var_seq2_pos.append([i_seq2, i_seq2+1])
            var_genome.append(genome)
            var_type.append("SNV")
            var_len.append(1)
            #all progress by 1
            i_anc += 1
            i_seq1 += 1
            i_seq2 += 1
        
        elif vartype == "del":
            #it's a deletion!!!!!!!!!!!!!!!!!!!!!!!
            if del_categories: del_len = int(np.random.choice(del_categories))
            else: del_len = int(np.ceil(np.random.exponential(del_mean)))
            
            _i_anc_ = i_anc + del_len
            if _i_anc_ >= l: continue # do not add this deletion as it extends beyond the end of the sequence
            
            if np.random.random() > 0.5:
                genome=1
                #deletion is in seq1, so only add to seq2
                seq2 += anc[i_anc:_i_anc_]
            else:
                genome=2
                #deletion is in seq2, so only add to seq1
                seq1 += anc[i_anc:_i_anc_]
            
            _i_seq1_ = i_seq1 + (del_len if genome == 2 else 0)
            _i_seq2_ = i_seq2 + (del_len if genome == 1 else 0)
            
            var_anc_pos.append([i_anc, _i_anc_])
            var_seq1_pos.append([i_seq1, _i_seq1_])
            var_seq2_pos.append([i_seq2, _i_seq2_])
            var_genome.append(genome)
            var_type.append("deletion")
            var_len.append(del_len)
            i_anc = _i_anc_
            i_seq1 = _i_seq1_
            i_seq2 = _i_seq2_
        
        elif vartype == "ins":
            #it's an insertion!
            if ins_categories: ins_len = int(np.random.choice(ins_categories))
            else: ins_len = int(np.ceil(np.random.exponential(ins_mean)))
            
            _i_anc_ = i_anc
            if np.random.random() > 0.5:
                genome = 1
                #insertion is in seq1, so only add to that
                seq1 += "".join(np.random.choice(["a","c","g","t"], ins_len))
            else:
                genome = 2
                #insertion is in seq2, so only add to that
                seq2 += "".join(np.random.choice(["a","c","g","t"], ins_len))
            
            _i_seq1_ = i_seq1 + (ins_len if genome == 1 else 0)
            _i_seq2_ = i_seq2 + (ins_len if genome == 2 else 0)
            
            var_anc_pos.append([i_anc, _i_anc_])
            var_seq1_pos.append([i_seq1, _i_seq1_])
            var_seq2_pos.append([i_seq2, _i_seq2_])
            var_genome.append(genome)
            var_type.append("insertion")
            var_len.append(ins_len)
            i_anc = _i_anc_
            i_seq1 = _i_seq1_
            i_seq2 = _i_seq2_
    
        elif vartype == "inv":
            #it's an inversion!!!!!!!!!!!!!!!!!!!!!!!!!!
            if inv_categories: inv_len = int(np.random.choice(inv_categories))
            else: inv_len = int(np.ceil(np.random.exponential(inv_mean)))
            
            _i_anc_ = i_anc + inv_len
            if _i_anc_ >= l: continue # do not add this inversion as it extends beyond the end of the sequence
            
            anc_subseq = anc[i_anc:_i_anc_]
            
            if np.random.random() > 0.5:
                genome = 1
                #inversion is in seq1
                seq1 += revComplement(anc_subseq).lower()
                seq2 += anc_subseq
            else:
                genome = 2
                #inversion is in seq2
                seq2 += revComplement(anc_subseq).lower()
                seq1 += anc_subseq
            
            _i_seq1_ = i_seq1 + inv_len
            _i_seq2_ = i_seq2 + inv_len
            
            var_anc_pos.append([i_anc, _i_anc_])
            var_seq1_pos.append([i_seq1, _i_seq1_])
            var_seq2_pos.append([i_seq2, _i_seq2_])
            var_genome.append(genome)
            var_type.append("inversion")
            var_len.append(inv_len)
            i_anc = _i_anc_
            i_seq1 = _i_seq1_
            i_seq2 = _i_seq2_
    
    return {"seq1":seq1, "seq2":seq2,
            "var_anc_pos":var_anc_pos, "var_seq1_pos":var_seq1_pos, "var_seq2_pos":var_seq2_pos,
            "var_genome":var_genome, "var_type":var_type, "var_len":var_len}


#subset list into smaller lists
def subset(things,subLen,asLists=False):
    starts = range(0,len(things),subLen)
    ends = [start+subLen for start in starts]
    if asLists: return [list(things[starts[i]:ends[i]]) for i in range(len(starts))]
    return [things[starts[i]:ends[i]] for i in range(len(starts))]

def fastaLineGen(names=None, seqs=None, lineLen=100):
    assert len(names) == len(seqs)
    for i in range(len(names)):
        yield ">" + names[i]
        for substr in subset("".join(seqs[i]),lineLen):
            yield substr


def parseFasta(string, makeUppercase=False):
    splitString = string.split(">")[1:]
    names = [s.split()[0] for s in splitString]
    seqs = [s[s.index("\n"):].replace("\n","").replace(" ","") for s in splitString]
    if makeUppercase: seqs = [s.upper() for s in seqs]
    return (names,seqs)

complementTrans = str.maketrans("ACGTKMRYVHBDNacgtkmryvhbdn", "TGCAMKYRBDVHNtgcamkyrbdvhn")

def revComplement(seq):
    return seq.translate(complementTrans)[::-1]

if __name__ == '__main__':

    ### parse arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta for ancestral sequence if you prefer not to use a randomly generated one.", action = "store")
    parser.add_argument("-l", "--seqLen", help="Length for randomly generated ancestral sequence", action = "store", type=int)
    parser.add_argument("-o", "--outPrefix", help="Output file prefix", action = "store", default = "sim")
    
    parser.add_argument("--SNV_rate", help="Expected SNVs per site.", action = "store", type=float, default=0.01)
    parser.add_argument("--ins_rate", help="Expected insertions per site.", action = "store", type=float, default=0.005)
    parser.add_argument("--del_rate", help="Expected deletions per site.", action = "store", type=float, default=0.005)
    parser.add_argument("--inv_rate", help="Expected inversions per site.", action = "store", type=float, default=0.0001)
    
    parser.add_argument("--ins_categories", help="Insertion len categories (means will be ignored)", action = "store", type=int, nargs="+")
    parser.add_argument("--del_categories", help="Deletion len categories (means will be ignored)", action = "store", type=int, nargs="+")
    parser.add_argument("--inv_categories", help="Inversion len categories (means will be ignored)", action = "store", type=int, nargs="+")
    
    parser.add_argument("--ins_mean", help="Insertion len is exponentially distributed with mean of INT (continuous mode).", action = "store", type=int, default = 10)
    parser.add_argument("--del_mean", help="Deletion len is exponentially distributed with mean of INT (continuous mode).", action = "store", type=int, default = 10)
    parser.add_argument("--inv_mean", help="Inversion len is exponentially distributed with mean of INT (continuous mode).", action = "store", type=int, default = 1000)
    
    parser.add_argument("--break_genome1", help="Break evolved genome 1 into contigs", action='store_true')
    parser.add_argument("--break_genome2", help="Break evolved genome 2 into contigs", action='store_true')
    
    parser.add_argument("--number_of_contigs", help="If breaking, number of contigs to break each chromosome into", action = "store", type=int, default = 3)
    
    parser.add_argument("--invertOddContigs", help="Reverse complement the first, third etc. contig in genome 2.", action='store_true')
    
    parser.add_argument("--contig_edge_pruning", help="Remove a length of sequence from edges of each contig (so these bits will be missing)", action = "store", type=int, nargs=2, default = [0,0], metavar=["mean","sd"])
    parser.add_argument("--contig_edge_padding", help="Add a length of random sequence to edges of each contig (unalignable noise)", action = "store", type=int, nargs=2, default = [0,0], metavar=["mean","sd"])
    
    args = parser.parse_args()
    
    if args.break_genome1 or args.break_genome2:
        if args.number_of_contigs <= 1:
            raise Exception("You must specify a number of contiges > 1")
    
    if args.fasta:
        print("\nLoading ancestral genome.", file=sys.stderr)
        with gzip.open(args.fasta,"rt") if args.fasta.endswith(".gz") else open(args.fasta,"rt") as fa:
            seq_names, sequences = parseFasta(fa.read())
    
    elif args.seqLen:
        print("\nGenerating ancestral genome.", file=sys.stderr)
        seq_names = ["chrom1"]
        sequences = ["".join(np.random.choice(["A","C","G","T"], args.seqLen))]
    else:
        raise Exception("Please provide a sequence length or an input fasta of your choice.")
    
genome_data_anc = {"seq_names":seq_names, "sequences": sequences, "var_pos":[], "var_type":[], "var_genome":[], "var_len":[]}

genome_data_1 = {"seq_names":seq_names, "sequences":[], "var_pos":[], "var_type":[], "var_genome":[], "var_len":[]}
genome_data_2 = {"seq_names":seq_names, "sequences":[], "var_pos":[], "var_type":[], "var_genome":[], "var_len":[]}

#evolve sequences
print("\nEvolving two descendant genomes genomes.", file=sys.stderr)
for i, sequence in enumerate(genome_data_anc["sequences"]):
    print(genome_data_anc["seq_names"][i], file=sys.stderr)
    evolved_data = evolve_sequences(sequence,
                                    SNV_rate=args.SNV_rate, ins_rate=args.ins_rate, del_rate=args.del_rate, inv_rate=args.inv_rate,
                                    ins_mean=args.ins_mean, del_mean=args.del_mean, inv_mean=args.inv_mean,
                                    ins_categories=args.ins_categories, del_categories=args.del_categories, inv_categories=args.inv_categories)
    
    genome_data_anc["var_pos"].append(evolved_data["var_anc_pos"])
    genome_data_anc["var_type"].append(evolved_data["var_type"])
    genome_data_anc["var_genome"].append(evolved_data["var_genome"])
    genome_data_anc["var_len"].append(evolved_data["var_len"])
    
    genome_data_1["sequences"].append(evolved_data["seq1"])
    genome_data_1["var_pos"].append(evolved_data["var_seq1_pos"])
    genome_data_1["var_type"].append(evolved_data["var_type"])
    genome_data_1["var_genome"].append(evolved_data["var_genome"])
    genome_data_1["var_len"].append(evolved_data["var_len"])
    
    genome_data_2["sequences"].append(evolved_data["seq2"])
    genome_data_2["var_pos"].append(evolved_data["var_seq2_pos"])
    genome_data_2["var_type"].append(evolved_data["var_type"])
    genome_data_2["var_genome"].append(evolved_data["var_genome"])
    genome_data_2["var_len"].append(evolved_data["var_len"])


#if breaking further into contigs (UNTESTED SINCE REWRITE 27 Sept 2023)
if args.break_genome1:
    print("\nBreaking genome 1 into contigs.", file=sys.stderr)
    genome_data_1B = {"seq_names":[], "sequences":[], "var_pos":[], "var_type":[], "var_genome":[], "var_len":[]}
    
    for i in range(len(genome_data_1["seq_names"])):
        contigs, contig_intervals, var_pos_by_contig, var_indices_by_contig, contig_pruned_lens, contig_padded_lens = break_into_contigs(genome_data_1["sequences"][i],
                                                                                            args.number_of_contigs,
                                                                                            genome_data_1["var_pos"][i],
                                                                                            pruning = args.contig_edge_pruning,
                                                                                            padding = args.contig_edge_padding)
        
        for j in range(len(contigs)):
            genome_data_1B["seq_names"].append(genome_data_1["seq_names"][i] + "." + str(j))
            genome_data_1B["sequences"].append(contigs[j])
            genome_data_1B["var_pos"].append(var_pos_by_contig[j])
            genome_data_1B["var_type"].append([genome_data_1["var_type"][i][k] for k in var_indices_by_contig[j]])
            genome_data_1B["var_genome"].append([genome_data_1["var_genome"][i][k] for k in var_indices_by_contig[j]])
            genome_data_1B["var_len"].append([genome_data_1["var_len"][i][k] for k in var_indices_by_contig[j]])
    
    genome_data_1 = genome_data_1B

if args.break_genome2:
    print("\nBreaking genome 2 into contigs.", file=sys.stderr)
    genome_data_2B = {"seq_names":[], "sequences":[], "var_pos":[], "var_type":[], "var_genome":[], "var_len":[]}
    
    for i in range(len(genome_data_2["seq_names"])):
        contigs, contig_intervals, var_pos_by_contig,
        var_indices_by_contig, contig_pruned_lens, contig_padded_lens = break_into_contigs(genome_data_2["sequences"][i],
                                                                                            args.number_of_contigs,
                                                                                            genome_data_2["var_pos"][i],
                                                                                            pruning = args.contig_edge_pruning,
                                                                                            padding = args.contig_edge_padding)
        
        for j in range(len(contigs)):
            genome_data_2B["seq_names"].append(genome_data_2["seq_names"][i] + "." + str(j+1))
            genome_data_2B["sequences"].append(contigs[j])
            genome_data_2B["var_pos"].append(var_pos_by_contig[j])
            genome_data_2B["var_type"].append([genome_data_2["var_type"][i][k] for k in var_indices_by_contig[j]])
            genome_data_2B["var_genome"].append([genome_data_2["var_genome"][i][k] for k in var_indices_by_contig[j]])
            genome_data_2B["var_len"].append([genome_data_2["var_len"][i][k] for k in var_indices_by_contig[j]])
    
    genome_data_2 = genome_data_2B

### Invert every other contig in the second genome (to simulate sequencing of reverse strand)
if args.invertOddContigs:
    print("\nInverting odd contigs in genome 2.", file=sys.stderr)
    for i in range(0, len(genome_data_2["seq_names"]), 2):
        genome_data_2["sequences"][i] = revComplement(genome_data_2["sequences"][i])
        l = genome_data_2["sequences"][i]
        genome_data_2["var_pos"][i] = [[l-var[1], l-var[0]] for var in reversed(genome_data_2["var_pos"][i])]
        #you also need to reverse the other info columns!
        genome_data_2["var_len"][i] = [var for var in reversed(genome_data_2["var_len"][i])]
        genome_data_2["var_genome"][i] = [var for var in reversed(genome_data_2["var_genome"][i])]
        genome_data_2["var_type"][i] = [var for var in reversed(genome_data_2["var_type"][i])]

### Write outputs
print("\nWriting output files", file=sys.stderr)

#if a random sequence was used, export it
if not args.fasta:
    with gzip.open(args.outPrefix + ".anc.fa.gz", "wt") as fa_anc:
        for line in fastaLineGen(names=genome_data_anc["seq_names"], seqs=genome_data_anc["sequences"], lineLen=100):
            fa_anc.write(line + "\n")

with gzip.open(args.outPrefix + ".1.fa.gz", "wt") as fa1:
    for line in fastaLineGen(names=genome_data_1["seq_names"], seqs=genome_data_1["sequences"], lineLen=100):
        fa1.write(line + "\n")

with gzip.open(args.outPrefix + ".2.fa.gz", "wt") as fa2:
    for line in fastaLineGen(names=genome_data_2["seq_names"], seqs=genome_data_2["sequences"], lineLen=100):
        fa2.write(line + "\n")

with gzip.open(args.outPrefix + ".1.variants.bed.gz", "wt") as bed1:
    for i in range(len(genome_data_1["seq_names"])):
        for j in range(len(genome_data_1["var_pos"][i])):
            bed1.write("\t".join([genome_data_1["seq_names"][i],
                                    str(genome_data_1["var_pos"][i][j][0]),
                                    str(genome_data_1["var_pos"][i][j][1]),
                                    genome_data_1["var_type"][i][j],
                                    str(genome_data_1["var_genome"][i][j]),
                                    str(genome_data_1["var_len"][i][j])]) + "\n")

with gzip.open(args.outPrefix + ".2.variants.bed.gz", "wt") as bed2:
    for i in range(len(genome_data_2["seq_names"])):
        for j in range(len(genome_data_2["var_pos"][i])):
            bed2.write("\t".join([genome_data_2["seq_names"][i],
                                    str(genome_data_2["var_pos"][i][j][0]),
                                    str(genome_data_2["var_pos"][i][j][1]),
                                    genome_data_2["var_type"][i][j],
                                    str(genome_data_2["var_genome"][i][j]),
                                    str(genome_data_2["var_len"][i][j])]) + "\n")

with gzip.open(args.outPrefix + ".anc.variants.bed.gz", "wt") as bed2:
    for i in range(len(genome_data_anc["seq_names"])):
        for j in range(len(genome_data_anc["var_pos"][i])):
            bed2.write("\t".join([genome_data_anc["seq_names"][i],
                                    str(genome_data_anc["var_pos"][i][j][0]),
                                    str(genome_data_anc["var_pos"][i][j][1]),
                                    genome_data_anc["var_type"][i][j],
                                    str(genome_data_anc["var_genome"][i][j]),
                                    str(genome_data_anc["var_len"][i][j])]) + "\n")
