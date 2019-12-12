import pandas as pd

bases = ["T", "C", "A", "G"]
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon_table = dict(zip(codons, amino_acids))

codon_data = []

for codon in codon_table:

    aa = codon_table[codon]

    syn_cnt = 0
    non_syn_cnt = 0

    #now do 9 mutations for each codon
    for i in range(0,3):
        base = codon[i]
        bases = ["T", "C", "A", "G"]
        alt_bases = [b for b in bases if b is not base]

        for base in alt_bases:

            new_codon_list = list(codon)
            new_codon_list[i] = base
            new_aa = codon_table["".join(new_codon_list)]

            if aa != new_aa:
                non_syn_cnt += 1
            else:
                syn_cnt += 1

    fraction_S_N = syn_cnt / non_syn_cnt
    codon_data.append([codon, codon_table[codon], syn_cnt, non_syn_cnt, fraction_S_N])

    print("codon: {}\taa: {}\tnon_syn_cnt: {}\tsyn_cnt: {}\tfraction_S_N: {}".format(codon, aa, non_syn_cnt, syn_cnt, round(fraction_S_N,4)))

df = pd.DataFrame(codon_data, columns=['codon', 'aa', 'syn', 'non_syn', 'fraction_S_N'])
df = df.round(decimals=4)
df.to_csv("codon_syn_non_syn_probabilities.txt", index=False)
