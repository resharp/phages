import matplotlib.pyplot as plt
import statistics as stat

bases = ["T", "C", "A", "G"]
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon_table = dict(zip(codons, amino_acids))

codon_syn_list = []

for codon in codon_table:

    aa = codon_table[codon]

    syncnt = 0
    nonsyncnt = 0

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
                # print("nonsyn: {}".format(new_aa))
                nonsyncnt += 1
            else:
                # print("syn: {}".format(new_aa))
                syncnt += 1

    fS_N = syncnt / nonsyncnt
    codon_syn_list += [fS_N]

    print("codon: {} aa: {} nonsyn: {} syn: {} fS_N: {}".format(codon, aa, nonsyncnt, syncnt, fS_N))

    print(stat.mean(codon_syn_list))

print(len(codon_table))

# aas = [a for a in amino_acids]

plt.plot(codons, codon_syn_list)
# plt.clf()
# plt.hist(codon_syn_list)
plt.show()