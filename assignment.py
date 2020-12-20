from itertools import product, combinations, chain
from collections import Counter

translation_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                     'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                     'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                     'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                     'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                     'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                     'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                     'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                     'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                     'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                     'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                     'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                     'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                     'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                     'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                     'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
                     }

# nomenclature for degenerate codons
expanded_code = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
                 'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'],
                 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                 'N': ['A', 'C', 'G', 'T']
                 }

# Precomputing inverse lookup for expanded code
inv_expanded_code = {str(sorted(v)): k for k, v in expanded_code.items()}

# Caching initial lookup results for faster calculations
global pos_dict
pos_dict = {}
global comb_dict
comb_dict = {}

# helpful for validating input
valid_nucleotides = 'ACGTWSMKRYBDHVN'
valid_aa = 'GAVLIMFWPSTCYNQDEKRH*'


def check_valid_aa(amino_acids):
    valid_aa_set = set(list(valid_aa))
    if len(amino_acids) == 0:
        print("*************invalid input***************")
        exit(1)

    for i in amino_acids:
        if i not in valid_aa_set:
            print("*************invalid input***************")
            exit(1)


def update_pos_dict_get_codon_list(amino_acids):
    triplet_list = []

    for a in amino_acids:
        pos_dict[a] = {1: set({}), 2: set({}), 3: set({})}
        for k, v in translation_table.items():
            if v == a:
                pos_dict[a][1].add(k[0])
                pos_dict[a][2].add(k[1])
                pos_dict[a][3].add(k[2])
                triplet_list.append(k)
    return triplet_list


def update_comb_dict(amino_acids):
    for k in amino_acids:
        if k not in comb_dict.keys():
            comb_dict[k] = {1: set({}), 2: set({}), 3: set({})}
            for i in range(1, 4):
                for j in range(1, len(pos_dict[k][i]) + 1):
                    for item in combinations(pos_dict[k][i], j):
                        comb_dict[k][i].add(item)


def get_efficiency(code_triplet_list, list_of_triplet_dict):
    count_observed_triplet = {}
    for p in code_triplet_list:
        total_number_of_codons = len(expanded_code[p[0]]) * len(expanded_code[p[1]]) * len(expanded_code[p[2]])
        count_observed_triplet[''.join(p)] = len(
            (set(list_of_triplet_dict[p[0]][0]).intersection(set(list_of_triplet_dict[p[1]][1]))).intersection(
                set(list_of_triplet_dict[p[2]][2]))) / total_number_of_codons

    return count_observed_triplet


def get_valid_triplet(expanded_base, triplet_list, pos):
    result = []
    for j in triplet_list:
        if j[pos] in expanded_code[expanded_base]:
            result.append(j)

    return result


def get_cartesian(comb_dict, pos, amino_acids):
    tmp_set = []
    result_set = []
    for k in amino_acids:
        if not tmp_set:
            tmp_set = (list(comb_dict[k][pos]))
            tmp_set = [' '.join(tups) for tups in tmp_set]
        else:
            comb = list(product(tmp_set, [' '.join(tups) for tups in list(comb_dict[k][pos])]))
            comb = [' '.join(tups) for tups in comb]
            tmp_set = comb

    for ele in tmp_set:
        spl = set((ele.split()))
        result_set.append(list(spl))

    result = []
    for i in result_set:
        if set(i) not in result:
            result.append(set(i))

    result_set = [list(i) for i in result]

    #print(result_set)
    return result_set


def get_all_possible_codon_combinations(amino_acids):
    update_comb_dict(amino_acids)
    set1 = get_cartesian(comb_dict, 1, amino_acids)
    set2 = get_cartesian(comb_dict, 2, amino_acids)
    set3 = get_cartesian(comb_dict, 3, amino_acids)
    all_possible_combinations = list(product(set1, set2, set3))
    #print(all_possible_combinations)
    return all_possible_combinations


def get_expanded_codon_list(all_possible_combinations):
    expanded_code_result_list = []
    for i in all_possible_combinations:
        res = []
        for j in i:
            code = inv_expanded_code[str(sorted(j))]
            res.append(code)
        expanded_code_result_list.append(res)
    return expanded_code_result_list


def get_expanded_base_codon_map(expanded_code_result_list, codon_list):
    expanded_base_codon_map = {}
    for p in expanded_code_result_list:
        for q in range(len(p)):
            if p[q] not in expanded_base_codon_map.keys():
                expanded_base_codon_map[p[q]] = {}
                expanded_base_codon_map[p[q]][q] = get_valid_triplet(p[q], codon_list, q)

            elif q not in expanded_base_codon_map[p[q]].keys():
                expanded_base_codon_map[p[q]][q] = get_valid_triplet(p[q], codon_list, q)

    return expanded_base_codon_map


def get_codon_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set, float
        returns two values the set of most efficient codons for the input set list, e.g.
        {'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'} and the achieved efficiency e.g. 0.75
    """
    check_valid_aa(amino_acids)
    codon_list = update_pos_dict_get_codon_list(amino_acids)
    all_possible_combinations = get_all_possible_codon_combinations(amino_acids)
    expanded_code_result_list = get_expanded_codon_list(all_possible_combinations)
    expanded_base_codon_map = get_expanded_base_codon_map(expanded_code_result_list, codon_list)

    count_observed_triplet = get_efficiency(expanded_code_result_list, expanded_base_codon_map)
    itemMaxValue = max(count_observed_triplet.items(), key=lambda x: x[1])
    setOfKeys = set({})
    for key, value in count_observed_triplet.items():
        if value == itemMaxValue[1]:
            setOfKeys.add(key)
    return setOfKeys, itemMaxValue[1]


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))


def truncate_list_of_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set
        the set of sets of amino acids that can be coded with 100% efficiency, i.e. {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
    """
    check_valid_aa(amino_acids)
    pset = ([set(i) for i in powerset(amino_acids)])

    answer_set = set({})
    for i in reversed(pset):
        (a, b) = get_codon_for_amino_acids(i)
        if b == 1.0:
            if len(answer_set) == 0:
                answer_set.add(frozenset(i))
                l = len(i)
            elif len(i) < l:
                break
            elif len(i) == l:
                answer_set.add(frozenset(i))

    return answer_set


if __name__ == "__main__":
    # using sets instead of lists throughout the code since the oder doesn't matter and all items should be unique
    assert get_codon_for_amino_acids({'A', 'I', 'V'}) == ({'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'}, 0.75)
    assert get_codon_for_amino_acids({'M', 'F'}) == ({'WTS', 'WTK', "WTB"}, 0.5)

    # "frozenset" here since this seems to be the only way to get a set of sets - see https://stackoverflow.com/questions/5931291/how-can-i-create-a-set-of-sets-in-python
    assert truncate_list_of_amino_acids({'A', 'V', 'I'}) == {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
