import math

E6_GROUPS = {
    'A': 'aliphatic',
    'C': 'aliphatic',
    'I': 'aliphatic',
    'L': 'aliphatic',
    'M': 'aliphatic',
    'V': 'aliphatic',
    'F': 'aromatic',
    'H': 'aromatic',
    'W': 'aromatic',
    'Y': 'aromatic',
    'N': 'polar',
    'Q': 'polar',
    'S': 'polar',
    'T': 'polar',
    'K': 'positive',
    'R': 'positive',
    'D': 'negative',
    'E': 'negative',
    'G': 'special',
    'P': 'special'
}

#test_counts2 = {'A': 5, 'I': 3, 'L': 2, 'K': 3, 'D': 2, 'G': 1, 'F': 4}
#sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLLN"
def calculate_e6(amino_acid_count):
    group_totals = {}
    for aa, count in amino_acid_count.items():
        if aa not in E6_GROUPS:
            continue
        group = E6_GROUPS[aa]
        if group not in group_totals:
            group_totals[group] = count
        else:
            group_totals[group] = group_totals[group] + count
    total = sum(group_totals.values())

    e6 = 0
    for group, count in group_totals.items():
        p = count / total
        e6 -= p * math.log2(p)
    return e6



def calculate_e20(amino_acid_count):
    total = sum(amino_acid_count.values())

    e20 = 0
    for aa, count in amino_acid_count.items():
        if aa not in E6_GROUPS:
            continue
        p = count / total
        e20 -= p * math.log2(p)
    return e20


def calculate_vkbat(predictions):
    helix_count = predictions.count("Helix")
    other_count = predictions.count("Other")
    sheet_count = predictions.count("Sheet")
    #disordered_count = predictions.count("Disordered")

    hos_list = [helix_count, other_count, sheet_count] #Add disorded_count

    max_val = 0
    k = 0
    for item in hos_list:
        if item > 0:
            k += 1
        if item > max_val:
            max_val = item
    if max_val == 0:
        return 0.0
    return ((k / max_val) * len(predictions))

#test_preds = ['Helix', 'Sheet', 'Other', 'Helix', 'Sheet', 'Helix']
#print(calculate_vkbat(test_preds))







