'''
Analyis of PROTACtor run and output
'''

import os

#Protein protein docking analysis
def pp_compatibility():
    """
    asses compatibiltiy of of protein complexes by lookiing at the rate of pocket to pocket binds compared to other outputs
    """
    score = []
    for file in os.listdir('.'):
        if file.startswith('complex.') and file.endswith('.pdb'):
            num = file.split('.')[1]
            score.append(num)
    avg = (sum(score) / len(score)) / 1000
    return avg


def min_max(text):
    """
    Extract the minimum and maximum values from a text file containing distances.
    :param text: Path to the text file.
    """
    min_val = float('inf')
    max_val = float('-inf')
    with open(text, 'r') as f:
        for line in f:
            if ":" not in line:
                continue
            try:
                # extract value after the colon, before "Å"
                value = float(line.split(":")[1].split()[0])
                min_val = min(min_val, value)
                max_val = max(max_val, value)
            except Exception as e:
                print(f"Skipping line: {line.strip()} — {e}")
    return [min_val, max_val]


def get_lysines():
    with open('lysines.txt', 'r') as f:
        out = f.readlines()
        return out


def main():
    with open('output.txt', 'w') as f:
        f.write(f'Analysis of PROTACtor Output\n'
                f'#Protein-Protein Docking\n'
                f'  compatibility score: {pp_compatibility()}\n'
                f'  Min and Max distance values: {min_max('lig_distances.txt')}\n'
                f'  Accessable Surface Lysines: {get_lysines}')
        
    


main()