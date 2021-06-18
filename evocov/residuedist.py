from Bio.PDB import PDBParser

# create parser
parser = PDBParser()

# read structure from file
structure = parser.get_structure('SPK', '7dk3.pdb')

model = structure[0]
chain = model['A']

# this example uses only the first residue of a single chain.
# it is easy to extend this to multiple chains and residues.
for residue1 in chain:
    for residue2 in chain:
        if residue1 != residue2:
            # compute distance between CA atoms
            try:
                distance = residue1['CA'] - residue2['CA']
            except KeyError:
                ## no CA atom, e.g. for H_NAG
                continue
            if distance < 6:
                print(residue1, residue2, distance)
        # stop after first residue
        break
