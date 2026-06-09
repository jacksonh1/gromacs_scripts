import fragforge.utils.pyrosetta_tools as pyrosetta_tools
import fragforge.utils.structure.bptools as bptools

s = bptools.basic.import_structure_biopython('./helix_fusion-V2.pdb')
seq = bptools.basic.extract_sequences_from_structure(s)['A']
new_seq_1 = ''.join(['A' for i in seq])
new_seq_2 = ''.join(['S' for i in seq])

pyrosetta_tools.thread_and_relax_full_structure('./helix_fusion.pdb', 'A', new_seq=new_seq_1, out_pdb='./helix_fusion.pdb')
pyrosetta_tools.thread_and_relax_full_structure('./helix_fusion.pdb', 'A', new_seq=new_seq_2, out_pdb='./helix_fusion-V2.pdb')


