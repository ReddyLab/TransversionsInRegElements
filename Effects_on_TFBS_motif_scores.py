from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio import Seq

f = open('jaspar_Ts_Tv.tab', 'w')
f.write("name\tpos\tAG\tCT\tAC\tAT\tCG\tGT\n")
with open("pfm_all.txt") as handle:
 for m in motifs.parse(handle, "jaspar"):
    counts = m.counts
    cons = m.consensus
    pssm = m.counts.normalize(pseudocounts=0.5).log_odds()
    cons_score = pssm.calculate(cons)
    cons_list = list(cons)
    cons_str =  str(cons)
    for i, c in enumerate(cons_list):
      new_cons_str_A = Seq.Seq("".join((cons_str[0:i], "A", cons_str[i+1:])), IUPAC.unambiguous_dna)
      new_cons_str_C = Seq.Seq("".join((cons_str[0:i], "C", cons_str[i+1:])), IUPAC.unambiguous_dna)
      new_cons_str_G = Seq.Seq("".join((cons_str[0:i], "G", cons_str[i+1:])), IUPAC.unambiguous_dna)
      new_cons_str_T = Seq.Seq("".join((cons_str[0:i], "T", cons_str[i+1:])), IUPAC.unambiguous_dna)
      new_score_A = pssm.calculate(new_cons_str_A)
      new_score_C = pssm.calculate(new_cons_str_C)
      new_score_G = pssm.calculate(new_cons_str_G)
      new_score_T = pssm.calculate(new_cons_str_T)
      Ts_deltas = (abs(new_score_A-new_score_G), abs(new_score_C-new_score_T))
      Tv_deltas = (abs(new_score_A-new_score_C), abs(new_score_A-new_score_T), 
                   abs(new_score_C-new_score_G), abs(new_score_T-new_score_G))
      f.write(m.name)
      f.write("\t")
      central_distance = abs(0.5 - float(i)/len(counts[1,:]))
      f.write(str(central_distance))
      f.write("\t")
      f.write("\t".join(map(str, Ts_deltas + Tv_deltas)))
      f.write("\n")

f.close()
