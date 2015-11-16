from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio import Seq

#
# Download the most up-to-date version of the Jaspar nonredundant pfms.
import urllib
urllib.urlretrieve ("http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt", "pfm_all.txt")

# There is a bug in Bio::motifs that throws an error if there are any blank lines in the jaspar database.
# This is a workaround to get rid of those lines in the downloaded pfm file.
#
fixed_pfm_file = open("pfm_all.fixed.txt", "w")

with open("pfm_all.txt") as f:
  for line in f.readlines():
    if line.strip():
      fixed_pfm_file.write(line)

fixed_pfm_file.close()

# Output is printed to stdout to enable pipes to other process, etc. The output gives the effect of each possible mutation 
# for each jaspar pfm. The file is tab delimited, with a header to make for easy reading into R or other downstream analyses.
# Each line of the output has the following fields:
#
# 1)   name -- the name of the pfm in JASPAR
# 2)   pos  -- the relative position within the matrix. The value is from -0.5 to 0.5, where 0 is the center of the motif.
# 3-8) NN   -- the change in pssm score associated with each possible mutation at that position in the motif.
#
print "name\tpos\tAG\tCT\tAC\tAT\tCG\tGT"
with open("pfm_all.fixed.txt") as handle:
 for m in motifs.parse(handle, "jaspar"):

#
# Get the counts and the consensus motif for the pfm
#
    counts = m.counts
    cons = m.consensus

#
# convert to pssm, adding a pseudocount of 0.5 to each base.
#
    pssm = m.counts.normalize(pseudocounts=0.5).log_odds()
    cons_score = pssm.calculate(cons)
    cons_list = list(cons)
    cons_str =  str(cons)

# 
# for each position, generate a new test sequence for each possible nucleotide
# at that position. Then score that test sequence relative to the original pssm.
# Next, evaluate the absolute value of the score difference between every pair of
# test sequence and classify each pair as either a transition or transversion.
#
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
      central_distance = abs(0.5 - float(i)/len(counts[1,:]))

      print "%(name)s\t%(pos)f\t%(AG)f\t%(CT)f\t%(AC)f\t%(AT)f\t%(CG)f\t%(GT)f" % \
            {'name': m.name, 'pos':  central_distance,  \
             'AG': abs(new_score_A-new_score_G),  \
             'CT': abs(new_score_C-new_score_T),  \
             'AC': abs(new_score_A-new_score_C),    \
             'AT': abs(new_score_A-new_score_T),  \
             'CG': abs(new_score_C-new_score_G),  \
             'GT': abs(new_score_G-new_score_T)}
