from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq

from fuzzywuzzy import fuzz


p1 = Seq('(((EM)|(I+))|((L+)|[^WSQECAPKRIGVYHMFL]))|(([^LWFPHQNTAEIGC]|(FP))|((K+)[YGARHFDVQIE]))')
p2 = Seq('(((K|M)|(R+))|((W|Q)|[WRN]))|(((R+)|C{2})|((K+)[KQLGSWMP]))')

# ([SP][AHP]{5})|(PH)

# alignement of the 2 patterns
alignments = pairwise2.align.globalxx(p1, p2)
print(format_alignment(*alignments[0]))
best_alignment = alignments[0]
identity = (best_alignment[2] / len(best_alignment[0])) * 100


print(identity)


x = fuzz.ratio(p1, p2)
print(x)
