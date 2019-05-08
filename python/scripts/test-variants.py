#!/usr/bin/env python

import MOODS.scan
import MOODS.tools
import MOODS.parsers

from operator import itemgetter

mat1 = [
       [1,1,1,1,1], # A
       [0,0,0,0,0], # C
       [0,0,0,0,0], # G
       [0,0,0,0,0], # T
      ]
      
mat2 = [
       [1,1,1,1,1,1,1], # A
       [0,0,0,0,0,0,0], # C
       [0,0,0,0,0,0,0], # G
       [0,0,0,0,0,0,0], # T
      ]
      
      
matrices = [mat1,mat2]
thresholds = [5,7]
bg = MOODS.tools.flat_bg(4)

scanner = MOODS.scan.Scanner(7)
scanner.set_motifs(matrices, bg, thresholds)


# test sequence
# W is either A or T (standard IUPAC symbol)
seq = "WWWWWWWWWWWWWWW"
#      01234567890123456789012345


snps = MOODS.tools.snp_variants(seq)

indels = []
indels.append((4,5,""))
indels.append((5,6,""))
indels.append((4,5,"AAA"))
indels.append((5,8,"AAA"))

variants = list(snps) + [MOODS.tools.variant(start, end, replacement) for (start, end, replacement) in indels]

results_variants = scanner.variant_matches(seq,variants)

errors = 0

for i in range(len(results_variants)):
    print("---------")
    print("Matrix {}".format(i))
    for hit in results_variants[i]:
        
        hit_variants = sorted([(variants[j].start_pos, variants[j].end_pos, variants[j].modified_seq) for j in hit.variants], key=itemgetter(0,1))

        modified_seq = list(seq.lower())
        hit_indicator = list(hit.pos*" " + len(matrices[i][0])*"^")
            
        offset = 0
        for var in hit_variants:
            if var[2] != "":
                modified_seq[var[0] + offset:var[1] + offset] = list(var[2])
                offset = offset - (var[1] - var[0]) + len(var[2])
            else:
                modified_seq[var[0]+offset:var[1]+offset] = ["-"]
                hit_indicator.insert(var[0]+offset, "-")
                offset = offset - (var[1] - var[0]) + 1
            # print(offset)
        modified_seq = "".join(modified_seq)
        hit_indicator = "".join(hit_indicator)
        
        for k in range(len(hit_indicator)):
            if (modified_seq[k] != "A" and hit_indicator[k] == "^") or (hit_indicator[k] == "-" and modified_seq[k] != "-") or (hit_indicator[k] != "-" and modified_seq[k] == "-"):
                print("Error?")
                errors = errors + 1
                break
    
        print("Position: {}; Score: {}; Variants: {}".format(str(hit.pos), str(hit.score), str(hit_variants)))
        print(modified_seq)
        print(hit_indicator)
    
    print("Total matches: "+str(len(results_variants[i])))

print("Spotted {} error(s)".format(errors))