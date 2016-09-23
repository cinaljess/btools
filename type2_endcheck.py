"""
Module to check whether TypeIIs end-sites are compatible. Looks for
3bp homology then 2 basepair edge-homology for each input sequence given.
Comparisons are made between each element in the list and to the reverse
complement of each element. Repeat elements are also validated along with
noncanonical basepairs.

Example:
    python2.7 type2_endcheck.py "GAGG, GAGG, TACT, GACT, PAPP"
"""
import sys
from string import maketrans
from collections import Counter


def main(ends_list):
    ends_list = ends_list.split(",")
    ends_list = [x.strip(" ") for x in ends_list]
    ncs = 'ATGC'
    silly_list = []
    for end in ends_list:
        for c in end:
            if c not in ncs:
                ter = ends_list.index(end)
                silly_list = ends_list.pop(ter)
                break
    notgood = False
    sim_list = set([])
    rc_list = [revcomplement(x) for x in ends_list]
    counts = Counter(ends_list)
    self_list = set([])
    # Check list for repeats
    for c, n in counts.items():
        if n >= 2:
            notgood = True
            self_list.add((c))
    for x in ends_list:
        # Validate no ends share homology to each other
        for g in ends_list:
            if g != x:
                score = align(x, g)
                if score >= 3:
                    notgood = True
                    sim_list.add((x, g))
        # Validate no reverse complements are equivalent to entry list
        if x in rc_list:
            notgood = True
            idx = rc_list.index(x)
            sim_list.add((x, rc_list[idx]))
        # Validate no ends share 3 max homology & 2bp edge homology of revers complement list
        for h in rc_list:
            revscore = align(x, h)
            if revscore >= 3:
                rrevset = [h, reverse_region(h)]
                for p in rrevset:
                    rpositionscore = align(x[:2], p[:2])
                    if rpositionscore == 2:
                        notgood = True
                        sim_list.add((x, p))
    if not notgood:
        print('Good to go!!!')
        if silly_list:
            print 'Bad entry: ', silly_list
    else:
        print('Not good!')
        if silly_list:
            print 'Bad entry: ', silly_list
        for x in sim_list:
            print 'Entry: ' + str(x[0]) + ' > (' + revcomplement(x[0]) + ') : ' + revcomplement(x[1]) + ' >  (' + \
                  reverse_region(x[1]) + ')'
        for x in self_list:
            print 'Entry: ' + x + ' appeared more than once'


def revcomplement(seq):
    """
    A quick reverse-complement routine that understands
    IUPAC ambiguity codes, and preserves case.
    """
    revcompTBL = maketrans('AGCTagctWSKMYRnN', 'TCGAtcgaWSMKTYnN')
    _t = list(seq.translate(revcompTBL))
    _t.reverse()
    rc = ''.join(_t)
    return rc


def reverse_region(region):
    return region[::-1]


def align(end, rcend):
    pairs = zip(end, rcend)
    match_score = 0
    for a, b in pairs:
        if a == b:
            match_score += 1
    return match_score


if __name__ == "__main__":
    main(sys.argv[1])
