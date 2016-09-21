"""
Module to check whether TypeIIs end-sites are compatible. Looks for
greater than 2 basepair similarity for each input sequence given.
Comparisons are made between each element in the list and to the reverse
complement of each element.
"""
import sys
import difflib
from string import maketrans


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
    for x in ends_list:
        for g in ends_list:
            if g != x:
                selfscore = align(x, g)
                if selfscore.size >= 3:
                    revset = [g, reverse_region(g)]
                    for k in revset:
                        positionscore = align(x, k[:3])
                    if positionscore.size == 2:
                        notgood = True
                        sim_list.add((x, g))
        for h in rc_list:
            revscore = align(x, h)
            if revscore.size >= 3:
                rrevset = [h, reverse_region(h)]
                for p in rrevset:
                    rpositionscore = align(x, p[:3])
                    if rpositionscore.size == 2:
                        notgood = True
                        sim_list.add((x, h))
    if not notgood:
        print('Good to go!!!')
        if silly_list:
            print 'Bad entry: ', silly_list
    else:
        print('Not good!')
        if silly_list:
            print 'Bad entry: ', silly_list
        for x in sim_list:
            print 'Entry: ' + str(x[0]) + ' > ' + str(x[1]) + '(' + revcomplement(x[1]) + ')'


revcomplement_memo = {'A':'T'}
revcompTBL = maketrans('AGCTagctWSKMYRnN', 'TCGAtcgaWSMKTYnN')
def revcomplement(seq):
    """
    revcomplement(seq)
    A quick reverse-complement routine that understands
    IUPAC ambiguity codes, and preserves case.
    """
    global revcomplement_memo
    try:
        rc = revcomplement_memo[seq]
    except KeyError:
        _t = list(seq.translate(revcompTBL))
        _t.reverse()
        rc = ''.join(_t)
        revcomplement_memo[seq] = rc
        revcomplement_memo[rc]  = seq
    return(rc)


def reverse_region(region):
    return region[::-1]


def align(end, rcend):
    matcher = difflib.SequenceMatcher(a = end, b = rcend)
    match_score = matcher.find_longest_match(0, len(matcher.a), 0, len(matcher.b))
    return match_score


if __name__ == "__main__":
    main(sys.argv[1])
