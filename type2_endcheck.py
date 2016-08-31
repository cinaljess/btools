import sys
import difflib
from string import maketrans

revcomplement_memo = {'A':'T'}
revcompTBL = maketrans("AGCTagctWSKMYRnN", "TCGAtcgaWSMKTYnN")
def revcomplement(seq):
    """
    revcomplement(seq)
    A quick reverse-complement routine that memo-izes queries, understands
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


def align(end, rcend):
    matcher = difflib.SequenceMatcher(a = end, b = rcend)
    match_score = matcher.find_longest_match(0, len(matcher.a), 0, len(matcher.b))
    return match_score


def main(ends_list):
    ends_list = ends_list.split(",")
    ends_list = [x.strip(" ") for x in ends_list]
    ncs = 'ATGC'
    rc_list = []
    silly_list = []
    for end in ends_list:
        for c in end:
            if c not in ncs:
                ter = ends_list.index(end)
                silly_list = ends_list.pop(ter)
                break

    notgood = False
    sim_list = []
    rc_list = [revcomplement(x) for x in ends_list]
    for x in ends_list:
        for d in rc_list:
            score = align(x, d)
            if score.size >= 3:
                notgood = True
                sim_list.append((x, revcomplement(d)))
    if not notgood:            
        print("Good to go!!!")
        if silly_list:
            print "Bad entry: ", silly_list
    else:
        print("Not good!")
        if silly_list:
            print "Bad entry: ", silly_list
        for x in sim_list:
            print(x)

if __name__ == "__main__":
    main(sys.argv[1])
