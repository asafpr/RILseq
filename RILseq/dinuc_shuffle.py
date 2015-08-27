from collections import defaultdict, Counter
import random

def difreq(rseq):
    counts = defaultdict(lambda: defaultdict(int))
    for a, b in zip(rseq, rseq[1:]):
        counts[a][b] += 1
    return dict((k, dict(v)) for k,v in counts.iteritems())

# <http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice>
weighted_choice = lambda s : random.choice(sum(([v]*wt for v,wt in s),[]))

def shuffle_difreq(rseq):
    freqs = difreq(rseq)

    # get the first base by the total frequency across the sequence
    shuff_seq = [weighted_choice(Counter(rseq).items())]

    for i in range(1, len(rseq)):
        # each following base is based of the frequency of the previous base
        # and their co-occurence in the original sequence.
        try:
            shuff_seq.append(weighted_choice(freqs[shuff_seq[-1]].items()))
        except KeyError:
            # If the nt is never the first of di-nuc (only one appearance at
            # the end of the sequence, choose the next nt as the first one
            shuff_seq.append(weighted_choice(Counter(rseq).items()))

    return "".join(shuff_seq)

if __name__ == "__main__":
    seq = 'AACAAAACAAAAAGTTGGGTAGTTTGAGAAC' * 20
    print difreq(seq)
    shuff_seqr = shuffle_difreq(seq)
    print shuff_seqr
    print difreq(shuff_seqr)
