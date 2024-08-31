import oxli

def test_simple():
    cg = oxli.KmerCountTable()
    kmer = "ATCG"

    assert cg.get(kmer) == 0
    assert cg.count(kmer) == 1
    assert cg.get(kmer) == 1
