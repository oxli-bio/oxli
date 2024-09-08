import pytest
import oxli

def create_sample_kmer_table(ksize, kmers):
    table = oxli.KmerCountTable(ksize)
    for kmer in kmers:
        table.count(kmer)
    return table


def test_simple():
    # yo dawg it works
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    assert cg.get(kmer) == 0
    assert cg.count(kmer) == 1
    assert cg.get(kmer) == 1


def test_wrong_ksize():
    # but only with the right ksize
    cg = oxli.KmerCountTable(3)
    kmer = "ATCG"

    with pytest.raises(ValueError):
        cg.count(kmer)

    with pytest.raises(ValueError):
        cg.get(kmer)


def test_consume():
    # test basic consume
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    assert cg.consume(kmer) == 1
    assert cg.get("ATCG") == 1


def test_consume_2():
    # test reverse complement
    cg = oxli.KmerCountTable(4)
    seq = "ATCGG"

    assert cg.consume(seq) == 2
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1 # reverse complement!


def test_consume_bad_DNA():
    # test an invalid base in last position
    cg = oxli.KmerCountTable(4)
    seq = "ATCGGX"
    with pytest.raises(ValueError,
                        match="bad k-mer encountered at position 2"):
        cg.consume(seq, allow_bad_kmers=False)


def test_consume_bad_DNA_2():
    # test an invalid base in first position
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    with pytest.raises(ValueError,
                        match="bad k-mer encountered at position 0"):
        cg.consume(seq, allow_bad_kmers=False)


def test_consume_bad_DNA_ignore():
    # we can ignore bad DNA
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    print(cg.consume(seq, allow_bad_kmers=True))
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1 # rc


def test_consume_bad_DNA_ignore_is_default():
    # ignoring bad DNA is default
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    print(cg.consume(seq))
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1 # rc


def test_count_get():
    # test a bug reported by adam taranto: count and get should work together!
    kmer = 'TAAACCCTAACCCTAACCCTAACCCTAACCC'

    cg = oxli.KmerCountTable(ksize=31)
    hashkey = cg.hash_kmer(kmer)

    assert cg.get(kmer) == 0
    assert cg.count(kmer) == 1
    assert cg.count(kmer) == 2

    x = cg.get(kmer)
    assert x == 2, x
    

def test_get_hash():
    ''' Retrieve counts using hash key. '''
    table = create_sample_kmer_table(3, ['AAA', 'TTT', 'AAC'])
    # Find hash of kmer 'AAA'
    hash_aaa = table.hash_kmer('AAA') # 10679328328772601858
    # Lookup counts for hash of 'AAA' and rc 'TTT'
    count_aaa = table.get_hash(hash_aaa)
    
    assert count_aaa == 2, "Hash count for 'AAA' should be 2"

    hash_aac = table.hash_kmer('AAC') # 6579496673972597301
    count_aac = table.get_hash(hash_aac)
    
    assert count_aac == 1, "Hash count for 'AAC' should be 1"

def test_get_hash_array():
    ''' 
    Get vector of counts corresponding to vector of hash keys.
    '''
    table = create_sample_kmer_table(3, ['AAA', 'AAC'])
    hash_aaa = table.hash_kmer('AAA')
    hash_aac = table.hash_kmer('AAC')
    hash_ggg = table.hash_kmer('GGG') # key not in table
    
    hash_keys = [hash_aaa, hash_aac, hash_ggg]
    counts = table.get_hash_array(hash_keys)
    
    assert counts == [1, 1, 0], "Hash array counts should match the counts of 'AAA' and 'AAC' and return zero for 'GGG'."

def test_union():
    table1 = create_sample_kmer_table(3, ['AAA', 'AAC'])
    table2 = create_sample_kmer_table(3, ['AAC', 'AAG'])
    
    union_set = table1.union(table2)
    expected_union = set(table1.hashes).union(table2.hashes)
    
    assert union_set == expected_union, "Union of hash sets should match"

def test_intersection():
    table1 = create_sample_kmer_table(3, ['AAA', 'AAC'])
    table2 = create_sample_kmer_table(3, ['AAC', 'AAG'])
    
    intersection_set = table1.intersection(table2)
    expected_intersection = set(table1.hashes).intersection(table2.hashes)
    
    assert intersection_set == expected_intersection, "Intersection of hash sets should match"

def test_difference():
    table1 = create_sample_kmer_table(3, ['AAA', 'AAC'])
    table2 = create_sample_kmer_table(3, ['AAC', 'AAG'])
    
    difference_set = table1.difference(table2)
    expected_difference = set(table1.hashes).difference(table2.hashes)
    
    assert difference_set == expected_difference, "Difference of hash sets should match"

def test_symmetric_difference():
    table1 = create_sample_kmer_table(3, ['AAA', 'AAC'])
    table2 = create_sample_kmer_table(3, ['AAC', 'AAG'])
    
    symmetric_difference_set = table1.symmetric_difference(table2)
    expected_symmetric_difference = set(table1.hashes).symmetric_difference(table2.hashes)
    
    assert symmetric_difference_set == expected_symmetric_difference, "Symmetric difference of hash sets should match"

def test_dunder_methods():
    table1 = create_sample_kmer_table(3, ['AAA', 'AAC'])
    table2 = create_sample_kmer_table(3, ['AAC', 'AAG'])
    
    assert table1.__or__(table2) == table1.union(table2), "__or__ method should match union()"
    assert table1.__and__(table2) == table1.intersection(table2), "__and__ method should match intersection()"
    assert table1.__sub__(table2) == table1.difference(table2), "__sub__ method should match difference()"
    assert table1.__xor__(table2) == table1.symmetric_difference(table2), "__xor__ method should match symmetric_difference()"

def test_hashes_attribute():
    table = create_sample_kmer_table(3, ['AAA', 'TTT', 'AAC'])
    hashes = table.hashes
    hash_aaa = table.hash_kmer('AAA') #10679328328772601858
    hash_ttt = table.hash_kmer('TTT') #10679328328772601858
    hash_aac = table.hash_kmer('AAC') #6579496673972597301
    
    expected_hashes = set([hash_aaa, hash_ttt, hash_aac]) #{10679328328772601858, 6579496673972597301}
    assert set(hashes) == expected_hashes, ".hashes attribute should match the expected set of hash keys"