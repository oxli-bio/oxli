import pytest

import oxli
from test_basic import create_sample_kmer_table


# Set operations
def test_union():
    table1 = create_sample_kmer_table(3, ["AAA", "AAC"])
    table2 = create_sample_kmer_table(3, ["AAC", "AAG"])

    union_set = table1.union(table2)
    expected_union = set(table1.hashes).union(table2.hashes)

    assert union_set == expected_union, "Union of hash sets should match"


def test_intersection():
    table1 = create_sample_kmer_table(3, ["AAA", "AAC"])
    table2 = create_sample_kmer_table(3, ["AAC", "AAG"])

    intersection_set = table1.intersection(table2)
    expected_intersection = set(table1.hashes).intersection(table2.hashes)

    assert intersection_set == expected_intersection, (
        "Intersection of hash sets should match"
    )


def test_difference():
    table1 = create_sample_kmer_table(3, ["AAA", "AAC"])
    table2 = create_sample_kmer_table(3, ["AAC", "AAG"])

    difference_set = table1.difference(table2)
    expected_difference = set(table1.hashes).difference(table2.hashes)

    assert difference_set == expected_difference, "Difference of hash sets should match"


def test_symmetric_difference():
    table1 = create_sample_kmer_table(3, ["AAA", "AAC"])
    table2 = create_sample_kmer_table(3, ["AAC", "AAG"])

    symmetric_difference_set = table1.symmetric_difference(table2)
    expected_symmetric_difference = set(table1.hashes).symmetric_difference(
        table2.hashes
    )

    assert symmetric_difference_set == expected_symmetric_difference, (
        "Symmetric difference of hash sets should match"
    )


def test_dunder_methods():
    table1 = create_sample_kmer_table(3, ["AAA", "AAC"])
    table2 = create_sample_kmer_table(3, ["AAC", "AAG"])

    assert table1.__or__(table2) == table1.union(table2), (
        "__or__ method should match union()"
    )
    assert table1.__and__(table2) == table1.intersection(table2), (
        "__and__ method should match intersection()"
    )
    assert table1.__sub__(table2) == table1.difference(table2), (
        "__sub__ method should match difference()"
    )
    assert table1.__xor__(table2) == table1.symmetric_difference(table2), (
        "__xor__ method should match symmetric_difference()"
    )
