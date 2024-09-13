import oxli
import pytest
from test_basic import create_sample_kmer_table


def test_len_dunder_method():
    '''__len__ should return number of keys in KmerCountTable.'''
    pass

def test_iter_dunder_method():
    '''KmerCountTable should be iterable, yield hash:count pairs'''
    pass

def test_next_dunder_method():
    '''Select next key in generator'''
    pass

def test_getitem_dunder_method():
    '''Query an object to using the indexing syntax (obj[key])'''
    # Same behaviour as .get()
    pass

def test_setitem_dunder_method():
    '''Set values using the indexing syntax (obj[key] = value)'''
    pass