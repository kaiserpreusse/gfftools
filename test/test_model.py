import pytest

from gfftools.model import Record

@pytest.fixture(scope='module')
def record_object():
    # ctg123	.	CDS	5000	5500	.	+	0	ID=cds00001;Parent=mRNA00001
    r = Record(seqid='ctg123', source='test', type='CDS', start=1, end=500, score='.', strand='+', phase='.',
               attributes={'ID': 'gene00001', 'Name': 'EDEN'})
    return r


class TestRecord:

    def test_to_dict(self, record_object):
        d = record_object.to_dict()

        expected_d = {'ID': 'gene00001', 'Name': 'EDEN', 'seqid': 'ctg123', 'source': 'test', 'type': 'CDS', 'start': 1, 'end': 500, 'score': '.', 'strand': '+', 'phase': '.'}

        assert isinstance(d, dict)

        for k,v in expected_d.items():
            assert d[k] == v
