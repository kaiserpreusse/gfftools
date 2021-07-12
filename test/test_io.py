import pytest

from gfftools.io import GffReader
from gfftools.model import Record, Sequence, SequenceRegion
import os


@pytest.fixture
def ensembl_gff_file_path(files_directory):
    return os.path.join(files_directory, 'Homo_sapiens.GRCh38.103.chromosome.10.gff3.gz')


@pytest.fixture
def gff_with_fasta_path(files_directory):
    return os.path.join(files_directory, 'gff_with_fasta.gff3')


class TestGffReaderStreamMode:

    def test_header(self, ensembl_gff_file_path):

        gff_reader = GffReader(ensembl_gff_file_path)

        parsed_header = gff_reader.header()

        expected_header = {'gff-version': '3', 'sequence-region': '10 1 133797422', 'genome-build': 'GRCh38.p13',
                           'genome-version': 'GRCh38', 'genome-date': '2013-12',
                           'genome-build-accession': 'NCBI:GCA_000001405.28', 'genebuild-last-updated': '2020-08'}

        for k, v in expected_header.items():
            assert parsed_header[k] == v

    def test_record_iterator(self, ensembl_gff_file_path):

        gff_reader = GffReader(ensembl_gff_file_path)

        # assert that we have records
        assert list(gff_reader.records)

        for r in gff_reader.records:
            assert isinstance(r, Record)

    def test_sequence_region_iterator(self, ensembl_gff_file_path):

        gff_reader = GffReader(ensembl_gff_file_path)

        # assert that we have sequence regions
        assert list(gff_reader.sequence_regions)

        for r in gff_reader.sequence_regions:
            assert isinstance(r, SequenceRegion)
            assert isinstance(r.start, int)
            assert isinstance(r.end, int)

    def test_sequence_iterator(self, gff_with_fasta_path):

        gff_reader = GffReader(gff_with_fasta_path)

        list_of_sequences = list(gff_reader.sequences)
        assert list_of_sequences
        assert len(list_of_sequences) == 2

        for r in gff_reader.sequences:
            assert isinstance(r, Sequence)
            assert isinstance(r.id, str)
            assert isinstance(r.sequence, str)
            assert ' ' not in r.sequence
            assert '\n' not in r.sequence
            assert '\t' not in r.sequence


class TestGffReaderParseMode:

    def test_parser_run(self, gff_with_fasta_path):
        gff_reader = GffReader(gff_with_fasta_path)

        gff_reader.parse()

        assert isinstance(gff_reader.records, list)
        assert isinstance(gff_reader.sequences, list)
        assert isinstance(gff_reader.sequence_regions, list)

        assert len(gff_reader.records) > 0
        assert len(gff_reader.sequences) > 0
        assert len(gff_reader.sequence_regions) > 0

        for r in gff_reader.records:
            assert isinstance(r, Record)

        for r in gff_reader.sequence_regions:
            assert isinstance(r, SequenceRegion)
            assert isinstance(r.start, int)
            assert isinstance(r.end, int)

        for r in gff_reader.sequences:
            assert isinstance(r, Sequence)
            assert isinstance(r.id, str)
            assert isinstance(r.sequence, str)
            assert ' ' not in r.sequence
            assert '\n' not in r.sequence
            assert '\t' not in r.sequence
