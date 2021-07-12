import gzip

from gfftools.attributes import read_gff_attributes, read_gtf_attributes
from gfftools.model import Record, SequenceRegion, Sequence


class GffReader(object):
    """
    Read GFF3 and GTF files and return a set of Records.
    """

    def __init__(self, gff_file_path, file_type=None):
        self.gff_file_path = gff_file_path

        # determine file type (gtf/gff) from file name if not supplied
        if not file_type:
            if '.gtf' in gff_file_path:
                self.file_type = 'gtf'
            elif '.gff' in gff_file_path:
                self.file_type = 'gff'
            else:
                raise IOError("Cannot guess file type from filename: %s" % gff_file_path)
        else:
            if file_type in ['gff', 'gtf']:
                self.file_type = file_type
            else:
                raise AttributeError("file_type must be 'gtf' or 'gff'")

        # other file feature flags
        self.has_fasta = False
        self.has_sequence_regions = False

    @property
    def records(self):
        return self.iterate_records()

    @property
    def sequence_regions(self):
        return self.iterate_sequence_regions()

    @property
    def sequences(self):
        return self.iterate_sequences()

    def header(self):
        """
        Most GFF files have additional information in the header. Some is encoded with pragmas (such as ##gff-version)
        other is just text or non-standard encodings such as `#!genome-build` in Gencode files.

        This functions tries to collect in a dictionary.
        """
        print(self.gff_file_path)
        if self.gff_file_path.endswith('.gz'):
            f = gzip.open(self.gff_file_path, 'rt')
        else:
            f = open(self.gff_file_path, 'rt')

        header_lines = []

        for line in f:
            # break on first record
            if not line.startswith('#'):
                break
            header_lines.append(line.strip())

        header_dict = {}
        comment_text = ''
        for line in header_lines:
            if line.startswith('##') or line.startswith('#!'):
                # remove comment
                line = line.replace('##', '').replace('#!', '')
                print(line)
                # split on first whitespace
                key, value = line.split(None, 1)
                # clean leftover whitespaces
                key = key.strip()
                value = value.strip()
                header_dict[key] = value

            if line.startswith('#') and not line.startswith('#!'):
                comment_text += line.strip()

        f.close()
        if comment_text:
            header_dict['comment'] = comment_text
        return header_dict

    def analyze_file(self):
        """
        Analyze the file and determine if it has sequence regions, a FASTA part or other elements that
        influence logic of the parser.
        """
        if self.gff_file_path.endswith('.gz'):
            f = gzip.open(self.gff_file_path, 'rt')
        else:
            f = open(self.gff_file_path, 'rt')

        for l in f:
            if l.startswith('##sequence-region'):
                self.has_sequence_regions = True
            if l.startswith('##FASTA'):
                self.has_fasta = True

        f.close()

    def iterate_records(self) -> Record:
        """
        :return: Record namedtuple
        """

        if self.gff_file_path.endswith('.gz'):
            f = gzip.open(self.gff_file_path, 'rt')
        else:
            f = open(self.gff_file_path, 'rt')

        for l in f:

            # break when reaching ##FASTA
            if l.lower().startswith('##fasta'):
                break

            if not l.startswith('#'):
                flds = l.split('\t')

                # get attribtues and parse depending on type
                # gff files contain "attribute=some_value; another=some_value"
                # gtf files contain "gene_id "ENSG00000223972"; gene_name "DDX11L1";"
                attributes = None
                if self.file_type == 'gff':
                    attributes = read_gff_attributes(flds[8])
                elif self.file_type == 'gtf':
                    attributes = read_gtf_attributes(flds[8])

                yield Record(flds[0], flds[1], flds[2], int(flds[3]), int(flds[4]), flds[5], flds[6], flds[7],
                             attributes)

        f.close()

    def iterate_sequence_regions(self) -> Record:

        if self.gff_file_path.endswith('.gz'):
            f = gzip.open(self.gff_file_path, 'rt')
        else:
            f = open(self.gff_file_path, 'rt')

        for line in f:

            # break when reaching ##FASTA
            if line.lower().startswith('##fasta'):
                break

            if line.startswith('##sequence-region'):
                flds = line.split()

                yield SequenceRegion(flds[1], int(flds[2]), int(flds[3]))

        f.close()

    def iterate_sequences(self):
        """
        Iterate until first occurence of ##FASTA and start parsing sequence records.
        """
        if self.gff_file_path.endswith('.gz'):
            f = gzip.open(self.gff_file_path, 'rt')
        else:
            f = open(self.gff_file_path, 'rt')

        parse_records = False
        current_sequence_record = None
        for line in f:
            # wait until ##FASTA was found
            if parse_records:

                # check if new Sequence begins
                if line.startswith('>'):
                    # if we have a Sequence already, yield it
                    if current_sequence_record:
                        yield current_sequence_record
                    # create a new Sequence
                    record_id = line.strip().replace('>', '')
                    current_sequence_record = Sequence(id=record_id)
                else:
                    # if not new record, just append to the sequence
                    current_sequence_record.sequence += line.strip()

            # start parsing when reaching ##FASTA
            if line.lower().startswith('##fasta'):
                parse_records = True

        # yield last Sequence if it exists
        if current_sequence_record:
            yield current_sequence_record

        f.close()
