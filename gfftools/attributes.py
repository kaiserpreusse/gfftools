def read_gff_attributes(attribute_column: str) -> dict:
    """
    Parse attributes for a GFF3 record. Attributes with pre-defined meaning are parsed according to their
    specification (e.g. Dbxref usually has multiple values which are split up: 'GeneID:1234,Genbank:NM_9283').


    :param attribute_column: Attribute column of a GFF3 file.
    :return: Dictionary of attributes.
    :rtype: dict
    """
    attributes = {}
    for a in attribute_column.strip().split(';'):
        # there is a leading space for some fields
        a = a.strip()
        # an attribute looks like 'Alias=MIMAT0027693'
        # some files contain additional whitespaces or empty elements
        # only continue if a has content
        if a:
            key = a.split('=')[0]
            value = a.split('=')[1]

            # handle pre-defined attributes
            if key == 'Dbxref':
                # a dbxref line looks like: GeneID:1234,Genbank:NM_9283
                dbxref_dict = {}
                for dbxref_entry in value.split(','):
                    dbxref_key, dbxref_value = dbxref_entry.split(':', 1)
                    dbxref_dict[dbxref_key] = dbxref_value
                value = dbxref_dict

            attributes[key] = value

    return attributes


def read_gtf_attributes(attribute_column: str):
    attributes = {}
    for a in attribute_column.split(';'):
        # there is a leading space for some fields
        a = a.strip()
        try:
            # an attribute looks like 'gene_id "ENSG00000274890"'
            key = a.split()[0].replace('"', '')
            value = a.split()[1].replace('"', '')
            attributes[key] = value
        except IndexError:
            pass

    return attributes


def get_attributes(file_type, flds):
    # get attribtues and parse depending on type
    # gff files contain "attribute=some_value; another=some_value"
    # gtf files contain "gene_id "ENSG00000223972"; gene_name "DDX11L1";"
    attributes = None
    if file_type == 'gff':
        attributes = read_gff_attributes(flds[8])
    elif file_type == 'gtf':
        attributes = read_gtf_attributes(flds[8])
    return attributes
