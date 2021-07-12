from dataclasses import dataclass
from enum import Enum

DEFINED_ATTRIBUTES = ('ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Derives_from', 'Note',
                      'Dbxref', 'Ontology_term', 'Is_circular')


@dataclass
class Record:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: dict

    def get_attribute(self, key):
        if key in self.attributes:
            return self.attributes[key]

    @property
    def ID(self):
        return self.get_attribute('ID')

    @property
    def Name(self):
        return self.get_attribute('Name')

    @property
    def Alias(self):
        return self.get_attribute('Alias')

    @property
    def Parent(self):
        return self.get_attribute('Parent')

    @property
    def Target(self):
        return self.get_attribute('Target')

    @property
    def Gap(self):
        return self.get_attribute('Gap')

    @property
    def Derives_from(self):
        return self.get_attribute('Derives_from')

    @property
    def Note(self):
        return self.get_attribute('Note')

    @property
    def Dbxref(self):
        return self.get_attribute('Target')

    @property
    def Ontology_term(self):
        return self.get_attribute('Ontology_term')

    @property
    def Is_circular(self):
        return self.get_attribute('Is_circular')


@dataclass
class SequenceRegion:
    seqid: str
    start: int
    end: int


@dataclass
class Sequence:
    id: str
    sequence: str = ''
