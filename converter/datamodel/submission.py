"""The submission class"""

from collections import OrderedDict
from typing import List, Union

from utils.common_utils import get_ontology_source_file
from converter.datamodel.project import Project
from converter.datamodel.study import Study
from converter.datamodel.protocol import Protocol
from converter.datamodel.sample import Sample
from converter.datamodel.assay import SeqAssay, SingleCellAssay, MicroarrayAssay
from converter.datamodel.data import AssayData, Analysis


class Submission:
    """
    The submission class holds all objects that make up an ArrayExpress experiment.

    :param info: dictionary, information about the submission (file names, submitter details, etc.)
    :param project: object, Project class object
    :param study: object, Study class object
    :param protocol: list, Protocol class objects
    :param sample: list, Sample class objects
    :param assay: list, SeqAssay or MicroarrayAssay class objects
    :param assay_data: list, AssayData class objects
    :param analysis: list, Analysis class objects
    """
    info: dict
    project: Project
    study: Study
    protocol: List[Protocol]
    sample: List[Sample]
    assay: List[Union[SeqAssay, SingleCellAssay, MicroarrayAssay]]
    assay_data: List[AssayData]
    analysis: List[Analysis]

    def __init__(self, info, project, study, protocol, sample, assay, assay_data, analysis):
        self.info = info
        self.project = project
        self.study = study
        self.protocol = protocol
        self.sample = sample
        self.assay = assay
        self.assay_data = assay_data
        self.analysis = analysis

        # Create indexes so we can return given objects based on alias
        self.sample_lookup = {s.alias: s for s in self.sample}
        self.assay_lookup = {a.alias: a for a in self.assay}
        self.protocol_lookup = {p.alias: p for p in self.protocol}
        self.data_lookup = {ad.alias: ad for ad in self.assay_data}

    def get_sample(self, sample_alias):
        """Return the sample object for a given sample alias."""
        return self.sample_lookup.get(sample_alias)

    def get_assay(self, assay_alias):
        """Return the assay object for a given assay alias."""
        return self.assay_lookup.get(assay_alias)

    def get_protocol(self, protocol_alias):
        """Return the protocol object for a given protocol alias."""
        return self.protocol_lookup.get(protocol_alias)

    def get_assay_data(self, ad_alias):
        """Return the assay data object for a given assay data alias."""
        return self.data_lookup.get(ad_alias)

    def get_term_sources(self):
        term_sources = OrderedDict()
        # Make sure we have at least EFO (used for protocol types etc.)
        term_sources["EFO"] = "https://www.ebi.ac.uk/efo.owl"
        ontologies = {a.term_source for s in self.sample for a in s.attributes.values() if a.term_source}
        for o in ontologies:
            if not o.upper() == "EFO":
                term_sources[o] = get_ontology_source_file(o)
        # MA submissions need ArrayExpress as Term Ref for array design accessions
        if self.info.get("submission_type") == "microarray":
            term_sources["ArrayExpress"] = "https://www.ebi.ac.uk/arrayexpress/"
        return term_sources
