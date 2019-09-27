import os
import re

from converter.datamodel.submittable import AccessionedSubmittable
from converter.datamodel.components import Attribute


class Study(AccessionedSubmittable):
    """
    Attributes of the study
    :param alias: string, unique study name (auto-generated from MAGE-TAB file name)
    :param accession: string, ArrayExpress experiment accession
    :param title: string, (optional)
    :param description: string, free-text description of study background and aim
    :param protocolrefs: list,  protocol accessions/name used in the study
    :param projectref: string, project alias or accession
    :param experimental_factor: list, Attribute class objects
    :param experimental_design: list, Attribute class objects
    :param experiment_type: list, experiment type ontology terms
    :param date_of_experiment: string, (optional) date of experiment
    :param study_type: string, submission type (microarray, sequencing, singlecell)
    :param secondary_accession: list, accession of the same study in another archive (e.g. ENA)
    :param related_experiment: list, accession of experiments related to this study
    :param comments: dictionary, (optional) known IDF comments
                    keys: string, comment category
                    values: list, comment values
    """
    def __init__(self, **kwargs):
        AccessionedSubmittable.__init__(self, **kwargs)
        self.title = kwargs.get("title")
        self.protocolrefs = kwargs.get("protocolrefs", [])
        self.projectref = kwargs.get("projectref")
        self.experimental_factor = kwargs.get("experimental_factor", [])
        self.experimental_design = kwargs.get("experimental_design", [])
        self.experiment_type = kwargs.get("experiment_type", [])
        self.date_of_experiment = kwargs.get("date_of_experiment")
        self.submission_type = kwargs.get("submission_type")
        self.secondary_accession = kwargs.get("secondary_accession", [])
        self.related_experiment = kwargs.get("related_experiment", [])
        self.comments = kwargs.get("comments", {})

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, accession={self.accession}, title={self.title}, " \
               "description={self.description}, protocolrefs={self.protocolrefs}, projectref={self.projectref}, " \
               "experimental_factor={self.experimental_factor}, experimental_design={self.experimental_design}, " \
               "experiment_type={self.experiment_type}, date_of_experiment={self.date_of_experiment}, " \
               "submission_type={self.submission_type}, secondary_accession={self.secondary_accession}, " \
               "related_experiment={self.related_experiment}, comments={self.comments})".format(self=self)

    @classmethod
    def from_magetab(cls, study_info):
        accession = study_info.get("accession")
        idf_file = os.path.basename(study_info.get("idf_filename", ""))
        alias = re.sub(r"\.idf\.txt$", "", idf_file)
        projectref = "project_" + alias
        title = study_info.get("title")
        description = study_info.get("description")
        ef = study_info.get("experimental_factor", [])
        ef_objects = [Attribute(value=d.get("experimental_factor"),
                                unit=None,
                                term_accession=d.get("term_accession"),
                                term_source=d.get("term_source")) for d in ef]
        ed = study_info.get("experimental_design", [])
        ed_objects = [Attribute(value=d.get("experimental_design"),
                                unit=None,
                                term_accession=d.get("term_accession"),
                                term_source=d.get("term_source")) for d in ed]
        protocolrefs = study_info.get("protocolRefs", [])
        date_of_experiment = study_info.get("date_of_experiment", None)
        comments = study_info.get("comments", {})
        experiment_type = comments.get("experiment_type", [])

        return cls(alias=alias,
                   accession=accession,
                   title=title,
                   description=description,
                   protocolrefs=protocolrefs,
                   projectref=projectref,
                   experimental_factor=ef_objects,
                   experimental_design=ed_objects,
                   experiment_type=experiment_type,
                   date_of_experiment=date_of_experiment)
