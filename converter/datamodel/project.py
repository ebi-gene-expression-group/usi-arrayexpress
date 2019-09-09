import re
import os

from converter.datamodel.submittable import AccessionedSubmittable
from converter.datamodel.components import Contact, Publication
from utils.converter_utils import get_controlled_vocabulary


class Project(AccessionedSubmittable):
    """
       Attributes of the project
       :param alias: string, unique project name (auto-generated from MAGE-TAB file name)
       :param accession: string, BioStudies accession
       :param title: (optional) string, submitter provided project description
       :param description: string, (optional) free-text project description
       :param releaseDate: string, date of public release
       :param publications: list, Publication class objects
       :param contacts: list, Contact class objects
       """
    def __init__(self, **kwargs):
        AccessionedSubmittable.__init__(self, **kwargs)
        self.title = kwargs.get("title")
        self.releaseDate = kwargs.get("releaseDate")
        self.publications = [Publication(**p) for p in kwargs.get("publications", [])]
        self.contacts = [Contact(**c) for c in kwargs.get("contacts", [])]

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, accession={self.accession}, " \
               "releaseDate={self.releaseDate}, publications={self.publications}, " \
               "contacts={self.contacts})".format(self=self)

    @classmethod
    def from_magetab(cls, study_info):
        idf_file = os.path.basename(study_info.get("idf_filename", ""))
        temp_project_name = re.sub("\.idf\.txt$", "", idf_file)
        alias = "project_" + temp_project_name
        # The project accession might be the BioStudies accession, which could be saved as IDF comment
        # Placeholder code which returns None for now
        comments = study_info.get("comments")
        accession = comments.get("biostudiesaccession", None)
        releaseDate = study_info.get("releaseDate", None)
        contact_terms = get_controlled_vocabulary("contact_terms")
        contacts_raw = study_info.get("contacts", [])
        contacts = [{contact_terms["Person First Name"]: c.get(contact_terms["Person First Name"]),
                     contact_terms["Person Last Name"]: c.get(contact_terms["Person Last Name"]),
                     contact_terms["Person Email"]: c.get(contact_terms["Person Email"]),
                     contact_terms["Person Affiliation"]: c.get(contact_terms["Person Affiliation"]),
                     contact_terms["Person Address"]: c.get(contact_terms["Person Address"]),
                     contact_terms["Person Phone"]: c.get(contact_terms["Person Phone"]),
                     contact_terms["Person Roles"]: re.split("\s*;\s*", c.get(contact_terms["Person Roles"], "")),
                     contact_terms["Person Mid Initials"]: c.get(contact_terms["Person Mid Initials"]),
                     contact_terms["Person Fax"]: c.get(contact_terms["Person Fax"])}
                    for c in contacts_raw]

        # Get publication terms and create list of
        publications_raw = study_info.get("publications", [])
        pub_terms = get_controlled_vocabulary("publication_terms")
        publications = [{pub_terms["Publication Title"]: pub.get(pub_terms["Publication Title"]),
                        pub_terms["Publication Author List"]: pub.get(pub_terms["Publication Author List"]),
                        pub_terms["PubMed ID"]: pub.get(pub_terms["PubMed ID"]),
                        pub_terms["Publication DOI"]: pub.get(pub_terms["Publication DOI"]),
                        pub_terms["Publication Status"]: pub.get(pub_terms["Publication Status"])}
                        for pub in publications_raw]
        title = study_info.get("title")
        description = study_info.get("description")

        return cls(alias=alias,
                   accession=accession,
                   title=title,
                   description=description,
                   releaseDate=releaseDate,
                   publications=publications,
                   contacts=contacts)
