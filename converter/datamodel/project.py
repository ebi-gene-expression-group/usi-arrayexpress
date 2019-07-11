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
        self.publications = kwargs.get("publications", [])
        self.contacts = kwargs.get("contacts", [])

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
        contacts = [Contact(c.get(contact_terms["Person First Name"]),
                            c.get(contact_terms["Person Last Name"]),
                            c.get(contact_terms["Person Email"]),
                            c.get(contact_terms["Person Affiliation"]),
                            c.get(contact_terms["Person Address"]),
                            c.get(contact_terms["Person Phone"]),
                            c.get(contact_terms["Person Roles"]),
                            c.get(contact_terms["Person Mid Initials"]),
                            c.get(contact_terms["Person Fax"])
                            ) for c in contacts_raw]
        # Transform roles to list
        for c in contacts:
            if c.roles:
                roles = c.roles
                c.roles = roles.split(';')
        # Get publication terms and create list of
        publications_raw = study_info.get("publications", [])
        pub_terms = get_controlled_vocabulary("publication_terms")
        publications = [Publication(pub.get(pub_terms["Publication Title"]),
                                    pub.get(pub_terms["Publication Author List"]),
                                    pub.get(pub_terms["PubMed ID"]),
                                    pub.get(pub_terms["Publication DOI"]),
                                    pub.get(pub_terms["Publication Status"])) for pub in publications_raw]

        title = study_info.get("title")
        description = study_info.get("description")

        return cls(alias=alias,
                   accession=accession,
                   title=title,
                   description=description,
                   releaseDate=releaseDate,
                   publications=publications,
                   contacts=contacts)

    @classmethod
    def from_dict(cls, project_dict):
        return cls(alias=project_dict.get("alias"),
                   accession=project_dict.get("accession"),
                   title=project_dict.get("title"),
                   description=project_dict.get("description"),
                   releaseDate=project_dict.get("releaseDate"),
                   publications=[Publication.from_dict(p) for p in project_dict.get("publications", [])],
                   contacts=[Contact.from_dict(p) for p in project_dict.get("contacts", [])])
