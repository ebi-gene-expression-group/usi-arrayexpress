"""The Project class"""

from converter.datamodel.submittable import AccessionedSubmittable
from converter.datamodel.components import Contact, Publication


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
