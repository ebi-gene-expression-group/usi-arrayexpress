
class DataFile:
    def __init__(self, **kwargs):
        self.name = kwargs.get("name")
        self.checksum = kwargs.get("checksum")
        self.checksum_method = kwargs.get("checksum_method")
        self.ftp_location = kwargs.get("ftp_location")
        self.read = kwargs.get("read")

    def __repr__(self):
        return "{self.__class__.__name__}(name={self.name}, checksum={self.checksum}, " \
               "checksum_method={self.checksum_method}, ftp_location={self.ftp_location}, " \
               "read={self.read})".format(self=self)


class Contact:
    def __init__(self, **kwargs):
        self.firstName = kwargs.get("firstName")
        self.lastName = kwargs.get("lastName")
        self.email = kwargs.get("email")
        self.affiliation = kwargs.get("affiliation")
        self.address = kwargs.get("address")
        self.phone = kwargs.get("phone")
        self.roles = kwargs.get("roles", [])
        self.middleInitials = kwargs.get("middleInitials")
        self.fax = kwargs.get("fax")
        self.orcidId = kwargs.get("orcidId")

    def __repr__(self):
        return "{self.__class__.__name__}(firstName={self.firstName}, lastName={self.lastName}, email{self.email}, " \
               "affiliation={self.affiliation}, address={self.address}, phone={self.phone}, roles={self.roles}, " \
               "middleInitials={self.middleInitials}, fax={self.fax}, orcidId={self.orcidId})".format(self=self)


class Publication:
    def __init__(self, **kwargs):
        self.articleTitle = kwargs.get("articleTitle")
        self.authors = kwargs.get("authors")
        self.pubmedId = kwargs.get("pubmedId")
        self.doi = kwargs.get("doi")
        self.publicationStatus = kwargs.get("publicationStatus")

    def __repr__(self):
        return "{self.__class__.__name__}(articleTitle={self.articleTitle}, authors={self.authors}, " \
               "pubmedId={self.pubmedId}, doi={self.doi}, status={self.publicationStatus})".format(self=self)


class Attribute:
    def __init__(self, **kwargs):
        self.value = kwargs.get("value")
        self.unit = kwargs.get("unit")
        self.term_source = kwargs.get("term_source")
        self.term_accession = kwargs.get("term_accession")

    def __repr__(self):
        return "{self.__class__.__name__}(value={self.value}, unit={self.unit}, " \
               "term_accession={self.term_accession}, term_source={self.term_source})".format(self=self)


class Unit:
    def __init__(self, **kwargs):
        self.value = kwargs.get("value")
        self.unit_type = kwargs.get("unit_type")
        self.term_source = kwargs.get("term_source")
        self.term_accession = kwargs.get("term_accession")

    def __repr__(self):
        return "{self.__class__.__name__}(value={self.value}, unit_type={self.unit_type}, " \
               "term_source={self.term_source}, term_accession={self.term_accession})".format(self=self)

