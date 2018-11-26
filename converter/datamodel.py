from parsing import read_json_file, get_controlled_vocabulary, remove_duplicates
from converting import is_accession
import os
import re


class Sample:

    def __init__(self, alias, accession, taxon, taxonId, characteristics, factors, material_type, description):
        pass


class Assay:

    def __init__(self, alias, accession, technology_type, protocolrefs, sampleref):
        # Attributes common to all assay types
        self.alias = alias
        self.accession = accession
        self.technology_type = technology_type
        self.protocolrefs = protocolrefs
        self.sampleref = sampleref

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.technology_type}, " \
               "{self.protocolrefs}, {self.sampleref})".format(self=self)

    def get_attributes(self):
        """This returns a list of all additional attributes that have values"""
        all_attributes = self.__dict__.keys()
        all_values = self.__dict__
        fixed_attributes = ("alias", "accession", "protocolrefs", "sampleref")
        other_attributes = list()
        for a in all_attributes:
            if a not in fixed_attributes and all_values[a]:
                other_attributes.append(a)
        return other_attributes


class MicroarrayAssay(Assay):
    def __init__(self, alias, accession, technology_type, protocolrefs, sampleref, label, arrayref):
        Assay.__init__(self, alias, accession, technology_type, protocolrefs, sampleref)
        self.label = label
        self.arrayref = arrayref

    @classmethod
    def from_magetab(cls, le_attributes, extract_attributes, assay_attributes):
        """Intialise assay attributes from MAGE-TAB data dicts.
        The central node for the microarray assay is the Labeled Extract Name."""

        alias = le_attributes.get("name")
        accession = None
        technology_type = remove_duplicates([a.get("technology_type") for a in assay_attributes])

        # Get all protocol refs
        protocolrefs = []
        for a in assay_attributes:
            protocolrefs.extend(a.get("protocol_ref"))
        protocolrefs.extend(extract_attributes.get("protocol_ref", []))
        protocolrefs.extend(le_attributes.get("protocol_ref", []))
        protocolrefs = remove_duplicates(protocolrefs)

        sampleref = extract_attributes.get("sample_ref")

        label = le_attributes.get("label")

        # Get Array design REF, we are only expecting one unique per extract
        arrayref = [a.get("array_design") for a in assay_attributes]

        return cls(alias, accession, technology_type[0], protocolrefs, sampleref, label, arrayref[0])


class SeqAssay(Assay):
    def __init__(self, alias, accession, technology_type, protocolrefs, sampleref, lib_attribs):
        Assay.__init__(self, alias, accession, technology_type, protocolrefs, sampleref)

        self.library_layout = lib_attribs.get("library_layout")
        self.library_selection = lib_attribs.get("library_selection")
        self.library_strategy = lib_attribs.get("library_strategy")
        self.library_strand = lib_attribs.get("library_strand")
        self.library_source = lib_attribs.get("library_source")
        self.orientation = lib_attribs.get("orientation")
        self.nominal_length = lib_attribs.get("nominal_length")
        self.nominal_sdev = lib_attribs.get("nominal_sdev")
        self.platform_type = lib_attribs.get("platform_type")
        self.instrument_model = lib_attribs.get("instrument_model")

    @classmethod
    def from_magetab(cls, extract_attributes, assay_attributes, protocols):
        """Intialise assay attributes from MAGE-TAB data dicts.
        The central node for the sequencing assay is the Extract Name."""

        alias = extract_attributes.get("name")

        # Get library attributes from extract comments
        comments = extract_attributes.get("comments")
        lib_attrib_cv = get_controlled_vocabulary("sdrf_comments_ena").keys()
        lib_attribs = {a.lower(): comments[a] for a in lib_attrib_cv if comments.get(a)}

        # Get technology type(s) from assay attributes
        technology_type = remove_duplicates([a.get("technology_type") for a in assay_attributes])

        # Get accession from ENA Experiment in assay comments
        accession = remove_duplicates([a.get('comments', {}).get('ENA_EXPERIMENT', "") for a in assay_attributes])
        if len(accession) == 1:
            accession = accession[0]
        #else: report ERROR

        # Get all protocol refs
        protocolrefs = []
        for a in assay_attributes:
            protocolrefs.extend(a.get("protocol_ref"))
        protocolrefs.extend(extract_attributes.get("protocol_ref", []))
        protocolrefs = remove_duplicates(protocolrefs)

        # Get sample refs
        samplerefs = extract_attributes.get("sample_ref")

        # Get platform and instrument from sequencing protocol
        for p in protocols:
            if p.get("title") in protocolrefs and p.get("protocol type") == "nucleic acid sequencing protocol":
                hardware = p.get("hardware")
                lib_attribs["instrument_model"] = hardware
                # TODO: need to look up the platform type from ENA's controlled vocab
                # Using ILLUMINA as placeholder for now as it fits 90% of cases
                lib_attribs["platform_type"] = "ILLUMINA"
                # Expecting only one sequencing protocol per assay
                break

        return cls(alias, accession, technology_type[0], protocolrefs, samplerefs, lib_attribs)

    @classmethod
    def from_json(cls, assay_object):

        pass


class Protocol:
    def __init__(self, alias, accession, description, protocol_type, hardware, software, parameters):
        self.alias = alias
        self.accession = accession
        self.description = description
        self.protocol_type = protocol_type
        self.hardware = hardware
        self.software = software
        self.parameters = parameters

    def get_ae_attributes(self):
        """Return a list of all AE attributes that have values."""
        all_attributes = self.__dict__.keys()
        all_values = self.__dict__
        ae_attributes = ("protocol_type", "hardware", "software", "parameters")
        found_attributes = list()
        for a in all_attributes:
            if a in ae_attributes and all_values[a]:
                found_attributes.append(a)
        return found_attributes

    @classmethod
    def from_magetab(cls, protocol_dict):
        alias = protocol_dict.get("title")
        if is_accession(alias):
            accession = alias
        else:
            accession = None
        description = protocol_dict.get("description")
        protocol_type = protocol_dict.get("protocol_type")
        hardware = protocol_dict.get("hardware")
        software = protocol_dict.get("software")
        parameters = protocol_dict.get("parameters")

        return cls(alias, accession, description, protocol_type, hardware, software, parameters)


class Project:
    def __init__(self, alias, accession, title, description, releaseDate, publications, contacts):
        self.alias = alias
        self.accession = accession
        self.title = title
        self.description = description
        self.releaseDate = releaseDate
        self.publications = publications
        self.contacts = contacts

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.releaseDate})," \
               "{self.publications}, {self.contacts}".format(self=self)

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
        publications = study_info.get("publications", [])
        contacts = study_info.get("contacts", [])
        # Transform rolse to list
        for c in contacts:
            if c.get("roles"):
                roles = c["roles"]
                c["roles"] = roles.split(';')
        title = study_info.get("title")
        description = study_info.get("description")

        return cls(alias, accession, title, description, releaseDate, publications, contacts)


class Study:
    def __init__(self, alias, accession, title, description, protocolrefs, projectref,
                 experimental_factor, experimental_design, experiment_type, date_of_experiment,
                 comments=None):
        self.alias = alias
        self.accession = accession
        self.title = title
        self.description = description
        self.protocolrefs = protocolrefs
        self.projectref = projectref
        self.experimental_factor = experimental_factor
        self.experimental_design = experimental_design
        self.experiment_type = experiment_type
        self.date_of_experiment = date_of_experiment
        self.comments = comments

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.title})," \
               "{self.description}, {self.protocolrefs}, {self.projectref}, {self.experimental_factor}, " \
               "{self.experimental_design}, {self.experiment_type}, {self.comments}".format(self=self)

    @classmethod
    def from_magetab(cls, study_info):
        accession = study_info.get("accession")

        idf_file = os.path.basename(study_info.get("idf_filename", ""))
        alias = re.sub("\.idf\.txt$", "", idf_file)
        projectref = "project_" + alias

        title = study_info.get("title")
        description = study_info.get("description")
        experimental_factor = study_info.get("experimental_factor", [])
        experimental_design = study_info.get("experimental_design", [])
        protocolrefs = study_info.get("protocolRefs", [])

        date_of_experiment = study_info.get("date_of_experiment", None)
        comments = study_info.get("comments", {})
        experiment_type = comments.get("experiment_type", [])

        return cls(alias, accession, title, description, protocolrefs, projectref,
                   experimental_factor, experimental_design, experiment_type, date_of_experiment,
                   comments)


class AssayData:
    def __init__(self, alias, files, datatype, assayrefs, protocolrefs, accession=None):
        self.alias = alias
        self.files = files
        self.assayrefs = assayrefs
        self.datatype = datatype
        self.protocolrefs = protocolrefs
        self.accession = accession


class DataFile:
    def __init__(self, name, ftp_location=None, checksum=None, checksum_method="md5"):
        self.name = name
        self.checksum = checksum
        self.checksum_method = checksum_method
        self.ftp_location = ftp_location

    def __repr__(self):
        return "{self.__class__.__name__}({self.name}, {self.ftp_location}, {self.checksum})".format(self=self)

    @classmethod
    def from_sdrf(cls, file_attributes):
        name = file_attributes.get("name")

        comments = file_attributes.get("comments", {})
        checksum = comments.get("MD5")

        if "ArrayExpress FTP file" in comments:
            ftp_location = comments.get("ArrayExpress FTP file")
        elif "FASTQ_URI" in comments:
            ftp_location = comments.get("FASTQ_URI")
        else:
            ftp_location = None

        return cls(name, ftp_location, checksum)


class Contact:
    def __init__(self, firstName, lastName, email, affiliation, address, phone, roles):
        self.firstName = firstName
        self.lastName = lastName
        self.email = email
        self.affiliation = affiliation
        self.address = address
        self.phone = phone
        self.roles = roles


class Publication:
    def __init__(self, articleTitle, authors, pubmedId, doi, status):
        self.articleTitle = articleTitle
        self.authors = authors
        self.pubmedId = pubmedId
        self.doi = doi,
        self.status = status

