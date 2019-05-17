
import os
import re
from collections import OrderedDict

from utils.converter_utils import is_accession, get_controlled_vocabulary, remove_duplicates, get_taxon


class Sample:
    def __init__(self, alias, accession, taxon, taxonId, attributes, material_type, description):
        """
        :param alias: string, unique sample name in the experiment
        :param accession: string, BioSamples accession
        :param taxon: string, latin species name
        :param taxonId: int, NCBI taxonomy identifier for the species
        :param attributes: dictionary of attribute categories as keys and Attribute class object as value
        :param material_type: string, (optional) one of: whole organism, organism part, cell, RNA, DNA
        :param description: string, (optional) free-text description of sample
        """
        self.alias = alias
        self.accession = accession
        self.taxon = taxon
        self.taxonId = taxonId
        self.material_type = material_type
        self.description = description
        self.attributes = attributes

    @classmethod
    def from_magetab(cls, sample_attributes):
        alias = sample_attributes.get("name")
        description = sample_attributes.get("description")
        material_type = sample_attributes.get("material_type")

        comments = sample_attributes.get("comments")
        accession = comments.get("BioSD_SAMPLE")

        characteristics = sample_attributes.get("characteristics")
        factors = sample_attributes.get("factors")
        organism = characteristics.get("organism", {})
        taxon = organism.get("value")
        taxonId = get_taxon(taxon)

        # Note this will overwrite the characteristics values if a factor is also a characteristics
        raw_attributes = characteristics.copy()
        raw_attributes.update(factors)

        attributes = OrderedDict()
        for c_name, c_attrib in raw_attributes.items():
            new_unit = None
            if "unit" in c_attrib:
                unit_attrib = c_attrib.get("unit")
                new_unit = Unit(unit_attrib.get("value"),
                                unit_attrib.get("unit_type"),
                                unit_attrib.get("term_accession"),
                                unit_attrib.get("term_source"))

            attributes[c_name] = Attribute(c_attrib.get("value"),
                                           new_unit,
                                           c_attrib.get("term_accession"),
                                           c_attrib.get("term_source"))

        return cls(alias, accession, taxon, taxonId, attributes, material_type, description)


class Assay:
    def __init__(self, alias, accession, technology_type, protocolrefs, sampleref):
        """
        Attributes common to all assay types
        :param alias: string, unique name of the assay
        :param accession: string, (not for microarray) ENA experiment accession
        :param technology_type: string, one of: array assay, sequencing assay
        :param protocolrefs: list, protocol accessions/name used to generate assay data
        :param sampleref: string, sample accession/name
        """
        self.alias = alias
        self.accession = accession
        self.technology_type = technology_type
        self.protocolrefs = protocolrefs
        self.sampleref = sampleref

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.technology_type}, " \
               "{self.protocolrefs}, {self.sampleref})".format(self=self)

    def get_attributes(self):
        """Return a list of all additional attributes that have values"""
        all_attributes = self.__dict__.keys()
        all_values = self.__dict__
        fixed_attributes = ("alias", "accession", "technology_type", "protocolrefs", "sampleref")
        other_attributes = [a for a in all_attributes if a not in fixed_attributes and all_values[a]]
        return other_attributes


class MicroarrayAssay(Assay):
    def __init__(self, alias, accession, technology_type, protocolrefs, sampleref, label, array_design):
        """
        Assay attributes specific to microarray assays

        :param label: string, the label type
        :param arrayref: string, ArrayExpress array design format accession
        """
        Assay.__init__(self, alias, accession, technology_type, protocolrefs, sampleref)
        self.label = label
        self.array_design = array_design

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
        array_design = [a.get("array_design") for a in assay_attributes]

        return cls(alias, accession, technology_type[0], protocolrefs, sampleref, label, array_design[0])

    @classmethod
    def from_json(self):
        pass


class SeqAssay(Assay):
    def __init__(self, alias, accession, technology_type, protocolrefs, sampleref, lib_attribs):
        """
        Assay attributes specific to sequencing assays

        :param lib_attribs: dictionary of sequencing specific attributes that are required for ENA brokering
                            keys: library attribute
                            values: values from SRA's controlled vocabulary
        """
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
            if p.get("title") in protocolrefs and p.get("protocol_type") == "nucleic acid sequencing protocol":
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
        """
        Attributes of a protocol
        :param alias: string, protocol accession or unique name in experiment
        :param accession: string, ArrayExpress protocol accession
        :param description: string, free-text description of protocol
        :param protocol_type: object, Attribute class object, experiment type ontology term
        :param hardware: string, free-text hardware description, sequencer model for sequencing experiments
        :param software: string, free-text software description
        :param parameters: list, parameter values
        """
        self.alias = alias
        self.accession = accession
        self.description = description
        self.protocol_type = protocol_type
        self.hardware = hardware
        self.software = software
        self.parameters = parameters

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.description}, " \
               "{self.protocol_type}, {self.hardware}, {self.software}, {self.parameters})".format(self=self)

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
        hardware = protocol_dict.get("hardware")
        software = protocol_dict.get("software")
        parameters = protocol_dict.get("parameters")

        protocol_type = Attribute(protocol_dict.get("protocol_type"), None,
                                  protocol_dict.get("term_accession"),
                                  protocol_dict.get("term_source"))

        return cls(alias, accession, description, protocol_type, hardware, software, parameters)


class Project:
    def __init__(self, alias, accession, title, description, releaseDate, publications, contacts):
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
        self.alias = alias
        self.accession = accession
        self.title = title
        self.description = description
        self.releaseDate = releaseDate
        self.publications = publications
        self.contacts = contacts

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.releaseDate}," \
               "{self.publications}, {self.contacts})".format(self=self)

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

        return cls(alias, accession, title, description, releaseDate, publications, contacts)


class Study:
    def __init__(self, alias, accession, title, description, protocolrefs, projectref,
                 experimental_factor, experimental_design, experiment_type, date_of_experiment,
                 comments=None):
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
        :param comments: dictionary, (optional) known IDF comments
                        keys: string, comment category
                        values: list, comment values
        """
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
        return "{self.__class__.__name__}({self.alias}, {self.accession}, {self.title}," \
               "{self.description}, {self.protocolrefs}, {self.projectref}, {self.experimental_factor}, " \
               "{self.experimental_design}, {self.experiment_type}, {self.comments})".format(self=self)

    @classmethod
    def from_magetab(cls, study_info):
        accession = study_info.get("accession")

        idf_file = os.path.basename(study_info.get("idf_filename", ""))
        alias = re.sub("\.idf\.txt$", "", idf_file)
        projectref = "project_" + alias

        title = study_info.get("title")
        description = study_info.get("description")
        ef = study_info.get("experimental_factor", [])

        ef_objects = [Attribute(d.get("experimental_factor"), None,
                                d.get("term_accession"),
                                d.get("term_source")) for d in ef]
        ed = study_info.get("experimental_design", [])
        ed_objects = [Attribute(d.get("experimental_design"), None,
                                d.get("term_accession"),
                                d.get("term_source")) for d in ed]
        protocolrefs = study_info.get("protocolRefs", [])

        date_of_experiment = study_info.get("date_of_experiment", None)
        comments = study_info.get("comments", {})
        experiment_type = comments.get("experiment_type", [])

        return cls(alias, accession, title, description, protocolrefs, projectref,
                   ef_objects, ed_objects, experiment_type, date_of_experiment,
                   comments)


class AssayData:
    def __init__(self, alias, files, data_type, assayrefs, protocolrefs, accession=None):
        """
        Attributes of the assay data
        :param alias: string, unique name in the experiment (auto-generated from Assay Name in SDRF)
        :param files: list, DataFile class objects
        :param data_type: string, raw or raw matrix
        :param assayrefs: list, assay name/accessions
        :param protocolrefs: list, protocol name/accessions used to generate data file
        :param accession: string, (not for microarray) ENA run accession
        """
        self.alias = alias
        self.files = files  # List of DataFile objects
        self.data_type = data_type  # "raw" or "raw matrix"
        self.assayrefs = assayrefs
        self.protocolrefs = protocolrefs
        self.accession = accession

    def __repr__(self):
        return "{self.__class__.__name__}({self.alias}, {self.files}, {self.data_type}, " \
               "{self.assayrefs}, {self.protocolrefs}, {self.accession})".format(self=self)

    @classmethod
    def from_magetab(cls, name, datafile_objects, file_attributes):

        alias = name

        # Assuming here that attributes of grouped files are the same
        common_file_attributes = file_attributes[0]

        # Try getting le_ref first for MA
        if common_file_attributes.get("le_ref"):
            assayrefs = remove_duplicates(common_file_attributes.get("le_ref"))
        else:
            assayrefs = remove_duplicates(common_file_attributes.get("extract_ref"))

        mage_data_type = common_file_attributes.get("data_type")
        if mage_data_type == "arraydatafile" or mage_data_type == "scanname":
            data_type = "raw"
        elif mage_data_type == "arraydatamatrixfile":
            data_type = "raw matrix"
        else:
            data_type = None

        comments = common_file_attributes.get("comments")
        accession = comments.get("ENA_RUN")

        protocolrefs = common_file_attributes.get("protocol_ref")

        return cls(alias, datafile_objects, data_type, assayrefs, protocolrefs, accession)


class Analysis:
    def __init__(self, alias, files, data_type, assaydatarefs, protocolrefs):
        """
        :param alias: string, unique name in the experiment (auto-generated from processed file name in SDRF)
        :param files: list, DataFile class objects
        :param data_type: string, processed or processed matrix
        :param assaydatarefs: list, assay data name/accessions
        :param protocolrefs: list, protocol name/accessions used to generate data file
        """
        self.alias = alias
        self.files = files
        self.data_type = data_type
        self.assaydatarefs = assaydatarefs
        self.protocolrefs = protocolrefs

    @classmethod
    def from_magetab(cls, datafile_objects, file_attributes):

        alias = file_attributes.get("name")

        mage_data_type = file_attributes.get("data_type")
        if mage_data_type == "derivedarraydatamatrixfile":
            data_type = "processed matrix"
        else:
            data_type = "processed"

        protocolrefs = file_attributes.get("protocol_ref")

        # Assay ref is generated based on le_ref (for MA) or extract_ref (for HTS)
        if file_attributes.get("le_ref"):
            assaydatarefs = remove_duplicates(file_attributes.get("le_ref"))
        else:
            assaydatarefs = remove_duplicates(file_attributes.get("extract_ref"))

        return cls(alias, datafile_objects, data_type, assaydatarefs, protocolrefs)


# Helper classes

class DataFile:
    def __init__(self, name, ftp_location=None, checksum=None, checksum_method=None):
        self.name = name
        self.checksum = checksum
        self.checksum_method = checksum_method
        self.ftp_location = ftp_location

    def __repr__(self):
        return "{self.__class__.__name__}({self.name}, {self.ftp_location}, {self.checksum})".format(self=self)

    @classmethod
    def from_magetab(cls, file_attributes):

        name = file_attributes.get("name")
        comments = file_attributes.get("comments", {})
        checksum = comments.get("MD5")

        if "ArrayExpress FTP file" in comments:
            ftp_location = comments.get("ArrayExpress FTP file")
        elif "FASTQ_URI" in comments:
            ftp_location = comments.get("FASTQ_URI")
        elif "Derived ArrayExpress FTP file" in comments:
            ftp_location = comments.get("Derived ArrayExpress FTP file")
        else:
            ftp_location = None

        return cls(name, ftp_location, checksum)


class Contact:
    def __init__(self, firstName, lastName, email, affiliation, address, phone, roles, middleInitials, fax):
        self.firstName = firstName
        self.lastName = lastName
        self.email = email
        self.affiliation = affiliation
        self.address = address
        self.phone = phone
        self.roles = roles
        self.middleInitials = middleInitials
        self.fax = fax

    def __repr__(self):
        return "{self.__class__.__name__}({self.firstName}, {self.lastName}, {self.email}, " \
               "{self.affiliation}, {self.address}, {self.phone}, {self.roles}, " \
               "{self.middleInitials}, {self.fax})".format(self=self)


class Publication:
    def __init__(self, articleTitle, authors, pubmedId, doi, status):
        self.articleTitle = articleTitle
        self.authors = authors
        self.pubmedId = pubmedId
        self.doi = doi
        self.status = status


class Attribute:
    def __init__(self, value, unit, term_accession, term_source="EFO"):
        self.value = value
        self.unit = unit
        self.term_source = term_source
        self.term_accession = term_accession

    def __repr__(self):
        return "{self.__class__.__name__}({self.value}, {self.unit}, {self.term_accession})".format(self=self)


class Unit:
    def __init__(self, value, unit_type, term_accession, term_source="EFO"):
        self.value = value
        self.unit_type = unit_type
        self.term_source = term_source
        self.term_accession = term_accession

    def __repr__(self):
        return "{self.__class__.__name__}({self.value}, {self.unit_type}, {self.term_accession})".format(self=self)


class Submission:
    def __init__(self, info, project, study, protocol, sample, assay, assay_data, analysis):
        """
        The submission class holds all objects that make up an ArrayExpress experiment

        :param info: dictionary, information about the submission (file names, submitter details, etc.)
        :param project: object, Project class object
        :param study: object, Study class object
        :param protocol: list, Protocol class objects
        :param sample: list, Sample class objects
        :param assay: list, SeqAssay or MicroarrayAssay class objects
        :param assay_data: list, AssayData class objects
        :param analysis: list, Analysis class objects
        """
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


