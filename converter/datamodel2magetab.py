"""Module to convert metadata in the submission data model to MAGE-TAB files."""

from collections import OrderedDict
from utils import converter_utils
import codecs
import csv


def generate_idf(sub):
    idf = OrderedDict([
        ("MAGE-TAB Version", "1.1"),
        ("Investigation Title", sub.study.title),
        ("Experiment Description", sub.study.description),
        ("Experimental Design", [d.value for d in sub.study.experimental_design]),
        ("Experimental Design Term Source REF", [d.term_source for d in sub.study.experimental_design]),
        ("Experimental Design Term Accession Number", [d.term_accession for d in sub.study.experimental_design]),
        ("Experimental Factor Name", [f.value for f in sub.study.experimental_factor]),
        ("Experimental Factor Type", [f.value for f in sub.study.experimental_factor]),
        ("Experimental Factor Term Source REF", [f.term_source for f in sub.study.experimental_factor]),
        ("Experimental Factor Term Accession Number", [f.term_accession for f in sub.study.experimental_factor]),
        ("Person Last Name", [p.lastName for p in sub.project.contacts]),
        ("Person First Name", [p.firstName for p in sub.project.contacts]),
        ("Person Mid Initials", [p.middleInitials for p in sub.project.contacts]),
        ("Person Email", [p.email for p in sub.project.contacts]),
        ("Person Phone", [p.phone for p in sub.project.contacts]),
        ("Person Fax", [p.fax for p in sub.project.contacts]),
        ("Person Address", [p.address for p in sub.project.contacts]),
        ("Person Affiliation", [p.affiliation for p in sub.project.contacts]),
        ("Person Roles", [';'.join(p.roles) for p in sub.project.contacts]),
        ("Date of Experiment", sub.study.date_of_experiment),
        ("Public Release Date", sub.project.releaseDate),
        ("PubMed ID", [p.pubmedId for p in sub.project.publications]),
        ("Publication DOI", [p.doi for p in sub.project.publications]),
        ("Publication Author List", [p.authors for p in sub.project.publications]),
        ("Publication Title", [p.articleTitle for p in sub.project.publications]),
        ("Publication Status", [p.status for p in sub.project.publications]),
        ("Publication Status Term Source REF", ["EFO" for p in sub.project.publications if p.status]),
        ("Publication Status Term Accession Number", []),
        ("Protocol Name", [p.alias for p in sub.protocol]),
        ("Protocol Type", [p.protocol_type.value for p in sub.protocol]),
        ("Protocol Term Source REF", [p.protocol_type.term_source for p in sub.protocol]),
        ("Protocol Term Accession Number", [p.protocol_type.term_accession for p in sub.protocol]),
        ("Protocol Description", [p.description for p in sub.protocol]),
        ("Protocol Hardware", [p.hardware for p in sub.protocol]),
        ("Protocol Software", [p.software for p in sub.protocol]),
        ("SDRF File", sub.info.get("accession") + ".sdrf.txt"),
        ("Term Source Name", ["EFO", "ArrayExpress"]),
        ("Term Source File", ["https://www.ebi.ac.uk/efo/", "https://www.ebi.ac.uk/arrayexpress/"]),
        ("Comment[AEExperimentType", [exptype for exptype in sub.study.experiment_type]),
        ("Comment[ArrayExpressAccession]", sub.info.get("accession"))
    ])

    if sub.info.get("submission_type") == "sequencing" or sub.info.submission_type == "singlecell":
        pass
        # TODO: Add comments for ENA accessions Comment[SecondaryAccession] and Comment[SequenceDataURI]

    return idf


def generate_sdrf(sub):
    """Transform metadata in data model to an SDRF table."""

    header = []
    rows = []

    if sub.info.get("submission_type") == "microarray":

        for sample in sub.sample:

            sample_values = [
                ("Source Name", sample.alias),
                ("Description", sample.description),
                ("Material Type", sample.material_type)
            ]
            print(sample_values)
            # Get all assay objects that belong to this sample
            assays = [assay for assay in sub.assay if assay.sampleref == sample.alias]

            for assay in assays:

                assay_values = [
                    ("Extract Name", sample.alias),
                    ("Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "nucleic acid labeling protocol"]),
                    ("Labelled Extract Name", assay.alias),
                    ("Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "nucleic acid hybridization to array protocol"])
                ]
                print(assay_values)

                # Get all assay data objects that belong to this assay
                data = [ad for ad in sub.assay_data if assay.alias in ad.assayrefs]

                for ad in data:
                    data_values = [
                        ("Assay Name", ad.alias),
                        ("Technology Type", assay.technology_type),
                        ("Array Design REF", assay.array_design),
                        ("Term Source REF", "Array Express"),
                        ("Protocol REF", [sub.get_protocol(pref).alias for pref in ad.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "array scanning and feature extraction protocol"]),
                    ]

                    # Get all data files
                    for f in ad.files:
                        if ad.data_type == "raw":
                            data_values.append(("Array Data File", f.name))
                        elif ad.data_type == "raw matrix":
                            data_values.append(("Array Data Matrix File", f.name))

                    print(sample_values + assay_values + data_values)

                    rows.append(sample_values + assay_values + data_values)

    return rows



