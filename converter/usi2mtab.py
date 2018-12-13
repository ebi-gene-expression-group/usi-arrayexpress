
import json
import re
from collections import OrderedDict


def create_idf(sdrf_file_name, project, study, protocols, translations):

    # Initialising an IDF dictionary
    idf = OrderedDict([("MAGE-TAB Version", u"1.1")])

    # Title and description (mandatory fields)
    idf["Investigation Title"] = study["title"]
    idf["Experiment Description"] = study["description"]

    # Optional fields about experimental design
    designs = study["attributes"].get("experimental design", "")
    if designs:  # a list of all design terms with ontology refs
        idf["Experimental Design"] = [t["value"] for t in designs]  # List comprehension picking up the value of terms
        idf["Experimental Design Term Source REF"] = ["EFO" for t in designs]  # This is currently not part of the json so hardcoding EFO here
        idf["Experimental Design Term Accession Number"] = [t["terms"][0]["url"] for t in designs]  # The schema allows >1 urls per term, expecting only one
        # Extracting EFO accession from URL
        url_to_efo_accession(idf["Experimental Design Term Accession Number"])
    # Experimental factors (mandatory)
    factors = study["attributes"].get("experimental factor", "")
    idf["Experimental Factor Name"] = [t["value"] for t in factors]
    idf["Experimental Factor Type"] = idf["Experimental Factor Name"]  # The factor type should is inferred from the factor name (should always be the same)

    # Optional experimental factor terms
    idf["Experimental Factor Term Source REF"] = ["EFO" if f["terms"] else "" for f in factors]  # Expecting all terms are from EFO
    idf["Experimental Factor Term Source Accession Number"] = [f["terms"][0]["url"] if f["terms"] else "" for f in factors]
    # Extracting EFO accession from URL
    url_to_efo_accession(idf["Experimental Factor Term Source Accession Number"])

    # Person details (mandatory fields)
    persons = project["contacts"]
    # Get IDF field dict
    idf_contact_terms = translations["contact_terms"]
    # Look up values for each term and create a list
    for idf_term, usi_term in idf_contact_terms.items():
        idf[idf_term] = [p[usi_term] if p.get(usi_term, "") else "" for p in persons]

    # Roles need to be concatenated to a single string and terms added
    list_to_semicol_string(idf["Person Roles"])

    # Optional date of experiment
    if "date of experiment" in study["attributes"].keys():
        idf["Date of Experiment"] = study["attributes"]["date of experiment"][0]["value"]

    # Release date (mandatory)
    idf["Public Release Date"] = project["releaseDate"]

    # Optional publication details
    if "publications" in project.keys():
        publications = project["publications"]
        for idf_term, usi_term in translations["publication_terms"].items():
            idf[idf_term] = [p[usi_term] if p.get(usi_term, "") else "" for p in publications]
    # TODO: Add "Publication Status Term Source REF" and "Publication Status Term Accession Number"

    # Protocols (mandatory)
    idf["Protocol Name"] = [p["title"] for p in protocols]

    idf["Protocol Type"] = [p["attributes"]["protocol_type"][0]["value"] for p in protocols]
    idf["Protocol Term Source REF"] = [p["attributes"]["protocol_type"][0]["terms"][0]["url"] for p in protocols]
    url_to_efo_accession(idf["Protocol Term Source REF"])

    protocol_terms = translations["protocol_terms"]
    del protocol_terms["Protocol Type"]  # Remove protocol_type because this is already done
    for idf_term, usi_term in protocol_terms.items():
        # Adding get value check here so that this works even if field is not in the json
        idf[idf_term] = [p["attributes"][usi_term][0]["value"] if p["attributes"].get(usi_term, "") else "" for p in protocols]

    # Ontology reference (mandatory if used)
    idf["Term Source Name"] = ["EFO", "ArrayExpress"]
    idf["Term Source File"] = ["http://www.ebi.ac.uk/efo/", "http://www.ebi.ac.uk/arrayexpress"]
    idf["Term Source Version"] = ["", ""]

    # SDRF reference (mandatory)
    idf["SDRF File"] = sdrf_file_name

    # ArrayExpress experiment type (mandatory)
    exp_type_term = translations["comment_terms"]["AEExperimentType"]
    idf["Comment[AEExperimentType]"] = study["attributes"][exp_type_term][0]["value"]

    # This is a dictionary of IDF fields as keys and the values as list (where there is more than 1 item)
    return idf


# Move to common?
def read_json_file(translations_file):
    with open(translations_file) as tf:
        translations_json = json.load(tf,  object_pairs_hook=OrderedDict)
    return translations_json


def url_to_efo_accession(term_list):
    # This takes in a list of OLS term urls and modifies it so that the list contains EFO accessions instead
    for i, url in enumerate(term_list):
        if re.search('efo', url, flags=re.IGNORECASE):
            efo_acc = url.split('/')[-1]
            term_list[i] = efo_acc
        else:
            # TODO: Look up EFO term via API
            pass  # For now we'll keep the whole URL


def list_to_semicol_string(term_list):
    for i, tl in enumerate(term_list):
        term_list[i] = ';'.join(tl)

