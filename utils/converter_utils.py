import json
import pkg_resources
import re
import requests
import urllib3.request
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def read_json_file(filename):
    try:
        with open(filename) as fh:
            data = json.load(fh, encoding="utf-8")
            return data
    except IOError as err:
        print("Cannot import file: %s" % err)
    except ValueError as j_err:
        print("Cannot read JSON file: %s" % j_err)
        raise


def ontology_lookup(category):
    """Read the json with expected EFO terms and return the dict for the given category."""
    resource_package = __name__
    resource_path = "ontology_terms.json"
    all_terms = json.loads(pkg_resources.resource_string(resource_package, resource_path))

    return all_terms[category]


def get_controlled_vocabulary(category):
    """Read the json with controlled vocab and return the dict for the given category."""
    resource_package = __name__
    resource_path = "term_translations.json"
    all_terms = json.loads(pkg_resources.resource_string(resource_package, resource_path))

    return all_terms[category]


def remove_duplicates(ref_list):
    """Return a copy of a list with all duplicated values removed."""
    return list(set(ref_list))


def is_accession(accession, archive=None):
    """Return True if the input is a valid accession format from specified EBI archives.
    With the optional parameter the test can be performed against a specific archive only."""

    regex_lookup = {
        "ARRAYEXPRESS": "^[A-Z]-[A-Z]{4}-[0-9]+",
        "BIOSAMPLES": "^SAM[END][AG]?[0-9]+",
        "ENA": "^ER[RXSP][0-9]+$",
        "BIOSTUDIES": "^S-[A-Z]+[0-9]+$"}

    regex_ebi_accession = "|".join(regex_lookup.values())

    if archive:
        try:
            regex = regex_lookup.get(archive)
            return re.match(regex, accession)
        except KeyError:
            print("Not a valid accession type: {}".format(archive))
    else:
        return re.match(regex_ebi_accession, accession)


def get_url(url):
    try:
        r = requests.get(url)
        data = json.loads(r.text)
        return data
    except requests.HTTPError as e:
        print('HTTPError when retrieving url: "' + url + '" : ' + str(e.errno))
    except requests.ConnectionError as e:
        print('ConnectionError when retrieving url: "' + url + '" : ' + str(e.errno))
    except json.JSONDecodeError:
        print('JSONDecodeError: The result is not a valid JSON response: ' + url)
    return None


# To store organisms that we have already looked-up in the taxonomy (this is slow...)
organism_lookup = {}


def get_taxon(organism):
    """Return the NCBI taxonomy ID for a given species name."""
    if organism and organism not in organism_lookup:
        url = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/suggest-for-search/' + organism.replace(
                ' ', '%20') + '?limit=200'
        organism_data = get_url(url)
        if not organism_data:
            if re.search(r" and | \+ ", organism):
                # It looks as if we have more than one organism mixed in one sample - in the case assign the 'mixed
                # sample' taxon_id (c.f. https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1427524) - as
                # per instructions on https://www.ebi.ac.uk/seqdb/confluence/display/GXA/Curation+Look-up
                return 1427524
            else:
                print("Failed to retrieve organism data from ENA taxonomy service for: " + organism)
        else:
            for organism_entry in organism_data:
                if organism_entry['scientificName'].lower() == organism.lower():
                    taxon_id = int(organism_entry['taxId'])
                    organism_lookup[organism] = taxon_id
                    return taxon_id
    else:
        return organism_lookup.get(organism)
