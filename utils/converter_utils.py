import json
import pkg_resources
import re
from pip._vendor import urllib3


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
    regex_lookup = {
        "ARRAYEXPRESS": "^[A-Z]-[A-Z]{4}-[0-9]+",
        "BIOSAMPLES": "^SAMEA[0-9]+$",
        "ENA": "^ER[RESP][0-9]+$",
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



# Cannot import "common" from fgsubs; it has too many dependencies, so copying this here
def get_url(url):
    data = None
    try:
        r = urllib3.urlopen(url)
        data = json.load(r)
    except urllib2.HTTPError as e:
        print('HTTPError when retrieving url: "' + url + '" : ' + str(e.code))
    except urllib2.URLError as e:
        print('URLError when retrieving url: "' + url + '" : ' + str(e.reason))
    else:
        return data


def get_taxon(organism):
    taxon_to_id = {}
    if organism != '':
        # TODO: Fix urllib.urlencode instead of organism.replace(' ','%20')
        organism_data = get_url(
            'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/suggest-for-search/' + organism.replace(
                ' ', '%20') + '?limit=200')
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
                    return taxon_id
