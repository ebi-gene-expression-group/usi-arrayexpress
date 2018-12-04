import json
import requests
import time

__author__ = 'Ahmed G. Ali'

BASE_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def esearch(db, term, history=False):
    """
    Calling ``esearch`` end point of `eutils`
    :param db: Database to search for.
    :type db: str
    :param term: Search keyword.
    :type term: str
    :param history: Flag whether to use history ``WebEnv`` in next request or not. Defaule: `False`
    :type history: bool
    :return: Json object containing results as collected from the endpoint.
    :rtype: dict or list
    """
    term = term.replace('(', ' ').replace(')', ' ')
    url = BASE_URL + 'esearch.fcgi'
    data = {'db': db, 'term': term, 'retmode': 'json'}
    if history:
        data['usehistory'] = 'y'
    r = requests.get(url, params=data)
    if 'Error 503' in r.text:
        print('eutils gave Error 503. Waiting 20 secs then trying again')
        time.sleep(20)
        return esearch(db, term)
    return json.loads(r.text)


def efetch(db, ids):
    """
    Calling the ``efetch`` endpoint `eutils`

    :param db: Database to search for.
    :type db: str
    :param ids: list of ids to be fetched.
    :type ids: :obj:`list` of :obj:`str`
    :return: Json object as returned from the endpoint
    :rtype: dict or list
    """
    url = BASE_URL + 'efetch.fcgi'
    data = {'db': db, 'id': ','.join(ids)}
    r = requests.get(url, params=data)
    if 'Error 503' in r.text:
        print('eutils gave Error 503. Waiting 20 secs then trying again')
        time.sleep(20)

    return json.loads(r.text)


def esummary(db, query_id, web_env, ret_start=0, ret_max=500):
    """
    Calling the ``esummary`` endpoint of `eutils`

    :param db: Database to search for.
    :type db: str
    :param query_id: ``query_key`` parameters required by the endpoint and can be retrieved from previous searches
    :type query_id: int
    :param web_env: Web environment obtained from previous searches with flag for use history.
    :type web_env: str
    :param ret_start: Index of first return result.
    :type ret_start: int
    :param ret_max: Number of returned objects.
    :type ret_max: int
    :return: Json object containing results as collected from the endpoint.
    :rtype: dict or list
    """
    url = BASE_URL + 'esummary.fcgi'
    data = {'db': db, 'query_key': query_id, 'WebEnv': web_env, 'retmode': 'json', 'retstart': ret_start,
            'retmax': ret_max}
    r = requests.get(url, params=data)
    if 'Error 503' in r.text:
        print('eutils gave Error 503. Waiting 20 secs then trying again')
        time.sleep(20)
        return esummary(db, query_id, web_env, ret_start, ret_max)
    return json.loads(r.text)


if __name__ == '__main__':
    organism = "Mycoplasma gallisepticum str. R(high)"
    db = 'taxonomy'
    a = esearch(db=db, term=organism)
    # print a
    # a = esearch(db='gds', term='GPL[ETYP]', history=True)
    print(a)
    if a['esearchresult']['idlist']:
        taxonomy_id = int(a['esearchresult']['idlist'][0])
        print(taxonomy_id)
    else:
        raise Exception('Taxonomy ID for %s was not found on eutils. '
                        'Please try to check the organism in the SDRF' % organism)
    # print esummary(
    #     db='gds',
    #     query_id=a['esearchresult']['querykey'],
    #     web_env=a['esearchresult']['webenv']
    # )
