import json


def read_json_file(filename):
    try:
        with open(filename) as fh:
            data = json.load(fh, encoding="utf-8")
            return data
    except IOError as err:
        print "Cannot import file: %s" % err
        # TODO: system exit
    except ValueError as j_err:
        print "Cannot read JSON file: %s" % j_err
        raise

