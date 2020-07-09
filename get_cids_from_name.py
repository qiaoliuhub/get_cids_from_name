#! /usr/bin/env python

try:
    from urllib.error import HTTPError
    from urllib.parse import quote, urlencode
    from urllib.request import urlopen
except ImportError:
    from urllib import urlencode
    from urllib2 import quote, urlopen, HTTPError

import argparse
import time
import pandas as pd
import json
import xml.etree.ElementTree as ET

# use this url to save cid list in eutils server
NAMES_LISTKEY_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/cids/JSON?list_return=listkey"

def xml2df(root):

    cols = 'ID,Name,PharmaActions,InChIKey'
    df = pd.DataFrame(columns=cols.split(","))
    for i, drug in enumerate(root):

        info = {}
        id = drug.find("./*[@Name='CID']")
        info['ID'] = id.text if id is not None else None
        name = drug.find("./*[@Name='MeSHHeadingList']/*[@Name='string']")
        info['Name'] = name.text if name is not None else None
        info['PharmaActions'] = ";".join(x.text for x in drug.findall("./*[@Name='PharmActionList']/*[@Name='string']"))
        inchikey = drug.find("./*[@Name='InChIKey']")
        info['InChIKey'] = inchikey.text if inchikey is not None else None
        df.loc[i] = pd.Series(info)

    return df

def get_id(name):

    namespace = 'name'
    post_body = urlencode([(namespace, name)]).encode('utf8')
    esummary = None
    try:
        response = urlopen(NAMES_LISTKEY_API, post_body)
        print("successfully get list key result")
        # Construct esummary retrieve url
        lsit_key_result = json.loads(response.read())
        esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi" \
                   "?db=pccompound&WebEnv=" + lsit_key_result['IdentifierList']['EntrezWebEnv']\
                   + "&query_key=" + str(lsit_key_result['IdentifierList']['EntrezQueryKey'])

    except HTTPError as e:
        print("Fail to retrieve messages for {0!r}, caused by {1!r}".format(name, e))

    
    try:
        time.sleep(5)
        if esummary is not None:
            summary_response = urlopen(esummary)
            print("successfully get summary result")
            # Parsing the downloaded esummary xml string
            xml_result = summary_response.read().decode('utf-8')
            root = ET.fromstring(xml_result)
            for i, drug in enumerate(root):
                if i == 0:
                    id = drug.find("./*[@Name='CID']")
                    return id.text
    except HTTPError as e:
        print("Fail to retrieve summaries for {0!r}, caused by {1!r}".format(name, e))

    

if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('name_ls_df', help = 'the data frame saving names to be searched')
    argparser.add_argument('result_df', help = 'directory to save final results')
    # argparser.add_argument('xml_file', help = "directory to save requested xml results from pubchem server")
    args = argparser.parse_args()

    names = pd.read_csv(args.name_ls_df)
    names = names[names['pert_type'] == 'trt_cp']
    names = names[~names['pert_iname'].str.startswith('BRD')]
    names_list = set(names['pert_iname'].astype(str))
    print(len(names_list))
    result = pd.DataFrame(columns = ['name', 'cid'])

    # construct apiurl and add cids into POST body
    for i, name in enumerate(names_list):
        cid = get_id(name)
        result.loc[i] = pd.Series({'name': name, 'cid': cid})
        if i%100 == 0:
            result.to_csv(args.result_df)
    result.to_csv(args.result_df)
