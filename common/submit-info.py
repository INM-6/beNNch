#!/usr/bin/env python
# encoding: utf8
'''
Usage: submit-info <logfile>

   Extract the block of essential information from JUBE run command output and
   write it to stdout.

'''
from pprint import pformat
import logging
log = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)
import sys
def main():
    filename = sys.argv[1]
    log.info("extracting JUBE info from %s", filename)

    info = {}
    with open(filename, 'r') as infile:
        #print("jube_run:")
        for line in infile:
            if not line.startswith(">>>> "): continue
            if line.strip().endswith(":"): continue
            line = line[5:].strip()
            #print("    " + line)
            key, val = line.split(":", 1)
            info[key.strip()] = val.strip()
    print("{handle} --id {id}".format(**info))
    log.info("identified handle '%s' (id %s)", info['handle'], info['id'])

if __name__ == '__main__':
    main()
