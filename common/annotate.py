#!/usr/bin/env python
# encoding: utf8
'''
Usage: annotate [options] [<basedir>]

   Genereate metadata annotations for all about/ directories. If <basedir> is
   not given the current directory is used as base directory.

Options:
    -v, --verbose       increase output
    -h, --help          print this text
'''
from docopt import docopt
import os, re
from fnmatch import fnmatch

from pprint import pformat
import logging
log = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)


def annotate(basedir):
    log.info("annotate %s", basedir)

    def wildcard(relpath, pattern):
        path = os.path.join(basedir, relpath)   # closure of basedir
        choice = [os.path.join(relpath, name) for name in os.listdir(path) if fnmatch(name, pattern)]
        assert len(choice) == 1, "wildcard matches must produce single result"
        return choice[0]

    submitlog = wildcard("../../../logs", "*.submit.log")
    annotations = {
        "cpu_model": ("about/cpuinfo.out", r'model name\s*:\s*(.*)'), # model of the compute node CPUs (assumes fully symmetric compute)
        "cpu_MHz":   ("about/cpuinfo.out", r'cpu MHz\s+:\s+([\d\.]+)'), # exemplary core frequency (dependant on cpufreq govenor)
        "cache_L1_size": ("about/lstopo", r'L1d? \((\d+)KB\)'), # single L1d cache size as reported by lstopo (usually per core)
        "cache_L2_size": ("about/lstopo", r'L2d? \((\d+)KB\)'), # single L2 cache size as reported by lstopo (may be shared among cores)
        "cache_L3_size": ("about/lstopo", r'L3d? \((\d+)KB\)'), # single L3 cache (usually shared among number of cores)
        "cache_L1_assoc": ("about/hwloc-topology.xml", r'depth="1" cache_linesize="\d+" cache_associativity="(\d+)"'), # cache associations *FIXME: more description*
        "cache_L2_assoc": ("about/hwloc-topology.xml", r'depth="2" cache_linesize="\d+" cache_associativity="(\d+)"'), # cache associations *FIXME: more description*
        "cache_L3_assoc": ("about/hwloc-topology.xml", r'depth="3" cache_linesize="\d+" cache_associativity="(\d+)"'), # cache associations *FIXME: more description*
        "openmpi_version": ("about/ompi_info.out", r'Open MPI: ([^\s]+)'), # if(mpi==openmpi) the version of the openmpi library
        "openmpi_revision": ("about/ompi_info.out", r'Open MPI repo revision: ([^\s]+)'),  # if(mpi==openmpi) the repository revision the openmpi library was built from
        "openmpi_api":      ("about/ompi_info.out", r'MPI API: ([^\s]+)'), # if(mpi==openmpi) the supported MPI standard
        "mpi": ("about/ldd-nest.out", r'lib((mpich|.*mpi)[^\s]+)'), # the mpi type used in the simulation (at runtime)
        "author": ("about/env-vars.out", r'USER=(\w+)'),    # author of the benchmark
        "hostname": ("about/hostname.out", r'(\w+)'),   # hostname where the job started to run
        "machine": (submitlog, r'machine/(\w+)'),   # machine configuration used for the run
        "benchmark": (submitlog, r'# benchmark: (\w+)'), # benchmark type configured to run
        "model": (submitlog, r'path: model/(\w+)'), # model configured to be benchmarked
        "node_memory": ("about/meminfo.out", r'MemTotal:\s*(\d+\s+kB)'), # total memory of the node
        "corespernode": ("about/nproc.out", r'\d+'), # nproc reported number of cores
    }
    cmd = "git annex metadata {basedir} --set '{field}={value}'"

    metadata = dict()
    for field, (filename, expr_re) in annotations.items():
        expr = re.compile(expr_re)
        try:
            with open(os.path.join(basedir, filename), 'r') as infile:
                match = expr.search(infile.read())  # get first match
                if match is None:
                    log.warning("no match found for %s in %s: regex was %s", field, filename, repr(expr_re))
                    continue
                metadata[field] = match.group(1)
        except FileNotFoundError as e:
            log.warning("extracting %s: %s", field, e)

    for key, val in metadata.items():
        print(cmd.format(
            basedir=basedir,
            field=key,
            value=val,
        ))


def annotate_all(basedir, callfolder="about"):
    '''
    calls annotate for all paths that contain a <callfolder> subdirectory.
    '''
    log.info("looking for annotation folders in %s", basedir)
    for path, dirs, files in os.walk(basedir):
        if callfolder in dirs:
            annotate(path)
            dirs.clear()    # do not traverse other directories here

def main():
    args = docopt(__doc__)
    if args['--verbose']:
        log.setLevel(logging.DEBUG)
    log.debug(pformat(args))

    log.info("Hello World")
    if not args['<basedir>']:
        args['<basedir>'] = os.curdir

    annotate_all(args['<basedir>'])

if __name__ == '__main__':
    main()
