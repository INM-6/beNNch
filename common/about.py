#!/usr/bin/env python

import os

from pprint import pformat
import logging
log = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)

recordables = {
    'date': 'date --iso=seconds',
    'lshw': 'lshw -json -quiet',
    'dmidecode': 'dmidecode',
    'lspci': 'lspci -v',
    'false': '/bin/false',
    'broken': 'nothing',
    'cpuinfo': 'cat /proc/cpuinfo',
    'meminfo': 'cat /proc/meminfo',
    'env-vars': '/usr/bin/env',
    'ldd-nest': 'ldd nest',
    'conda-environment': 'conda env export',
    'hostname': 'hostname -f',
    'ompi_info': 'ompi_info',
    'ip-r': 'ip r',
    'ip-l': 'ip l',
}

from subprocess import Popen, PIPE, STDOUT, CalledProcessError, TimeoutExpired
import shlex

class Recorder(object):
    def __init__(self, outdir = "about", timeout=5, errors_fatal=False):
        self.errors_fatal = errors_fatal
        self.timeout = timeout
        self.outdir = outdir
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
            log.warning("created output directory %s", outdir)

    def record(self, recordables):
        for name, command in recordables.items():
            log.info("recording %s...", name)
            outname = os.path.join(self.outdir, name) + ".out"
            errname = os.path.join(self.outdir, name) + ".err"

            try:
                with Popen(shlex.split(command), stdout=PIPE, stderr=PIPE) as infile:
                    infile.wait(timeout=self.timeout)
                    if infile.returncode != 0:
                        log.warning("%s: returned %s (non-zero)!", name, infile.returncode)
                    with open(outname, 'wb') as outfile:
                        outfile.write(infile.stdout.read())
                    errs = infile.stderr.read()
                    if errs:
                        with open(errname, 'wb') as errfile:
                            log.warning("ERRORS recorded for %s", name)
                            errfile.write(errs)
                            if self.errors_fatal:
                                log.fatal("ERRORS are configured to be fatal.")
                                raise ValueError("Process wrote errors to STDERR!")
            except TimeoutExpired as e:
                log.warning("%s: process did not finish in time! Output will be incomplete!", name)
            except CalledProcessError as e:
                log.error("%s: called process failed! retrun code: %d", name, e.return_code)
            except FileNotFoundError as e:
                log.error("%s: %s", name, e)

def main():
    recorder = Recorder()
    recorder.record(recordables)

if __name__ == '__main__':
    main()