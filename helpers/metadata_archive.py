#!/usr/bin/env python
# encoding: utf8

# beNNch - Unified execution, collection, analysis and
# comparison of neural network simulation benchmarks.
# Copyright (C) 2021 Forschungszentrum Juelich GmbH, INM-6

# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

# SPDX-License-Identifier: GPL-3.0-or-later

import os
import sys
import time
import logging
import logging.config
from subprocess import Popen, PIPE, DEVNULL, CalledProcessError, TimeoutExpired
import shlex


save_path = sys.argv[1]
log = logging.getLogger()

try:
    import yaml
    basepath = os.path.dirname(__file__)
    with open(os.path.join(basepath, "../logging.yaml"), 'r') as infile:
        logging.config.dictConfig(yaml.safe_load(infile))
except Exception as e:
    logging.basicConfig(level=logging.DEBUG)
    log.warning("using basic logging config due to exception %s", e)

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
    'nproc': 'nproc',
    'hwloc-info': 'hwloc-info',
    'hwloc-ls': 'hwloc-ls',
    # 'hwloc-topology': 'hwloc-gather-topology {outdir}/hwloc-topology',
    'lstopo': 'lstopo --of ascii {outdir}/{name}',
    'getconf': 'getconf -a',
    'ulimit': 'ulimit -a',
    'modules': 'module list',
}


class Recorder(object):
    def __init__(self, outdir="about", timeout=3, errors_fatal=False):
        self.errors_fatal = errors_fatal
        self.timeout = timeout
        self.logtimethres = 10  # seconds
        self.outdir = outdir or '.'
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
            log.warning("created output directory %s", outdir)

    def record(self, recordables):
        for name, command in recordables.items():
            log.info("recording %s...", name)
            outname = os.path.join(self.outdir, name) + ".out"
            errname = os.path.join(self.outdir, name) + ".err"

            parameters = {
                'outdir': self.outdir,
                'name': name,
                'command': command,
            }
            try:
                starttime = time.time()
                stoptime = None
                iotime = None
                with Popen(shlex.split(command.format(**parameters)),
                           stdout=PIPE, stderr=PIPE, stdin=DEVNULL) as infile:
                    try:
                        (stdout_data, stderr_data) = infile.communicate(
                            timeout=self.timeout)
                    except TimeoutExpired:
                        log.warning(
                            "%s: process did not finish in time! Output will be"
                            "incomplete!", name)
                        infile.kill()
                        outs, errs = infile.communicate()
                        log.error("Final words on stdout:\n%s", outs)
                        log.error("Final words on stderr:\n%s", errs)
                    stoptime = time.time()
                    if infile.returncode != 0:
                        log.warning("%s: returned %s (non-zero)!",
                                    name, infile.returncode)
                    with open(outname, 'wb') as outfile:
                        outfile.write(stdout_data)
                    if stderr_data:
                        with open(errname, 'wb') as errfile:
                            log.warning("ERRORS recorded for %s", name)
                            errfile.write(stderr_data)
                            if self.errors_fatal:
                                log.fatal("ERRORS are configured to be fatal.")
                                raise ValueError(
                                    "Process wrote errors to STDERR!")
                    iotime = time.time()
            except CalledProcessError as e:
                log.error("%s: called process failed! retrun code: %d",
                          name, e.return_code)
            except FileNotFoundError as e:
                log.error("%s: %s", name, e)
            finally:
                if (stoptime and starttime and stoptime - starttime
                        > self.logtimethres):
                    log.info("%s execution took %s seconds",
                             name, stoptime - starttime)
                if (iotime and stoptime and iotime - stoptime > self.logtimethres):
                    log.info("%s io took %s seconds",
                             name, stoptime - starttime)


def main():
    recorder = Recorder(outdir=save_path)
    recorder.record(recordables)


if __name__ == '__main__':
    main()
