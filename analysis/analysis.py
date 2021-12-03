"""
NEST Benchmarking Framework - Unified execution, collection, analysis and
comparison of neural network simulation benchmarks.
Copyright (C) 2021 Forschungszentrum Juelich GmbH, INM-6

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

SPDX-License-Identifier: GPL-3.0-or-later
"""

import os
import sys

from analysis_helper import shell, shell_return, load, git_annex
from analysis_config import scaling_type, jube_bench_path
from plot_helper import plot

jube_id = str(sys.argv[1])
base_path = os.path.join(jube_bench_path, jube_id.zfill(6))
uuidgen_hash = shell_return('uuidgen')
shell(
    f'module load JUBE; jube analyse {jube_bench_path} --id {jube_id};'
    + f' jube result {jube_bench_path} --id {jube_id} > '
    + f'{base_path}/{uuidgen_hash}.csv')

cpu_info = load(os.path.join(base_path, '000000_bench/work', 'cpu.json'))
job_info = load(os.path.join(base_path, '000000_bench/work', 'job.json'))

git_annex(cpu_info=cpu_info,
          job_info=job_info,
          uuidgen_hash=uuidgen_hash,
          base_path=base_path)

plot(
    scaling_type=scaling_type,
    timer_hash=uuidgen_hash,
    timer_file=f'{jube_bench_path}/{jube_id.zfill(6)}/{uuidgen_hash}.csv',
    save_path=f'{jube_bench_path}/{jube_id.zfill(6)}'
)
