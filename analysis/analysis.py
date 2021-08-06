import os
import sys

from analysis_helper import shell, shell_return, load, git_annex
from analysis_config import model, jube_bench_path
from plot_helper import plot

jube_id = str(sys.argv[1])
base_path = os.path.join(jube_bench_path, jube_id.zfill(6))
uuidgen_hash = shell_return('uuidgen')
shell(
    f'module load JUBE; jube analyse {jube_bench_path} --id {jube_id};'
    + f' jube result {jube_bench_path} --id {jube_id} > '
    + f'{base_path}/{uuidgen_hash}.csv')

cpu_info = load(os.path.join(base_path, '000000_bench/work', 'cpu.pkl'))
job_info = load(os.path.join(base_path, '000000_bench/work', 'job.pkl'))

git_annex(cpu_info=cpu_info,
          job_info=job_info,
          uuidgen_hash=uuidgen_hash,
          base_path=base_path)

plot(
    model=model,
    timer_hash=uuidgen_hash,
    timer_file=f'{jube_bench_path}/{jube_id.zfill(6)}/{uuidgen_hash}.csv',
    save_path=f'{jube_bench_path}/{jube_id.zfill(6)}'
)
