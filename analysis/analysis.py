import os
import sys

from analysis_helper import shell, shell_return, update_catalogue, load
from analysis_config import model, jube_bench_path, catalogue_path
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

update_catalogue(catalogue_path=catalogue_path,
                 uuidgen_hash=uuidgen_hash,
                 cpu_info=cpu_info,
                 job_info=job_info)

try:
    # plotting only works if run goes across nodes or virtual processes
    plot(
        model=model,
        timer_hash=uuidgen_hash,
        timer_path=f'{jube_bench_path}/{jube_id.zfill(6)}',
        catalogue_path=catalogue_path
    )
except ValueError:
    print('plotting only works if run goes across nodes or virtual processes')
