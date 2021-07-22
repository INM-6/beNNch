import sys

from analysis_helper import shell, shell_return, update_catalogue
from plot_helper import plot


model = str(sys.argv[1])
jube_bench_path = str(sys.argv[2])
jube_id = str(sys.argv[3])
# overwrite = str(sys.argv[4])
catalogue_path = './catalogue.yaml'
catalogue_path = str(sys.argv[4])

uuidgen_hash = shell_return('uuidgen')
shell(
    f'module load JUBE; jube analyse {jube_bench_path} --id {jube_id};'
    + f' jube result {jube_bench_path} --id {jube_id} > '
    + f'{jube_bench_path}/{jube_id.zfill(6)}/{uuidgen_hash}.csv')

update_catalogue(catalogue_path=catalogue_path,
                 uuidgen_hash=uuidgen_hash)

plot(
    model=model,
    timer_hash=uuidgen_hash,
    timer_path=f'{jube_bench_path}/{jube_id.zfill(6)}',
    catalogue_path=catalogue_path
)
