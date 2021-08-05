import glob
import sys
import os

from nbformat import v4 as nbf
import nbformat
from IPython.display import Image
from IPython.display import HTML

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from analysis.plot_helper import plot


model = sys.argv[1]
attributes_to_display = sys.argv[2:]

csv_files = glob.glob('results/*.csv')
timer_hashes = []

os.system('mkdir -p results/plots')
for csv_file in csv_files:
    timer_hash = csv_file.split('/')[-1].split('.')[0]
    timer_hashes.append(timer_hash)
    plot(model=model,
         timer_hash=timer_hash,
         timer_path='results',
         save_path='results/plots'
         )

plot_files = glob.glob('results/plots/*.png')


def display_plot(timer_hash, timer_path, plot_path, attributes):
    display_list = '<ul>\n'
    for attribute in attributes:
        value = os.popen(
            f'git annex metadata {timer_path}/{timer_hash}.csv '
            + f'--get {attribute}').read().strip()
        display_list += f'  <li>{attribute}: {value}</li>\n'
    display_list += '</ul>'

    display(Image(filename=os.path.join(plot_path, timer_hash + '.png')))
    display(HTML(display_list))


def make_notebook(outPath: str):
    nb = nbf.new_notebook()
    cells = []
    codes = {
        'skip': ["from slideshow.slideshow import display_plot"],
        'slide': [f"display_plot('{timer_hashes[0]}', "
                  + "'results', "
                  + "'results/plots', "
                  + f"{attributes_to_display})"],
        'subslide': []
    }

    for timer_hash in timer_hashes[1:]:
        codes['subslide'].append(f"display_plot('{timer_hash}', "
                                 + "'results', "
                                 + "'results/plots', "
                                 + f"{attributes_to_display})")
    for key, code_list in codes.items():
        for code in code_list:
            cells.append(nbf.new_code_cell(code, metadata={
                         "slideshow": {"slide_type": key}}))

    nb['cells'] = cells

    fname = 'slideshow.ipynb'
    with open(os.path.join(outPath, fname), 'w') as _:
        nbformat.write(nb, _)


if __name__ == '__main__':

    make_notebook('./')
    os.system("jupyter nbconvert --inplace --execute slideshow.ipynb")
    os.system("jupyter nbconvert --to slides slideshow.ipynb "
              + "--TemplateExporter.exclude_input=True "
              + "--SlidesExporter.reveal_transition='none'")
    # os.system("rm slideshow.ipynb")
