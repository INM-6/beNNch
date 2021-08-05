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
os.system('mkdir -p results/plots')
for csv_file in csv_files:
    timer_hash = csv_file.split('/')[-1].split('.')[0]
    plot(model=model,
         timer_hash=timer_hash,
         timer_path='results',
         save_path='results/plots'
         )

plot_files = glob.glob('results/plots/*.png')


def display_attributes(file, attributes):
    display_list = '<ul>\n'

    for attribute in attributes:
        value = os.popen(
            f'git annex metadata {file} '
            + f'--get {attribute}').read().strip()
        display_list += f'  <li>{attribute}: {value}</li>\n'

    display_list += '</ul>'

    display(Image(filename=file))
    display(HTML(display_list))
    return


def make_notebook(outPath: str):
    nb = nbf.new_notebook()
    cells = []
    codes = {
        'skip': ["from slideshow.slideshow import display_attributes"],
        'slide': [f"display_attributes('{plot_files[0]}',"
                  + f"{attributes_to_display})"],
        'subslide': []
    }

    for file in plot_files[1:]:
        codes['subslide'].append(f"display_attributes('{file}',"
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
