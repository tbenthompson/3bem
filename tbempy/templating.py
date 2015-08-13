import mako.template
from tbempy.setup import get_tbempy_srces, get_tbempy_headers
import os

def process_templated_code():
    templating_map = create_templating_map()
    perform_templating(templating_map)
    create_gitignore(templating_map.values())

def create_templating_map():
    templating_map = dict()
    for f in get_tbempy_srces() + get_tbempy_headers():
        filename_pieces = f.split('.')
        indices = [
            i for i, piece in enumerate(filename_pieces) if piece == 'template'
        ]
        if len(indices) > 0:
            del filename_pieces[indices[0]]
            new_filename = '.'.join(filename_pieces)
            templating_map[f] = new_filename
    return templating_map

def perform_templating(templating_map):
    for in_filename, out_filename in templating_map.iteritems():
        print(
            'Processing template ' +
            in_filename +
            ' to source file ' +
            out_filename
        )
        t = mako.template.Template(filename = in_filename)
        with open(out_filename, 'w') as out_file:
            out_file.write(t.render())

def create_gitignore(out_filenames):
    print('Creating gitignore for templating results')
    gitignore_filepath = os.path.join('cpp', '.gitignore')
    with open(gitignore_filepath, 'w') as gitignore_file:
        lines = [os.path.split(f_path)[-1] + '\n' for f_path in out_filenames]
        lines.append('.gitignore')
        gitignore_file.writelines(lines)
