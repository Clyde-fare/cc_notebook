__author__ = 'clyde'

import io
import re
import random
import os
import copy
import numpy as np
import ConfigParser
from IPython.core import display
from IPython import nbformat
from IPython.core.getipython import get_ipython
from .install import enable_notebook

# activates javascript dependencies
enable_notebook()

config = ConfigParser.RawConfigParser()
config.read(os.path.expanduser('~/.cc_notebook.ini'))

#version of IPython - controls notebook reading and writing functions
ipy_ver=3

#python functions to used by cc_notebook.js to enable a smart-log button, actual programs used set in ~/.cc_notebook.ini

def click(fn):
    command = config.get('smart logs', 'click') + ' ' + fn + ' &'
    os.system(command)

def s_click(fn):
    command = config.get('smart logs', 'shift_click') + ' ' + fn + ' &'
    os.system(command)

def c_click(fn):
    command = config.get('smart logs', 'cntrl_click') + ' ' + fn + ' &'
    os.system(command)

def cs_click(fn):
    command = config.get('smart logs', 'cntrl_shift_click') + ' ' + fn + ' &'
    os.system(command)

def view_ipython_jmol(files, width=300, height=300, sync=False, label=False, title=True, vib=0, delta=None,
                      params=None, script=None, **kwargs):
    """views a file with jsmol from within an ipython notebook.

    label (bool) shows atom no's,
    vib (int) shows the xth vibration is shown,
    title (bool/string, +/- list) titles the applets,
    delta (ase.atom/file) colours the molecules bonds by the difference between the bond lengths given in the files objects and those given by this keyword,
    kwargs commands and require subcommands and selection criteria e.g. color=['red', [0,1,2,3,4]]
    params (list) allows different kwarg dicts to be applied to the individual applets"""

    files = copy.deepcopy(files)

    #this means we can pass a single file string or a list of file_strings
    if not isinstance(files, (list, tuple)):
        files = [files]

    #means we can pass Ase molecules instead of files (need to be careful if the ase object geometry
    #and the .log file don't match  e.g. after a geometry optimisation, we opt to look at the log
    #file geometry if it's defined

    for i in range(len(files)):
        try:
            #if we have a Gaussian calculator with a log file use that
            files[i] = files[i].calc.log
        except AttributeError:
            #if we have an ase object write an xyz file and use that
            try:
                #if the ase object has a single-point calculator with a label variable use that as the filename
                fn = files[i].calc.label + '.xyz'
            except AttributeError:
                #otherwise use a temporary name
                fn = 'temp_atoms_{}.xyz'.format(i)

            try:
                files[i].write(fn)
                files[i] = fn
            except AttributeError:
                pass

    notebook_path = os.path.realpath(get_ipython().starting_dir)
#    current_abs_path = os.getcwd()
#    current_rel_path = os.path.relpath(current_abs_path, notebook_path)

    abs_fs = [os.path.realpath(f) for f in files]
    rel_fs = [os.path.relpath(abs_f, notebook_path) for abs_f in abs_fs]

    #we auto label by atom_no
    init_string = 'select all; '
    if label:
        label_size = int(float(width)/40)
        init_string += 'label %[atomNo]; font label {s};'.format(s=label_size)

    #this allows us to pass in keywords that will be applied to the atoms specified (in all applets if there are multiple)
    if kwargs:
        atom_str = ""
        for arg in kwargs:
            atom_str += 'select ' + ','.join(['atomno={no}'.format(no=num+1) for num in kwargs[arg][1]]) + '; {command} {sub_command};'.format(command = arg, sub_command = kwargs[arg][0])
        init_string += atom_str

    #vibrations in a gaussian file start from frame 2 corresponding to the first vibration (frame 1 is the static molecule)
    #so to specify the first frequency (which is the negative frequency in a transition state) we set vib=1 which loads frame 2
    #for more info see http://jmol.sourceforge.net/demo/vibration/
    vib_strs = []
    if vib and not isinstance(vib, (list, tuple)):
        vib_strs = [' vectors on; color vectors yellow; move 10 -20 10 0 0 0 0 0 1; delay 1; vibration on; frame {n};'.format(n=vib+1) for i in range(len(rel_fs))]
    elif vib and len(vib) != len(rel_fs):
        raise RuntimeError('Must supply a vibrational mode for every applet')
    elif vib:
        vib_strs = [' vectors on; color vectors yellow; move 10 -20 10 0 0 0 0 0 1; delay 1; vibration on; frame {n};'.format(n=v+1) for v in vib]
    else:
        vib_strs = ['' for i in range(len(rel_fs))]

    #this includes a title/titles, default is to use the filename
    font_size = int(float(width)/20)
    if title and isinstance(title, (list, tuple)):
        title_strs = ['set echo echoname 50% 90%; font echo {s}; echo {t};'.format(t=title[i], s=font_size) for i in range(len(files))]
    elif title and isinstance(title, basestring):
        title_strs = ['set echo echoname 50% 90%; font echo {s}; echo {t};'.format(t=title, s=font_size) for i in range(len(files))]
    elif title:
        title_strs = ['set echo echoname 50% 90%; font echo {s}; echo {t};'.format(t=files[i].split('.')[0],s=font_size) for i in range(len(files))]
    else:
        title_strs = ['' for i in range(len(files))]

    #this allows us to pass in a list of dictionaries each specifying different keywords for a different applet
    # - i.e. multiple kwargs one for each applet
    atom_strs = []
    if params and len(params) != len(rel_fs):
            raise RuntimeError('Must supply parameters for every applet')
    elif params:
        for p in params:
            atom_str = ""
            for arg in p:
                atom_str += 'select ' + ','.join(['atomno={no}'.format(no=num+1) for num in p[arg][1]]) + '; {command} {sub_command};'.format(command = arg, sub_command = p[arg][0])
            atom_strs.append(atom_str)
    else:
        atom_strs = ['' for i in range(len(rel_fs))]

    #this allows us to pass in a jmol script directly, could abandon params/kwargs in favour of this
    script_strs = []
    if script and not isinstance(script, (list, tuple)):
        script_strs = [script for i in range(len(files))]
    elif script and len(script) != len(rel_fs):
        raise RuntimeError("Must supply a script for every applet")
    elif script:
        script_strs = script
    else:
        script_strs = ['' for i in range(len(rel_fs))]

    if delta:
        #set delta to the log file so that we can check if any of the files being displayed is the structure we are using as a comparison (if so we do not alter its colouring
        try:
            delta = delta.calc.label +'.log'
        except AttributeError:
            pass

        for i, f in enumerate(files):
            if f != delta:
                script_strs[i] += color_by_delta(files[i], delta)


    id_strs = ["id_" + str(int(random.random()*10000000000000)) for f in rel_fs]

    #load the file initialise, turn off logo display and set default synchronisation to slave mode
    script_command_strs = ["load {f}; {init}; set frank off; sync . SLAVE;".format(f=rel_f, init=init_string + vib_strs[i] + title_strs[i] + atom_strs[i] + script_strs[i]) for i, rel_f in enumerate(rel_fs)]

    html_strs = ["""<div style=float:left><div id='applet_div_{id}'></div></div>\n""".format(id=id_str) for id_str in id_strs]

    html_str = "".join(html_strs)

    js_init_str = """//initialisation
    var Info = {{
        color: "#FFFFFF",
        width: {w},
        height: {h},
        serverURL: "/nbextensions/jsmol/jsmol.php ",
        use: "HTML5",
        j2sPath: "/nbextensions/jsmol/j2s",
        console: "jmolApplet0_infodiv"
    }}
    //prevent applet from immediately generating html - allows us to place the applet in the div we wish
    Jmol.setDocument(0)""".format(w=width, h=height)

    js_applet_strs = ["""
    //start applet
    applet_{id} = Jmol.getApplet("applet_{id}", Info)
    //insert applet html into applet_div_{id}
    $("#applet_div_{id}").html(Jmol.getAppletHtml(applet_{id}))
    //execute scripts
    Jmol.script(applet_{id},"{c}")""".format(id=id_str, c=script_command_str) for id_str, script_command_str in zip(id_strs, script_command_strs)]

    js_applet_str = "".join(js_applet_strs)

    #so slow in jsmol it's almost unuseable
    if sync:
        js_script_str= '\nJmol.script(applet_{id}, "sync . ON; sync * set syncMouse off;set syncScript off")'.format(id=id_strs[0])
    else:
        js_script_str = ''

    jsmol_str = html_str + '<script type="text/Javascript">' + js_init_str + js_applet_str + js_script_str + '</script>'

    return display.HTML(jsmol_str)


#select all; connect delete; select atomno=500, atomno=199; connect single; color bonds green;select atomno=500, atomno=198; connect single; color bonds red;
def color_by_delta(atoms1, atoms2):
    """returns the jmol script required to colour the bonds in atoms1 according to the difference in bond length between atoms1 and atoms2"""
    from ase.io import read
    from ase_extensions import ase_utils

    if isinstance(atoms1, basestring):
        atoms1 = read(atoms1)
    if isinstance(atoms2, basestring):
        atoms2 = read(atoms2)

    #bond_dist_delta returns two lists of the same length, the first is a list of bond indices, the second is a list of the changes in the bond length corresponding to that index
    bond_inds_delta = ase_utils.bond_dist_delta(atoms1, atoms2)

    #creates a list of bonds (bond_ind, bond_delta)
    z_bond_inds_delta = zip(*bond_inds_delta)

    #separates out bonds which have grown vs. those which have shrunk
    try:
        p_bond_inds, p_bond_delta = zip(*[e for e in z_bond_inds_delta if e[1]>0])
    except ValueError:
        p_bond_inds, p_bond_delta = [],[]

    try:
        n_bond_inds, n_bond_delta = zip(*[e for e in z_bond_inds_delta if e[1]<=0])
    except ValueError:
        n_bond_inds, n_bond_delta = [],[]

    #zipping and unzipping turns the original numpy arrays into tuples so we turn them back
    p_bond_delta, n_bond_delta = np.array(p_bond_delta), np.array(n_bond_delta)

    max1 = max(p_bond_delta) if any(p_bond_delta) else 0
    max2 = abs(min(n_bond_delta)) if any(n_bond_delta) else 0

    #scale bond_delta, strange that we have to convert to a float for the division - seems to be a numpy bug
    if max1 or max2:
        p_bond_delta /= float(max([max1,max2]))
        n_bond_delta /= float(max([max1,max2]))

    delete_bonds_str = 'select all; connect delete;'
    #jmol assumes that rgb values either go from 0. to 1. or from 1.1 (not 1.0) so we use the scaled bond_delta values
    color_p_bonds_str = ''.join(["select atomno={atom1_ind}, atomno={atom2_ind}; connect single; color bonds [{r} 0 0];".format(atom1_ind = p_bond_inds[i][0]+1, atom2_ind = p_bond_inds[i][1]+1, r=p_bond_delta[i]) for i in range(len(p_bond_delta))])
    color_n_bonds_str = ''.join(["select atomno={atom1_ind}, atomno={atom2_ind}; connect single; color bonds [0 0 {b}];".format(atom1_ind = n_bond_inds[i][0]+1, atom2_ind = n_bond_inds[i][1]+1, b=abs(n_bond_delta[i])) for i in range(len(n_bond_delta))])

    script_str = delete_bonds_str + color_p_bonds_str + color_n_bonds_str

    return script_str

def color_by_curvature(atoms, colorise=None):
    """Returns a script to color a jsmol molecule by the dihedral of an atom with it's three neighboring atoms"""
    from ase_extensions import ase_utils

    dihedral_data = ase_utils.sp2_dihedrals(atoms)
    ignored_atoms = [i for i in range(len(atoms)) if i not in zip(*dihedral_data)[0]]
    ignored_atom_str = "select " +  ','.join(['atomno={no}'.format(no=ind+1) for ind in ignored_atoms]) + "; color atoms [255 255 255]; "
    abs_dihedral_f = lambda a: abs(a) if abs(a) < 90 else abs(abs(a)-180)
    abs_dihedrals = [abs_dihedral_f(d) for i, d in dihedral_data]

    if not colorise:
        colorise = lambda abs_d: abs_d/max(abs_dihedrals)

    color_atoms_str = ''.join(["select atomno={no}; color atoms [{r} 0 0];".format(no=ind+1, r=colorise(abs_dihedral_f(dihedral))) for (ind, dihedral) in dihedral_data])
    script_str = ignored_atom_str + color_atoms_str
    return script_str

def color_bonds_by_mag(set_bond_inds, mags, inv=False, colorise=None):
    """Colors bonds specified by set_bond_inds by values in mags"""

    mags = np.array(mags)
    if not colorise:
        colorise = lambda m: abs(int(inv)-((m-mags.min())/(mags.max()-mags.min()))**3)

    bond_mags = [list(b) + [i] for b,i in zip(set_bond_inds, mags)]
    color_bonds_str = ''.join(["select atomno={atom1_ind}, atomno={atom2_ind}; connect single; color bonds [{r} 0 0];".format(atom1_ind = b1+1, atom2_ind = b2+1, r=colorise(m)) for b1,b2,m in bond_mags])

    return  color_bonds_str

def view_jmol(atoms, *args, **kwargs):
    atoms.write('temp_atoms.xyz')
    return view_ipython_jmol('temp_atoms.xyz', *args, **kwargs)


def gen_good_movie(fname, atoms, frame_rate=1, bbox=None, pov=False):
    """Generates movie from a list of atoms by generating a collection of pngs and then using ffmpeg to stitch the together

    Requires ffmpeg to be installed on the system

    Optionally uses pov ray for much higher quality (note this makes generating the movie much slower)"""

    from ase.io import write

    for i,mol in enumerate(atoms):
        if pov:
            write('temp_frame_{f:05d}.pov'.format(f=i),mol, bbox=bbox,run_povray=True,display=False)
        else:
            write('temp_frame_{f:05d}.png'.format(f=i),mol, bbox=bbox)

    stitching_png = 'ffmpeg -y -framerate {fr} -i temp_frame_%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p {fn}.mp4; rm temp_frame_*'.format(fn=fname,fr=frame_rate)
    os.system(stitching_png)


def gen_movie(title, list_atoms):
    from ase.io import xyz

    with open(title+'.xyz', 'w') as movie_f:
        xyz.write_xyz(movie_f, list_atoms)

def movie_from_frames(mol, use_fchk=False):
    import ase

    movie_fn = mol.calc.label +'_movie.xyz'
    if not os.path.exists(movie_fn):
        atom_nos = mol.calc.max_data['atomnos']
        coord_frames = mol.calc.max_data['atomcoords']
        atoms_list = [ase.Atoms(numbers=atom_nos, positions=coord_frames[i]) for i in range(len(coord_frames))]
        gen_movie(movie_fn.rstrip('.xyz'), atoms_list)

    view_ipython_jmol(movie_fn)


def mols_to_html(list_mols, data_func=None, sort=None, colour=None):
    """Takes a list of ASE Atoms objects associated with Gaussian calculations and returns an HTML table with a smart log button for each calculation"""

    if any([not mol.calc or 'Gaussian' not in mol.calc.get_name() for mol in list_mols]):
        raise RuntimeError('currently only implemented for Gaussian')

    if not colour:
        colour = 'white'

    if sort and isinstance(sort, bool):
        list_mols.sort(key=lambda e: e.calc.label)
    elif sort:
        list_mols.sort(key=sort)

    #defines 'view_function'
    #html_script = jscript

    names = [mol.calc.label for mol in list_mols]
    log_links = ['<input type="button" value="Smart Log" id="{n}" onclick="IPython.cc_notebook.view_function(event, this.id)" />'.format(n=mol.calc.label) for mol in list_mols]
    status = [mol.calc.status for mol in list_mols]
    notes = [mol.calc.notes for mol in list_mols]

    #used to extract data from the list of molecules and display it in the subsequent table
    if data_func:
        ad_data = []
        for mol in list_mols:
            try:
                ad_data.append(data_func(mol))
            except IOError:
                ad_data.append(None)
    else:
        ad_data = [None for i in range(len(list_mols))]

    init = """
    <style>
    table
    {{
        border-collapse:collapse;
    }}
    td
    {{
        padding:15px;
    }}
    </style>
    <body>
    <table bgcolor="{c}">
    <col/>""".format(c=colour)

    switch_dict = {0: names, 1: log_links, 2: status, 3: notes, 4: ad_data}
    lst_table_eles = []

    table_eles = ''
    for i in range(len(list_mols)):
        temp_data = [names[i], log_links[i], status[i], notes[i], ad_data[i]]
        row_data = [data for n, data in enumerate(temp_data) if any(switch_dict[n]) or n==4 and data_func]
        row_element = '<tr>' + ''.join(['<td>{d}</td>'.format(d=datum) for datum in row_data]) + '</tr>'
        lst_table_eles.append(row_element)
        table_eles = ''.join(lst_table_eles)

    if not len(list_mols):
        table_eles = ''

    final = """
    </table>
    </body>"""

    return display.HTML(init + table_eles + final)


def unwind(lst_of_lsts):
    return [e for l in lst_of_lsts for e in l]


def copy_cells(from_notebook='', from_cells=(0, 0), to_notebook='', at_cell=0):

    if '.ipynb' not in from_notebook:
        from_notebook += '.ipynb'
    if '.ipynb' not in to_notebook:
        to_notebook += '.ipynb'

    with open(from_notebook) as nb1_f:
        nb1_raw = nb1_f.read()

    with open(to_notebook) as nb2_f:
        nb2_raw = nb2_f.read()

    nb1 = nbformat.reads(nb1_raw, as_version=ipy_ver)
    nb2 = nbformat.reads(nb2_raw, as_version=ipy_ver)

    start_id, end_id = from_cells
    copied_cells = nb1['worksheets'][0]['cells'][start_id:end_id]

    active_cells = nb2['worksheets'][0]['cells']
    nb2['worksheets'][0]['cells'] = active_cells[:at_cell] + copied_cells + active_cells[at_cell:]

    nb2_modified_raw = nbformat.writes(nb2, as_version=ipy_ver)

    with open(to_notebook, 'w') as nb2_f:
        nb2_f.write(nb2_modified_raw)


def extract_linked_nb(mkdn_cell):
    """Extracts url links to other notebooks from a markdown cell"""
    nb_link_extractor = re.compile("(?<=\]\().+ipynb")
    linked_nbs = nb_link_extractor.findall(mkdn_cell['source'])
    return linked_nbs


def inherit_upto_cell(nb_name, cell_id=-1, silent=False):
    """inherits from another notebook up to cell number cell_id"""

    ip = get_ipython()

    if '.ipynb' not in nb_name:
        nb_name += '.ipynb'

    with io.open(nb_name) as f:
        ancestor_nb = nbformat.read(f, as_version=ipy_ver)

    code_cells = [cell for cell in ancestor_nb.worksheets[0].cells if cell.cell_type == 'code']

    if cell_id == -1:
        cell_id = len(code_cells)
        print('Inheriting entire notebook')

    ancestor_code_cells = [cell for i, cell in enumerate(code_cells) if i <= cell_id]

    for cell in ancestor_code_cells:
        ip.run_cell(cell.input, silent=silent)

    return


def inherit_from(nb_name, uptolink='', silent=False):
    """Inherits code from notebook 'nb_name' up to the point where a link to 'uptolink' is found in a markdown cell of nb_name"""

    if '.ipynb' not in nb_name:
        nb_name += '.ipynb'

    if '.ipynb' not in uptolink:
        uptolink += '.ipynb'

    with io.open(nb_name) as f:
        ancestor_nb = nbformat.read(f, as_version=ipy_ver)

    cell_id = -1

    for cell in ancestor_nb.worksheets[0].cells:
        if cell.cell_type == 'code':
            cell_id+=1
        if cell.cell_type == 'markdown' and uptolink in extract_linked_nb(cell):
            break
    else:
        inherit_upto_cell(nb_name, cell_id=-1, silent=silent)
        return

    inherit_upto_cell(nb_name, cell_id=cell_id, silent=silent)
