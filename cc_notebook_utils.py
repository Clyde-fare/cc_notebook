__author__ = 'clyde'

from IPython.core import display
import IPython.nbformat.current as nb_current
import tempfile
import os
import copy
import ASE_utils
import numpy as np
import pandas


JMOL_PATH = '/home/clyde/Dropbox/Project Stuff/Notebooks/jmol'
NOTEBOOK_PATH = '/home/clyde/Dropbox/Project Stuff/Notebooks/'

#this function opens the linked file in gaussSum/gaussview/avogadro/vi depending on the command key held down
##depreciated since IPython 2.0 - javascript utility functions are now defined in custom.js
#jscript = """
           # <script type="text/javascript">
           # function view_function (event, file_n) {
           #     if (event.ctrlKey && event.shiftKey){
           #         IPython.notebook.kernel.execute("cc_notebook_utils.pygview('" + file_n + ".log')");
           #     }
           #     else if (event.ctrlKey){
           #         IPython.notebook.kernel.execute("cc_notebook_utils.pygausssum('" + file_n + ".log')");
           #     }
           #     else if (event.shiftKey){
           #         IPython.notebook.kernel.execute("cc_notebook_utils.pyvogadro('" + file_n + ".log')");
           #     }
           #     else {
           #         IPython.notebook.kernel.execute("cc_notebook_utils.pyvim('" + file_n + ".log')");
           #         IPython.notebook.kernel.execute("print(file_n)");
           #     }
           # }
           # </script>"""

#from gaussian_job_manager import extract_status
#def mols_to_html(list_mols, colour=None):
#    import os
#
#    if any([not mol.calc or mol.calc.get_name() != 'Gaussian' for mol in list_mols]):
#        raise RuntimeError('currently only implemented for Gaussian')
#
#    if not colour:
#        colour = 'white'
#
#    #defines 'view_function'
#    html_script = jscript
#
#    names = [mol.calc.label for mol in list_mols]
#    log_links = ['<input type="button" value="Smart Log" id="{n}" onclick="view_function(event, this.id)" />'.format(n=mol.calc.label) for mol in list_mols]
#    com_links = ['<a href = "files/{p}/{n}.com" target = "_blank">com</a>'.format(p=os.path.relpath(os.getcwd(),'/home/clyde/Dropbox/Project Stuff/Notebooks/'), n=mol.calc.label) for mol in list_mols]
#    status = [extract_status(name+'.log') for name in names]
#
#
#    init =  """
#    <style>
#    table
#    {{
#        border-collapse:collapse;
#    }}
#    td
#    {{
#        padding:15px;
#    }}
#    </style>
#    <body>
#    <table bgcolor="{c}">
#    <col/>""".format(c=colour)
#
#    table_eles = "".join(['<tr><td>{n}</td><td>{c}</td><td>{l}</td><td>{s}</td></tr>'.format(n=names[i], c=com_links[i], l=log_links[i], s=status[i]) for i in range(len(names))])
#
#    final = """
#    </table>
#    </body>"""
#
#    return display.HTML(html_script + init + table_eles + final)

def pyvim(fn):
    os.system('gvim {f} &'.format(f=fn))
def pyvogadro(fn):
    os.system('avogadro {f} &'.format(f=fn))
def pygausssum(fn):
    os.system('GaussSum.py {f} &'.format(f=fn))
def pygview(fn):
    if '.log' in fn:
        fn.replace('.log', '.chk')

    scratch_dir = ASE_utils.get_equiv_scratch_dir()
    os.system("ssh -Y cjf05@login.cx1.hpc.ic.ac.uk 'source /etc/profile.d/module.sh;module load gaussview;gview {d}{f}' &".format(d=scratch_dir,f=fn))


#todo currently only works from Notebooks directory rather than inside subdirectories
#can't call IPython.core.display.HTML from a function but can return the object (or indeed any object with a _repr_html_ method that returns an html string
#this has the downside that we can't delete the temporary file that we are using to hold jmol_basic_str, so instead we simply delete the temporary files left over
#by previous calls at the beginning
def view_ipython_jmol(files, width=400, height=300, label=False, title=True, vib=0, delta=None, params=None, script=None, **kwargs):
    """views a file with jmol from within an ipython notebook, label (bool) shows atom no's, vib (int) shows the xth vibration is shown, title (bool/string, +/- list) titles the applets, delta (ase.atom/file) colours the molecules bonds by the difference between the bond lengths given in the files objects and those given by this keyword, kwargs commands and require subcommands and selection criteria e.g. color=['red', [0,1,2,3,4]], params (list) allows different kwarg dicts to be applied to the individual applets"""
    files = copy.deepcopy(files)

    jmol_path = JMOL_PATH
    notebook_path = NOTEBOOK_PATH

    #this means we can pass a single file string or a list of file_strings
    if not isinstance(files, (list, tuple)):
        files = [files]

    #means we can pass Ase molecules instead of files
    for i in range(len(files)):
        try:
            files[i] = files[i].calc.label +'.log'
        except AttributeError:
            pass

    abs_fs = [os.path.abspath(f) for f in files]
    rel_jmol_path = os.path.relpath(jmol_path, notebook_path)
    rel_fs = [os.path.relpath(abs_f, jmol_path) for abs_f in abs_fs]

    #need another way of doing this because if we immediately delete we cannot have multiple molecules open at a time
  # delete previous temp_files
  # for fl in glob.glob('{pth}/*ase-*.html'.format(pth=jmol_path)):
  #     os.remove(fl)

    #we auto label by atom_no and set background to white
    init_string = 'select all; '
    if label:
        init_string += 'label %[atomNo]; '
    init_string += 'background [1 1 1];'

    #this allows us to pass in keywords that will be applied to the atoms specified (in all applets if there are multiple)
    if kwargs:
        atom_str = ""
        for arg in kwargs:
            atom_str += 'select ' + ','.join(['atomno={no}'.format(no=num+1) for num in kwargs[arg][1]]) + '; {command} {sub_command};'.format(command = arg, sub_command = kwargs[arg][0])
        init_string += atom_str

    #vibrations in a gaussian file start from frame 2 corresponding to the first vibration (frame 1 is the static molecule)
    #so to specify the first frequency (which is the negative frequency in a transition state) we set vib=1 which loads frame 2, for more info see http://jmol.sourceforge.net/demo/vibration/
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
    font_size = int(width/50)
    title_strs = []
    if title and isinstance(title, (list, tuple)):
        title_strs = ['set echo echoname 50% 90%; font echo {s}; echo {t};'.format(t=title[i],s=font_size) for i in range(len(files))]
    elif title and isinstance(title, basestring):
        title_strs = ['set echo echoname 50% 90%; font echo {s}; echo {t};'.format(t=title,s=font_size) for i in range(len(files))]
    elif title:
        title_strs = ['set echo echoname 50% 90%; font echo {s}; echo {t};'.format(t=files[i].split('.')[0],s=font_size) for i in range(len(files))]
    else:
        title_strs = ['' for i in range(len(files))]

    #this allows us to pass in a list of dictionaries each specifying different keywords for a different applet - i.e. multiple kwargs one for each applet
    atom_strs = []
    if params and len(params) != len(rel_fs):
            raise RuntimeError('Must supply parameters for every applet')
    elif params:
        for p in params:
            atom_str = ""
            for arg in p:
                atom_str+= 'select ' + ','.join(['atomno={no}'.format(no=num+1) for num in p[arg][1]]) + '; {command} {sub_command};'.format(command = arg, sub_command = p[arg][0])
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

    applet_strings = ['jmolApplet({h}, "load {xyz}; {init}; sync . on");'.format(h=height -50, init=init_string + vib_strs[i] + title_strs[i] + atom_strs[i] + script_strs[i], xyz= rel_f) for i,rel_f in enumerate(rel_fs)]
    applet_string = '\n'.join(applet_strings)
#    applet_string = 'jmolApplet({h}, load FILES '.format(h=height -50) + ' '.join(['"{xyz}"'.format(xyz= rel_f) for rel_f in rel_fs]) + '; {init}; sync . on");'.format(init=init_string)

    jmol_basic_str="""<head><script type="text/javascript" src="Jmol.js"></script></head>
    <body>
    <script type="text/javascript">
      jmolInitialize("");
      jmolSetCallback("UseCommandThread","true");
      {app_str}
    </script>
    </body>""".format(app_str = applet_string)

    fd, jmol_html = tempfile.mkstemp('.' + 'html', 'ase-', dir=jmol_path)
    fd = os.fdopen(fd, 'w')
    fd.write(jmol_basic_str)
    fd.close()

    max_h = int(1680 / width)
    horizontal_stack = max_h if len(files) >=max_h else len(files)
    vertical_stack = 1+int(len(files)/max_h) if len(files) % max_h else int(len(files)/max_h)

    jmol_str = """
    <iframe
      width="{w}"
      height="{h}"
      src=/files/{r_jmol_path}/{src}
      frameborder="0"
      allowfullscreen
    ></iframe>""".format(w=50 + (width-50) * horizontal_stack, h=110 + (height-110) * vertical_stack, r_jmol_path = rel_jmol_path, src=os.path.basename(jmol_html))

    return display.HTML(jmol_str)

#select all; connect delete; select atomno=500, atomno=199; connect single; color bonds green;select atomno=500, atomno=198; connect single; color bonds red;
def color_by_delta(atoms1, atoms2):
    """returns the jmol script required to colour the bonds in atoms1 according to the difference in bond length between atoms1 and atoms2"""
    import ase
    if isinstance(atoms1, basestring):
        atoms1 = ase.io.read(atoms1)
    if isinstance(atoms2, basestring):
        atoms2 = ase.io.read(atoms2)

    #bond_dist_delta returns two lists of the same length, the first is a list of bond indices, the second is a list of the changes in the bond length corresponding to that index
    bond_inds_delta = ASE_utils.bond_dist_delta(atoms1, atoms2)

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


def view_jmol(atoms, width=400, height=300, label=True):
    atoms.write('temp_atoms.xyz')
    return view_ipython_jmol('temp_atoms.xyz',width, height, label)

def gen_movie(title, list_atoms):
    from ase.io import xyz

    with open(title+'.xyz', 'w') as movie_f:
        xyz.write_xyz(movie_f, list_atoms)


#views the difference between two molecules, better to use view_ipython_jmol using the delta keyword
def view_delta(atoms1, atoms2, **kwargs):
    #bond_dist_delta returns two lists of the same length, the first is a list of bond indices, the second is a list of the changes in the bond length corresponding to that index
    bond_inds_delta = ASE_utils.bond_dist_delta(atoms1, atoms2)

    #creates a list of bonds (bond_ind, bond_delta)
    z_bond_inds_delta = zip(*bond_inds_delta)

    #separates out bonds which have grown vs. those which have shrunk
    try:
        p_bond_inds, p_bond_delta = zip(*[e for e in z_bond_inds_delta if e[1]>0])
    except ValueError:
        p_bond_inds, p_bond_delta = [], []

    try:
        n_bond_inds, n_bond_delta = zip(*[e for e in z_bond_inds_delta if e[1]<0])
    except ValueError:
        n_bond_inds, n_bond_delta = [],[]

    #zipping and unzipping turns the original numpy arrays into tuples so we turn them back
    p_bond_delta = np.array(p_bond_delta)
    n_bond_delta = np.array(n_bond_delta)

    #scale bond_delta, strange that we have to convert to a float for the division - seems to be a numpy bug
    max1 = max(p_bond_delta)
    max2 = abs(min(n_bond_delta))
    p_bond_delta /= float(max([max1,max2]))
    n_bond_delta /= float(max([max1,max2]))

    delete_bonds_str = 'select all; connect delete;'
    #jmol assumes that rgb values either go from 0. to 1. or from 1.1 (not 1.0) so we use the scaled bond_delta values
    color_p_bonds_str = ''.join(["select atomno={atom1_ind}, atomno={atom2_ind}; connect single; color bonds [{r} 0 0];".format(atom1_ind = p_bond_inds[i][0]+1, atom2_ind = p_bond_inds[i][1]+1, r=p_bond_delta[i]) for i in range(len(p_bond_delta))])
    color_n_bonds_str = ''.join(["select atomno={atom1_ind}, atomno={atom2_ind}; connect single; color bonds [0 0 {b}];".format(atom1_ind = n_bond_inds[i][0]+1, atom2_ind = n_bond_inds[i][1]+1, b=abs(n_bond_delta[i])) for i in range(len(n_bond_delta))])

    script_str = delete_bonds_str + color_p_bonds_str + color_n_bonds_str

    return view_ipython_jmol(atoms1, script=script_str, **kwargs)

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
    import os

    if any([not mol.calc or mol.calc.get_name() != 'Gaussian' for mol in list_mols]):
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
    #since IPython 2.0, javascript has to be defined in the profile_{}/static/custom/custom.js which is where IPython.cc_notebook.* is defined
    log_links = ['<input type="button" value="Smart Log" id="{n}" onclick="IPython.cc_notebook.view_function(event, this.id)" />'.format(n=mol.calc.label) for mol in list_mols]
    com_links = ['<a href = "files/{p}/{n}.com" target = "_blank">com</a>'.format(p=os.path.relpath(os.getcwd(),'/home/clyde/Dropbox/Project Stuff/Notebooks/'), n=mol.calc.label) for mol in list_mols]
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

    switch_dict = {0: names, 1: com_links, 2: log_links, 3: status, 4: notes, 5: ad_data}
    lst_table_eles = []

    for i in range(len(list_mols)):
        temp_data = [names[i], com_links[i], log_links[i], status[i], notes[i], ad_data[i]]
        row_data = [data for n, data in enumerate(temp_data) if any(switch_dict[n]) or n==5 and data_func]
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

def copy_cells(from_notebook='', from_cells=(0,0), to_notebook='', at_cell=0):

    if '.ipynb' not in from_notebook:
        from_notebook += '.ipynb'
    if '.ipynb' not in to_notebook:
        to_notebook += '.ipynb'

    with open(from_notebook) as nb1_f:
        nb1_raw = nb1_f.read()

    with open(to_notebook) as nb2_f:
        nb2_raw = nb2_f.read()

    nb1 = nb_current.reads_json(nb1_raw)
    nb2 = nb_current.reads_json(nb2_raw)

    start_id, end_id = from_cells
    copied_cells = nb1['worksheets'][0]['cells'][start_id:end_id]

    active_cells = nb2['worksheets'][0]['cells']
    nb2['worksheets'][0]['cells'] = active_cells[:at_cell] + copied_cells + active_cells[at_cell:]

    nb2_modified_raw = nb_current.writes_json(nb2)

    with open(to_notebook, 'w') as nb2_f:
        nb2_f.write(nb2_modified_raw)

def gen_energy_table(products, reactants, delta=None):
    from ase import atoms
    rxn_Es = []
    names = []

    if delta:
        delt_Es = []

    for i in range(len(products)):
        if isinstance(products[i], atoms):
            main_product = products[i]
            try:
                product_E = products[i].calc.energy_zero
            except AttributeError:
                product_E = float('nan')
        else:
            main_product = products[i][0]
            try:
                product_E = sum([p.calc.energy_zero for p in products[i]])
            except AttributeError:
                product_E = float('nan')

        if isinstance(products[i], atoms):
            try:
                reactant_E = reactants[i].calc.energy_zero
            except AttributeError:
                reactant_E = float('nan')
        else:
            try:
                reactant_E = sum([r.calc.energy_zero for r in reactants[i]])
            except AttributeError:
                reactant_E = float('nan')

        rxn_E = 23.060542301388647 * (product_E - reactant_E)

        if delta:
            delt_E = rxn_E - delta
            delt_Es.append(delt_E)

        rxn_Es.append(rxn_E)
        names.append(main_product.calc.label)

    o_data = pandas.Series(rxn_Es ,names)
    d = {'Rxn Energy': o_data}

    if delta:
        o_dft_delt_data = pandas.Series(delt_Es, names)
        d.update({'Rxn Energy Delta': o_dft_delt_data})
    return pandas.DataFrame(d)

#pymol magic
# see: http://ipython.org/ipython-doc/dev/interactive/reference.html
# and
# http://proj.badc.rl.ac.uk/cedaservices/browser/ipython/IPython/extensions/cythonmagic.py?rev=233076612d2815f8bfe098230e82fb9e3ca3749e
import sys, os
from IPython.core.magic import (Magics, magics_class, cell_magic)

#import pymol
## The class MUST call this class decorator at creation time
#@magics_class
#class MyMagics(Magics):
#    @cell_magic
#    def pyMol(line, cell):
#        """Doc_String"""
#        t_script_contents = ["import {md}".format(md=mod) for mod in sys.modules]
#        t_script_contents += ["import pickle"]
#        t_script_contents += ["{vr} = unpickle {vr_file}".format(vr=var, vr_file=var_file) for var, var_file in list_vs_vfiles]
#        t_script_contents += [cell]
#
#        script_contents = "/n".join(t_script_contents)
#        fname = 'temp_pymol_script' + '.py'
#
#        with open(fname,'w') as f:
#            f.write(script_contents)
#
#        exitcode= os.system("pymol temp_pymol_script.py")
#
#        return 'Success' if not exitcode else 'Fail'
#
## In order to actually use these magics, you must register them with a
## running IPython.  This code must be placed in a file that is loaded once
## IPython is up and running:
##ip = get_ipython()
## You can register the class itself without instantiating it.  IPython will
## call the default constructor on it.
##ip.register_magics(MyMagics
