__author__ = 'clyde'

from IPython.core import display
from IPython.core.getipython import get_ipython
import IPython.nbformat.current as nb_current
import random
import os
import copy
import ASE_utils
import numpy as np
import pandas

#python functions to open vim/avogadro/gaussum/gaussview that are used by cc_notebook.js to enable a smart-log button
def pyvim(fn):
    """Opens default editor"""
    os.system('editor {f} &'.format(f=fn))
def pyvogadro(fn):
    """Opens avogadro"""
    os.system('avogadro {f} &'.format(f=fn))
def pygausssum(fn):
    """Opens Gaussum"""
    os.system('GaussSum.py {f} &'.format(f=fn))
def pygview(fn):
    """Opens gaussview over remote server"""
    if '.log' in fn:
        fn.replace('.log', '.chk')

    scratch_dir = ASE_utils.get_equiv_scratch_dir()
    os.system("ssh -Y cjf05@login.cx1.hpc.ic.ac.uk 'source /etc/profile.d/module.sh;module load gaussview;gview {d}{f}' &".format(d=scratch_dir,f=fn))

#requires jsmol to be present in profile/static/custom
def view_ipython_jmol(files, width=300, height=300, sync=False, label=False, title=True, vib=0, delta=None, params=None, script=None, **kwargs):
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

    #means we can pass Ase molecules instead of files
    for i in range(len(files)):
        try:
            files[i] = files[i].calc.label +'.log'
        except AttributeError:
            pass

    notebook_path = get_ipython().starting_dir
#    current_abs_path = os.getcwd()
#    current_rel_path = os.path.relpath(current_abs_path, notebook_path)

    abs_fs = [os.path.abspath(f) for f in files]
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
        serverURL: "/static/custom/jsmol/jsmol.php ",
        use: "HTML5",
        j2sPath: "/static/custom/jsmol/j2s",
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
    if isinstance(atoms1, basestring):
        atoms1 = read(atoms1)
    if isinstance(atoms2, basestring):
        atoms2 = read(atoms2)

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


def view_jmol(atoms, *args, **kwargs):
    atoms.write('temp_atoms.xyz')
    return view_ipython_jmol('temp_atoms.xyz', *args, **kwargs)

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
    """Takes a list of ASE Atoms objects associated with Gaussian calculations and returns an HTML table with a smart log button for each calculation"""
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

def copy_cells(from_notebook='', from_cells=(0, 0), to_notebook='', at_cell=0):

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
        if isinstance(products[i], atoms.Atoms):
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

        if isinstance(reactants[i], atoms.Atoms):
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

    o_data = pandas.Series(rxn_Es, names)
    d = {'Rxn Energy': o_data}

    if delta:
        o_dft_delt_data = pandas.Series(delt_Es, names)
        d.update({'Rxn Energy Delta': o_dft_delt_data})
    return pandas.DataFrame(d)
