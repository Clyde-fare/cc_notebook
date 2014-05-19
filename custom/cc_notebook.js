//set cc_notebook javascript functions
IPython.cc_notebook = {}
IPython.cc_notebook.view_function = function (event, file_n) {
    if (event.ctrlKey && event.shiftKey){
        IPython.notebook.kernel.execute("import cc_notebook_utils; cc_notebook_utils.pygview('" + file_n + ".log')");
    }
    else if (event.ctrlKey){
        IPython.notebook.kernel.execute("import cc_notebook_utils; cc_notebook_utils.pygausssum('" + file_n + ".log')");
    }
    else if (event.shiftKey){
        IPython.notebook.kernel.execute("import cc_notebook_utils; cc_notebook_utils.pyvogadro('" + file_n + ".log')");
    }
    else {
        IPython.notebook.kernel.execute("import cc_notebook_utils; cc_notebook_utils.pyvim('" + file_n + ".log')");
        IPython.notebook.kernel.execute("print(file_n)");
    }
}
