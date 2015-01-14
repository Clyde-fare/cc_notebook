//set cc_notebook javascript functions
IPython.cc_notebook = {}
IPython.cc_notebook.view_function = function (event, file_n) {
    if (event.ctrlKey && event.shiftKey){
        IPython.notebook.kernel.execute("from cc_notebook import pygview; pygview('" + file_n + ".log')");
    }
    else if (event.ctrlKey){
        IPython.notebook.kernel.execute("from cc_notebook import pygausssum; pygausssum('" + file_n + ".log')");
    }
    else if (event.shiftKey){
        IPython.notebook.kernel.execute("from cc_notebook import pyvogadro; pyvogadro('" + file_n + ".log')");
    }
    else {
        IPython.notebook.kernel.execute("from cc_notebook import pyvim; pyvim('" + file_n + ".log')");
        IPython.notebook.kernel.execute("print(file_n)");
    }
}
