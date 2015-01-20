//set cc_notebook javascript functions
IPython.cc_notebook = {}
IPython.cc_notebook.view_function = function (event, file_n) {
    if (event.ctrlKey && event.shiftKey){
        IPython.notebook.kernel.execute("from cc_notebook import cs_click; cs_click('" + file_n + ".log')");
    }
    else if (event.ctrlKey){
        IPython.notebook.kernel.execute("from cc_notebook import c_click; c_click('" + file_n + ".log')");
    }
    else if (event.shiftKey){
        IPython.notebook.kernel.execute("from cc_notebook import s_click; s_click('" + file_n + ".log')");
    }
    else {
        IPython.notebook.kernel.execute("from cc_notebook import click; click('" + file_n + ".log')");
    }
}
