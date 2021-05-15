# pyrealm

This forked project is based on pyrealm package devoloped by David Orme and its application on a global scale.

The main code that simulation global GPP is 'rda_bulk_run.py', which use .rda file as data inputs (see example input data 'tmp_0119.rda'). The 'data_preparation.py' reads .rda file in python into arrays so that P model can simulate in bulk. 'pmodel_function.py' is a simple function to calculate vapour pressure deficit using temperature and vapour pressure data.
You can also easily transform the code to calculate GPP at a point or site.

For any questions please contact w.cai17@imperial.ac.uk
