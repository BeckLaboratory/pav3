Command-line interface
======================

PAV is normally run through its command-line interface. The options below are generated from the
PAV argument parser, so they always match the version of PAV documented here.

For a guide to setting up the configuration files that ``pav3 batch`` reads, see the
`README <https://github.com/BeckLaboratory/pav3/blob/main/README.md>`_ in the PAV repository.

.. argparse::
   :module: pav3.cli
   :func: build_parser
   :prog: pav3
