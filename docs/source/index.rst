PAV documentation
=================

The Phased Assembly Variant (PAV) caller generates variant calls from genome assemblies.

PAV 3 is a reimplementation of earlier PAV versions. It uses many of the same algorithms and
includes new features for identifying large and complex structural variants. Improved
low-confidence alignment prediction helps to filter variant calls from erroneous alignment records.

PAV 3 is fast. By utilizing modern data science tools, specifically Polars, PAV produces callsets
using a fraction of the memory and CPU resources as earlier versions.

.. warning::

   PAV 3 is under active development and has not had a final release. See
   :doc:`API stability <stability>` before depending on it.

Install
-------

PAV requires Python 3.12 or newer.

.. code-block:: bash

   pip install pav3

Run
---

PAV is run through its built-in command-line interface:

.. code-block:: bash

   pav3 batch

In some Python environments the command must be run through ``python``:

.. code-block:: bash

   python -m pav3 batch

See the :doc:`command-line interface <cli>` for all subcommands and options.

Configuration
-------------

``pav3 batch`` reads a configuration file (``pav.json``) and a table of assemblies
(``assemblies.tsv``). For details about setting these up, along with guidance on choosing a
reference and a description of the output, see the
`README <https://github.com/BeckLaboratory/pav3/blob/main/README.md>`_ in the PAV repository.

This site is the API and command-line reference; the README is the usage guide.

Contributions to the project are welcome. See the
`PAV repository <https://github.com/BeckLaboratory/pav3>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   cli
   schemas
   stability
   changelog
   license
