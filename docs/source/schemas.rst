Table schemas
=============

PAV writes its tables in Parquet format. The columns and types below are read directly from
:mod:`pav3.schema` when this page is built, so they match the version documented here.

Coordinates are 0-based and half-open (BED-like). Query coordinates refer to the original query
sequence even when the alignment reverse-complemented it, which differs from SAM/BAM; this allows a
variant locus to be extracted from the assembly without first translating coordinates.

Alignment tables
----------------

Alignment records, excluding features added later by :mod:`pav3.align.features`.

.. schema-table:: pav3.schema.ALIGN

Variant tables
--------------

Variant calls. See :py:data:`pav3.schema.VARIANT` for notes on coordinates, inversion inner and
outer coordinates, and merge fields.

.. schema-table:: pav3.schema.VARIANT

Casting
-------

Use :func:`pav3.schema.cast` to apply a schema to a Polars ``DataFrame`` or ``LazyFrame``.
