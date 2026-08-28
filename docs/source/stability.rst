API stability
=============

.. warning::

   PAV 3 is a development release. The Python API is not stable and may change without a
   deprecation period until a 3.0.0 release is published.

PAV is primarily a command-line tool. The supported interface is the
:doc:`command-line interface <cli>` and the :doc:`table schemas <schemas>` of its Parquet output;
these are the parts intended to stay stable across releases.

The Python API documented under :doc:`API Reference <autoapi/index>` is published to support
development and inspection of PAV internals. Names beginning with an underscore are private.
Anything else may still be reorganized before 3.0.0, so pin an exact version if you import PAV as
a library.

Pre-release versions carry a suffix: ``.devN`` for development, ``aN`` for alpha, ``bN`` for beta,
and ``rcN`` for release candidates. All of these precede the corresponding final release and should
not be considered production-ready.

Please discuss with the authors prior to publishing results generated with PAV 3; see
`contact information <https://www.jax.org/people/peter-audano>`_. Report problems on the
`issues page <https://github.com/BeckLaboratory/pav3/issues>`_.
