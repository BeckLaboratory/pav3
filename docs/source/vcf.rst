VCF output
==========

PAV writes one VCF per sample alongside its Parquet :doc:`variant tables <schemas>`. The header is
built by :func:`pav3.vcf.get_headers` and the records by :func:`pav3.vcf.reformat_vcf_table`.

PAV targets VCF 4.5. Insertions, deletions, and SNVs are written with sequence-resolved ``REF`` and
``ALT``; inversions, duplications, and complex variants are written with symbolic alleles.

ALT alleles
-----------

.. vcf-fields:: pav3.vcf.ALT_FIELDS

INFO fields
-----------

.. vcf-fields:: pav3.vcf.INFO_FIELDS

FILTER fields
-------------

.. vcf-fields:: pav3.vcf.FILTER_FIELDS

FORMAT fields
-------------

.. vcf-fields:: pav3.vcf.FORMAT_FIELDS

Lengths and coordinates
-----------------------

``INFO/SVLEN`` and the reference span implied by ``POS`` and ``INFO/END`` measure different
quantities, and for some variant types they do not agree. Neither is derived from the other.

``INFO/END`` and ``INFO/SEQ`` are written only for symbolic ALT alleles. A sequence-resolved ``REF``
and ``ALT`` pair already gives both the span and the sequence, and VCF 4.4 and later treat ``END`` as
a computed field. The table below describes the quantities wherever they are written.

.. list-table::
   :header-rows: 1
   :widths: 10 30 30 30

   * - Type
     - ``SVLEN``
     - ``END`` minus ``POS``
     - Length of ``INFO/SEQ``
   * - ``SNV``
     - Not written
     - Not written
     - Not written
   * - ``INS``
     - Inserted bases (query)
     - Not written; spans no reference bases
     - Inserted bases (query)
   * - ``DEL``
     - Deleted bases (reference)
     - Same as ``SVLEN``
     - Deleted bases (reference)
   * - ``INV``
     - Inverted bases (reference)
     - Same as ``SVLEN``
     - Length of the inverted allele (query)
   * - ``DUP``
     - Duplicated bases
     - Same as ``SVLEN``
     - Duplicated bases (query)
   * - ``CPX``
     - Inserted bases (query)
     - Reference breakpoint span
     - Inserted bases (query)

Two cases deserve attention:

``INV``
   ``SVLEN`` and ``END`` both give the inverted reference span, but ``INFO/SEQ`` holds the inverted
   allele taken from the query. An inversion carrying indels has an allele whose length differs from
   the reference span, so ``len(SEQ)`` is not a size field.

``CPX``
   A complex variant has a reference breakpoint part and an inserted sequence part, and the two are
   reported separately. ``POS`` and ``END`` span the reference bases deleted or duplicated at the
   complex site, and ``SVLEN`` counts the ectopic sequence inserted there. A complex variant that
   deletes 500 bp and inserts 200 bp is written with ``SVLEN=200`` and an ``END`` 500 bases past
   ``POS``. The two must not be summed or compared.

   ``SVLEN`` counts tandem-duplicated reference bases as inserted, because when the flanking
   alignments overlap on the reference those bases appear as additional query sequence. Bases
   deleted at the complex site are not counted, as they contribute nothing to the query.

   A complex variant whose reference breakpoint part is a plain insertion site has ``END`` one base
   past ``POS``, which overstates its zero-length reference span by one base.

Tandem duplications
   Bases duplicated at a variant site fall inside both flanking alignments, so they are absent from
   the query region between them. They are counted by ``SVLEN`` and are written into ``ALT`` and
   ``INFO/SEQ`` by reading the reference over the duplicated region, the same region reported by
   ``INFO/DUP``. Sequence divergence between the two copies is therefore not represented. An
   insertion whose inserted bases are entirely a tandem duplication would otherwise carry no
   sequence at all.

Breakpoint homology
-------------------

``INFO/HOMLEN`` reports the total length of breakpoint homology, summing the upstream and downstream
sides. It is taken from the reference breakpoints for ``DEL`` and ``DUP`` and from the query
breakpoints for ``INS`` and ``CPX``, whichever describes that variant type better. SNVs have no
breakpoint homology, and PAV does not refine inversion breakpoints finely enough to report it.

The variant tables keep the two sides separately in the ``hom_ref`` and ``hom_qry`` columns, so the
upstream and downstream split is available there. PAV does not record the homologous sequence, so
``INFO/HOMSEQ`` is not written.

Structural variant claims
-------------------------

VCF requires ``INFO/SVCLAIM`` on symbolic ``DEL`` and ``DUP`` alleles to say whether a call rests on
a change in DNA abundance, a non-reference adjacency, or both. PAV writes ``J`` for adjacency.

PAV calls deletions and duplications from breakpoints between aligned segments and never measures
DNA abundance, so it cannot tell whether a deleted segment is absent from the genome or relocated
elsewhere. Claiming ``DJ`` would assert a copy-number observation PAV has not made. ``INS`` and
``INV`` carry an implicit adjacency claim and are left unannotated, as are sequence-resolved records,
where ``REF`` and ``ALT`` leave nothing to disambiguate.

Deviations from the specification
---------------------------------

``<CPX>``
   VCF 4.5 requires the first-level type of a symbolic structural variant allele to be one of
   ``DEL``, ``INS``, ``DUP``, ``INV``, or ``CNV``, and allows implementations to define subtypes of
   those. PAV writes ``<CPX>`` for complex variants, which falls outside that set. No conformant
   allele describes a complex rearrangement, and the same extension is in wide use.

``INFO/SVLEN`` on ``<CPX>``
   The specification defines ``SVLEN`` only for ``INS``, ``DEL``, ``DUP``, ``INV``, and ``CNV``
   alleles and directs writers to use the missing value for all others. PAV writes ``SVLEN`` on
   complex variants because the inserted length is more useful to consumers than a missing value.

``INFO/SUBTYPE`` for duplications
   Tandem and inverted duplications are written as the ``<DUP:TANDEM>`` and ``<DUP:INVDUP>`` alleles
   rather than ``INFO/SUBTYPE``, following the VCF convention of writing a symbolic allele as a
   colon-separated list of types and subtypes. ``INVDUP`` is a PAV-defined subtype; the
   specification allows implementations to define their own. Other variant types keep their subtype
   in ``INFO/SUBTYPE``.

``INFO/SVTYPE``
   Not written. VCF 4.5 deprecates it as redundant with ``ALT``. The variant type is carried by the
   symbolic allele for ``INV``, ``DUP``, and ``CPX``, and by the ``REF`` and ``ALT`` lengths for
   ``INS``, ``DEL``, and ``SNV``.

Reserved characters
-------------------

VCF reserves the comma as the delimiter between values within an INFO field, and forbids literal
semicolons and equals signs. PAV writes no INFO value containing any of these. Fields holding more
than one value separate them with ``/``: ``INFO/DUP`` writes its regions that way.

``INFO/CPX_REF_TRACE`` and ``INFO/CPX_QRY_TRACE`` report segment lengths without thousands
separators. The ``trace_ref`` and ``trace_qry`` columns of the variant tables keep the separators
for readability, so the trace strings in the VCF and in the Parquet tables differ.
