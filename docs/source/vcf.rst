VCF output
==========

PAV writes one VCF per sample alongside its Parquet :doc:`variant tables <schemas>`. The header is
built by :func:`pav3.vcf.get_headers` and the records by :func:`pav3.vcf.reformat_vcf_table`.

PAV targets VCF 4.5. Insertions, deletions, SNVs, and complex variants are written with
sequence-resolved ``REF`` and ``ALT``; inversions and duplications are written with symbolic alleles.
Config parameter ``vcf_cpx_symbolic`` writes complex variants with symbolic ``<CPX>`` alleles instead.

Records are sorted in reference order: chromosomes in the order of the ``##contig`` header lines,
which is the order of sequences in the reference FASTA, then by position. ``REF`` and ``ALT`` are
upper case, and IUPAC ambiguity codes are reduced to the base that is first alphabetically (e.g. ``R``
is written as ``A``), as VCF requires. A complex variant whose sequence matches the reference bases it
replaces describes no change and is left out of the VCF with a warning.

Multi-allelic records
---------------------

When haplotypes carry different alleles starting at the same position, PAV writes one record with an
``ALT`` per distinct allele (e.g. ``1|2``) rather than one record per allele (``1|0`` and ``0|1``),
which comparison tools and ``bcftools consensus`` do not reliably resolve. ``REF`` is the longest ``REF`` of the combined alleles, and shorter alleles are
extended with the reference bases it adds. ``ID`` lists the variant IDs separated by semicolons, and
per-allele INFO fields have one value per ``ALT``.

Records are combined only if they have the same ``FILTER``, so a passing allele is never filtered by
an allele on another haplotype. Records are not combined when one haplotype carries more than one of
them (e.g. an SNV and an insertion after the same base), and overlapping records starting at different
positions are not combined. Symbolic records are never combined.

Complex variants
----------------

A sequence-resolved complex variant replaces the reference bases skipped between its anchoring
alignments with the variant sequence (variant table columns ``seq_ref_pos``, ``seq_ref_end``, and
``seq``), so applying the record to the reference reconstructs the haplotype. When the anchoring
alignments overlap on the reference (a duplication at the outer breakpoint), nothing is deleted and the
variant sequence, which includes the duplicated bases, is inserted. ``INFO/SUBTYPE``,
``INFO/CPX_REF_TRACE``, ``INFO/CPX_QRY_TRACE``, and ``INFO/DUP`` describe the structure.

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

``INFO/SVLEN``, ``INFO/END``, and ``INFO/SEQ`` are written only for symbolic ALT alleles (``INV``,
``DUP``, and ``CPX`` with ``vcf_cpx_symbolic``). A sequence-resolved ``REF`` and ``ALT`` pair already
gives the length, the span, and the sequence. VCF defines ``SVLEN`` for symbolic structural variant
alleles and directs writers to use the missing value for all others, and VCF 4.4 and later treat
``END`` as a computed field. The table below describes the quantities wherever they are written.

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
     - Not written; duplications are derived variants (``FILTER`` ``DERIVED``)
   * - ``CPX``
     - Inserted bases (query, symbolic only)
     - Reference breakpoint span (symbolic only)
     - Inserted bases (query, symbolic only)

Three cases deserve attention:

``INV``
   ``SVLEN`` and ``END`` both give the inverted reference span, but ``INFO/SEQ`` holds the inverted
   allele taken from the query. An inversion carrying indels has an allele whose length differs from
   the reference span, so ``len(SEQ)`` is not a size field.

``CPX``
   For a symbolic ``<CPX>`` allele, a complex variant has a reference breakpoint part and an inserted
   sequence part, and the two are reported separately. ``POS`` and ``END`` span the reference bases deleted or duplicated at the
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
   the query region between them. They are written into ``ALT`` by reading the reference over the
   duplicated region, the same region reported by ``INFO/DUP``. Sequence divergence between the two copies is therefore not represented. An
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
semicolons and equals signs. Contig names embedded in INFO values (``INFO/QUERY``, ``INFO/DUP``, and
the complex variant traces) are percent-encoded where they contain one of these characters, a percent
sign, a tab, or a line break (e.g. ``;`` is written as ``%3B``). Fields holding more than one value
separate them with ``/``: ``INFO/DUP`` writes its regions that way.

``INFO/CPX_REF_TRACE`` and ``INFO/CPX_QRY_TRACE`` report segment lengths without thousands
separators. The ``trace_ref`` and ``trace_qry`` columns of the variant tables keep the separators
for readability, so the trace strings in the VCF and in the Parquet tables differ.
