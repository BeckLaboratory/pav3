"""Sphinx directives rendering PAV table schemas and VCF header fields as references.

Schemas are dictionaries mapping column names to Polars types (see :mod:`pav3.schema`). VCF header
fields are lists of field objects (see :mod:`pav3.vcf`). Rendering them at build time keeps the
documented columns and fields in step with the code instead of duplicating them into the
documentation by hand.

Usage::

    .. schema-table:: pav3.schema.ALIGN
    .. vcf-fields:: pav3.vcf.INFO_FIELDS
"""

__all__ = [
    'SchemaTableDirective',
    'VcfFieldsDirective',
    'setup',
]

import importlib
from typing import Any

from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.application import Sphinx


def _type_name(dtype: Any) -> str:
    """Render a Polars data type as a readable string.

    Polars types are sometimes classes and sometimes instances; nested types (List, Struct) carry
    their inner types in the repr.

    :param dtype: Polars data type, either a class or an instance.

    :return: Readable type name.
    """
    if isinstance(dtype, type):
        return dtype.__name__

    return str(dtype)


class SchemaTableDirective(Directive):
    """Render a schema dictionary as a two-column table of names and types."""

    required_arguments = 1
    has_content = False

    def run(self) -> list[nodes.Node]:
        """Build the table node.

        :return: List containing the table node.

        :raises RuntimeError: If the schema cannot be imported.
        """
        target = self.arguments[0]
        module_name, _, attr_name = target.rpartition('.')

        try:
            schema = getattr(importlib.import_module(module_name), attr_name)
        except (ImportError, AttributeError) as ex:
            raise RuntimeError(f'schema-table: cannot resolve "{target}": {ex}') from ex

        table = nodes.table()
        group = nodes.tgroup(cols=2)
        table += group
        group += nodes.colspec(colwidth=30)
        group += nodes.colspec(colwidth=70)

        header = nodes.thead()
        group += header
        header += self._row('Column', 'Type')

        body = nodes.tbody()
        group += body

        for name, dtype in schema.items():
            body += self._row(name, _type_name(dtype), literal=True)

        return [table]

    @staticmethod
    def _row(*cells: str, literal: bool = False) -> nodes.row:
        """Build a table row.

        :param cells: Cell text.
        :param literal: Whether to render cells as inline literals.

        :return: Row node.
        """
        row = nodes.row()

        for text in cells:
            entry = nodes.entry()
            node: nodes.Node = nodes.literal(text=text) if literal else nodes.Text(text)
            entry += nodes.paragraph('', '', node)
            row += entry

        return row


class VcfFieldsDirective(Directive):
    """Render a list of VCF header field objects as a table.

    Handles the field classes in :mod:`pav3.vcf`, which carry a name and description and, for INFO
    and FORMAT fields, a VCF ``Number`` and ``Type``. Columns absent from a field class are omitted.
    """

    required_arguments = 1
    has_content = False

    #: Field attribute names paired with their column headings.
    COLUMNS = (
        ('name', 'ID'),
        ('number', 'Number'),
        ('type_', 'Type'),
        ('description', 'Description'),
    )

    def run(self) -> list[nodes.Node]:
        """Build the table node.

        :return: List containing the table node.

        :raises RuntimeError: If the field list cannot be imported or is empty.
        """
        target = self.arguments[0]
        module_name, _, attr_name = target.rpartition('.')

        try:
            field_list = getattr(importlib.import_module(module_name), attr_name)
        except (ImportError, AttributeError) as ex:
            raise RuntimeError(f'vcf-fields: cannot resolve "{target}": {ex}') from ex

        if not field_list:
            raise RuntimeError(f'vcf-fields: "{target}" is empty')

        columns = [
            (attr, heading)
            for attr, heading in self.COLUMNS
            if hasattr(field_list[0], attr)
        ]

        table = nodes.table()
        group = nodes.tgroup(cols=len(columns))
        table += group

        for attr, _ in columns:
            group += nodes.colspec(colwidth=60 if attr == 'description' else 15)

        header = nodes.thead()
        group += header
        header += SchemaTableDirective._row(*[heading for _, heading in columns])

        body = nodes.tbody()
        group += body

        for field in field_list:
            row = nodes.row()

            for attr, _ in columns:
                text = str(getattr(field, attr))
                entry = nodes.entry()
                node: nodes.Node = (
                    nodes.Text(text) if attr == 'description' else nodes.literal(text=text)
                )
                entry += nodes.paragraph('', '', node)
                row += entry

            body += row

        return [table]


def setup(app: Sphinx) -> dict[str, Any]:
    """Register the directives.

    :param app: Sphinx application.

    :return: Extension metadata.
    """
    app.add_directive('schema-table', SchemaTableDirective)
    app.add_directive('vcf-fields', VcfFieldsDirective)

    return {
        'version': '1.0',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
