"""PAV Python library routines."""

import importlib
import logging

__version__ = '3.0.0.dev26'

__all__ = [
    'align',
    'anno',
    'call',
    'cli',
    'const',
    'errors',
    'fig',
    'format',
    'io',
    'inv',
    'kde',
    'lgsv',
    'maketable',
    'params',
    'pipeline',
    'region',
    'schema',
    'seq',
    'util',
    'vcf',
    'workflow',
]

__license__ = 'MIT License'
__author__ = 'Peter Audano'
__copyright__ = 'Copyright (c) 2024 Peter Audano'
__credits__ = ['Peter Audano']
__maintainer__ = 'Peter Audano'
__email__ = 'peter.audano@jax.org'
__status__ = 'Development'

__license_spdx__ = 'MIT'

__license_text__ = """MIT License

Copyright (c) 2024 Peter Audano

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

for name in __all__:
    globals()[name] = importlib.import_module(f'.{name}', package=__name__)


# Configure logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
