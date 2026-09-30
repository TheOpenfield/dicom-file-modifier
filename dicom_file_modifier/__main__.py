"""``python -m dicom_file_modifier <befehl> ...`` = ``dfm <befehl> ...`` (siehe ``cli``)."""

import sys

from .cli import main

sys.exit(main())
