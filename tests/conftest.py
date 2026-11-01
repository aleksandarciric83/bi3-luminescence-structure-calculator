import mimetypes
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Some locked-down systems forbid reading the Windows registry, which Streamlit's media handling does
# through mimetypes. The tests do not need the registry, so it is skipped when it cannot be read.
try:
    mimetypes.guess_extension("image/png")
except PermissionError:
    mimetypes.inited = True
    mimetypes._db = mimetypes.MimeTypes()
