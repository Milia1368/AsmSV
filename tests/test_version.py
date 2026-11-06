"""Keep installed, source and command-line version identities consistent."""
from pathlib import Path
import os
import re
import subprocess
import sys
import unittest

from asmsv_next import __version__


class VersionTests(unittest.TestCase):
    def test_project_metadata_and_cli_agree_with_imported_package(self):
        root = Path(__file__).resolve().parents[1]
        metadata = (root / "pyproject.toml").read_text()
        self.assertEqual(re.search(r'^version = "([^"]+)"', metadata, re.M)[1], __version__)
        env = dict(os.environ, PYTHONPATH=str(root / "src"))
        result = subprocess.run(
            [sys.executable, "-B", "-m", "asmsv_next.cli", "--version"],
            env=env, text=True, capture_output=True, check=True,
        )
        self.assertIn(__version__, result.stdout)


if __name__ == "__main__":
    unittest.main()
