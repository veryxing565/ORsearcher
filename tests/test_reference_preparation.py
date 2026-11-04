"""Pure-Python reference preparation safety checks; no network or private data."""
from importlib.util import module_from_spec, spec_from_file_location
import hashlib
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch
import uuid
import zipfile

SPEC = spec_from_file_location("prepare_default_references", Path(__file__).resolve().parents[1] / "tools" / "prepare_default_references.py")
PREPARE = module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


class ReferencePreparationTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(__file__).resolve().parent / (".reference_test_" + uuid.uuid4().hex)
        self.folder.mkdir()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def test_fasta_count_uses_headers_not_description_characters(self):
        records = list(PREPARE.records(">one annotation=x>y\nMA?\n>two\nMCT\n"))
        self.assertEqual(records, [("one annotation=x>y", "MA?"), ("two", "MCT")])

    def test_refuse_different_destination(self):
        target = self.folder / "reference.faa"
        target.write_bytes(b"existing")
        new = b"new"
        with self.assertRaises(FileExistsError):
            PREPARE.install_bytes(target, new, hashlib.sha256(new).hexdigest())
        self.assertEqual(target.read_bytes(), b"existing")

    def test_repeated_identical_install_is_idempotent(self):
        target = self.folder / "reference.faa"
        data = b">id\nMAA\n"
        checksum = hashlib.sha256(data).hexdigest()
        PREPARE.install_bytes(target, data, checksum)
        PREPARE.install_bytes(target, data, checksum)
        self.assertEqual(target.read_bytes(), data)
        self.assertFalse(target.with_name(target.name + ".partial").exists())

    def test_archive_selection_hash_verified_and_no_path_extraction(self):
        raw = b">GR1\nMA?\n>OR1\nMNN\n"
        archive = self.folder / "source.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("supp/sequences/unaligned_GROR.fasta", raw)
            handle.writestr("../../should_not_exist", "not extracted")
        selection = [dict(source_id="GR1", output_header="refGR-GR1 family=GR", sequence_sha256=hashlib.sha256(b"MAX").hexdigest())]
        with patch.object(PREPARE, "SUPPLEMENT_SHA256", hashlib.sha256(archive.read_bytes()).hexdigest()), patch.object(PREPARE, "RAW_GROR_SHA256", hashlib.sha256(raw).hexdigest()):
            self.assertEqual(PREPARE.extract_gr(archive, selection), b">refGR-GR1 family=GR\nMAX\n")
            selection[0]["sequence_sha256"] = "wrong"
            with self.assertRaises(ValueError):
                PREPARE.extract_gr(archive, selection)
        self.assertEqual(sorted(p.name for p in self.folder.iterdir()), ["source.zip"])


if __name__ == "__main__":
    unittest.main()
