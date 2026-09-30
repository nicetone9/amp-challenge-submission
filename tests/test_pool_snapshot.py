import json
import tempfile
import unittest
from pathlib import Path
from amp_submission.prepare_reference import expansion_record


class PoolSnapshotTests(unittest.TestCase):
    def record(self):
        return {"status": "RAW_BATCHES_VERIFIED_CLOUD_SYNC_PENDING",
                "per_branch": {"vq": 300000, "dima": 300000}, "raw_requested": 600000,
                "pool_unique": 597555, "batch_hash_list_sha256": "fixed"}

    def read_snapshot(self, record):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text(json.dumps(record))
            result, _ = expansion_record(directory, path)
            self.assertEqual(json.loads(path.read_text()), record)
            return result

    def test_explicit_snapshot_keeps_pending_cloud_status(self):
        record = self.read_snapshot(self.record())
        self.assertEqual(record["additional_counts"], {"vq": 250000, "dima": 250000})
        self.assertEqual(record["status"], "RAW_BATCHES_VERIFIED_CLOUD_SYNC_PENDING")

    def test_unbalanced_and_bad_count_rejected(self):
        for counts, raw in [({"vq": 300000, "dima": 299999}, 599999),
                            ({"vq": 300000, "dima": 300000}, 500000),
                            ({"vq": 1, "dima": 1}, 2)]:
            record = self.record()
            record.update(per_branch=counts, raw_requested=raw)
            with self.assertRaises(ValueError):
                self.read_snapshot(record)

    def test_missing_hash_and_wrong_status_rejected(self):
        for key, value in [("batch_hash_list_sha256", ""), ("status", "SUCCESS")]:
            record = self.record()
            record[key] = value
            with self.assertRaises(ValueError):
                self.read_snapshot(record)

    def test_default_uses_original_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "complete.json"
            record = {"status": "RAW_POOL_ONLY_NOT_QUALIFIED",
                      "additional_counts": {"vq": 100000, "dima": 100000}}
            path.write_text(json.dumps(record))
            self.assertEqual(expansion_record(directory), (record, path))


if __name__ == "__main__":
    unittest.main()
