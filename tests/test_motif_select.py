import unittest
from amp_submission.motif_select import motif_top

class MotifTopTests(unittest.TestCase):
    def rows(self):
        return [dict(sequence=str(i), cluster=str(i), tier="high" if i%2 else "low",
                     ranking=i/30, weakest=.3, motif_support=int(i<8),
                     short_motif_support=int(i<8)) for i in range(20)]

    def test_quota_and_tiers(self):
        result=motif_top(self.rows(),10,4)
        self.assertEqual(len({r["sequence"] for r in result}),10)
        self.assertEqual(sum(r["selection_reason"]=="short_motif_quota" for r in result),4)
        self.assertEqual(sum(r["tier"]=="high" for r in result),5)
        self.assertGreaterEqual(sum(bool(r["short_motif_support"]) for r in result),4)

    def test_input_order_determinism(self):
        self.assertEqual(motif_top(self.rows(),10,4),motif_top(self.rows()[::-1],10,4))

    def test_missing_support_fails(self):
        rows=self.rows()
        for r in rows:r["short_motif_support"]=0
        with self.assertRaises(ValueError):motif_top(rows,10,4)

    def test_bad_or_duplicate_input_fails(self):
        with self.assertRaises(ValueError):motif_top(self.rows(),10,3)
        rows=self.rows();rows.append(rows[0])
        with self.assertRaises(ValueError):motif_top(rows,10,4)
