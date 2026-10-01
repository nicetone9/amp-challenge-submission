# Current submission checklist

Current candidate: the new seed-42 mixed 600,000-attempt pool, selected into a
50,000-member library and ranked Top100. Historical selection reports remain
separate. No Kaggle submission or organizer acceptance is claimed.

| Minimum benchmark requirement | Current evidence |
| --- | --- |
| Method abstract | [ABSTRACT.md](ABSTRACT.md) |
| 50,000 designs and ranked Top100 | Source-finalization FASTAs and independent hard/metric checks pass |
| Selection, ranking, filters, training and external-resource summary | [WRITEUP.md](WRITEUP.md), [data disclosure](../DATA.md) |
| GitHub inference code and model weights | Public repository and seed-42 release |
| Organizer read access | Public repository/assets readable anonymously |

| Additional full requirement | Current evidence / gap |
| --- | --- |
| Public code, weights and usage docs | Published; [public cold-clone filtering certificate passes](../../reports/frozen-pool-seed42-public-verification.json) |
| Permissive OSI-approved code license | Original code MIT; full upstream MIT/BSD notices retained |
| Python version and uv.lock | Python 3.11; locked default installation |
| Default uv run generate with seed 42 | Disclosed frozen-pool selection route; all arguments have defaults |
| Identical repeated outputs | Two public cold clones pass with byte-identical, source-hash-matched outputs |
| Full training-data disclosure / permissive nonpublic-data release | Incomplete; no redistribution authority assumed |

[The complete unchanged official validator passed](../../reports/official-template-seed42-public-verification.json)
on public revision `aea6ce7ab7b50bd517b7788ed30037f83bafb914`: a third fresh
clone, dependency installation, two default generations, sequence/reference checks
and byte identity. Subsequent documentation-only changes do not alter its runtime inputs.
A validator pass does not assess training-data rights, experimental quality or
organizer acceptance of a fixed-pool workflow. Fresh-model resampling identity
remains uncertified. Full co-authorship eligibility is not claimed.
