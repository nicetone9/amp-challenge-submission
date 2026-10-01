# Training release licensing and attribution

The [training archive](../training-data.json) contains mixed-origin data.
Project owner redistribution authority was confirmed in this project chat on
2026-09-30 (America/New_York); no private conversation transcript is published.

| Material / rights | Release terms and attribution |
| --- | --- |
| Owner-controlled ZYL used-positive annotations and original curation rights | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); credit the AMP Challenge VQ-VAE/DiMA project dataset contributors and preserve original per-row source/evidence attribution. Applies only to rights the owner can grant. |
| Prepared mixed database rights | [ODbL 1.0](https://opendatacommons.org/licenses/odbl/1-0/), retaining potential Peptipedia-derived database share-alike obligations. Individual contents are not blanket relicensed. |
| OmegAMP supplied files | Repository [MIT notice at pinned c9500d3](https://github.com/szczurek-lab/OmegAMP/blob/c9500d3ca749f7db9f56409fc9868c2d09a62b3e/LICENSE) retained. Repository MIT does not override upstream component database/content rights. |
| UniProt-labelled sequences | Credit UniProt Consortium; [UniProt CC BY 4.0 terms](https://www.uniprot.org/help/license). Historical accession/version mapping was not retained in all inputs. |
| Peptipedia-labelled/composite sequences | Credit Cabas-Mora et al., Peptipedia v2.0, Database 2024, baae113, [doi:10.1093/database/baae113](https://doi.org/10.1093/database/baae113). Retain [database ODbL notice](https://github.com/kren-ai-lab/Peptipedia/blob/main/LICENSE-DATABASE.md); exact historical version/per-row component map is unknown. Web-platform license is distinct and is not substituted for database terms. |
| Internally AMPSphere-named project package | Four input filenames/hashes and per-row component membership are preserved. This composite includes UniProt background; the name does not establish an official AMPSphere catalog version or blanket license. Owner-controlled rights follow the first row; third-party contents retain upstream terms. |
| Original processing/training code | Project MIT notice. External model/code licenses remain separately scoped. |

Contains information from Peptipedia, which is made available under the
Open Database License. The prepared database is offered under ODbL 1.0 for
database rights; users must also preserve individual content attribution/terms.
The archive includes this source attribution, source/right records, OmegAMP MIT,
Peptipedia ODbL, project-code MIT and the CC BY legal-code reference.

Changes from inputs: canonical/intact length filtering, protected exact/H30/cluster/
noncanonical/motif exclusions, de-duplication, provenance/component joins, AMP
strata and cluster-wise generator splits. No measured-activity claim is added.

Owner confirmation resolves the previously unapproved nonpublic release; it does
not independently establish all third-party rights or recover missing historical
accessions. This is a transparent mixed-license release, not a representation that
every component has one permissive license or that organizers accepted the data
licensing/workflow. See [exact scope and remaining limitations](DATA.md).
