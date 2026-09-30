# Competition alignment (checked 2026-09-29)

Sources: [official template](https://github.com/szczurek-lab/amp-challenge-2027),
[category definitions](https://szczurek-lab.github.io/amp-challenge-website/#categories),
[Kaggle overview](https://www.kaggle.com/competitions/amp-challenge/overview).

## Five experimental categories

| Official category | Primary criterion | Tie-break / gate | Current evidence |
| --- | --- | --- | --- |
| Broad-spectrum | Success fraction across20strains | MIC90 over20 | EC,PA,SA species-level proxies only; no20-strain MIC panel |
| Gram-positive | Success fraction across5Gram+strains | MIC50 overGram+ | SA proxy only; no5-strain measurements |
| Gram-negative | Success fraction across15Gram−strains | MIC50 overGram− | EC,PA proxies only; no15-strain measurements |
| MDR ESKAPE | Success fraction across8MDRisolates | MIC50 overMDR | NoMDR-specific oracle; worst-of-three proxy only |
| Optimal selectivity | HC50/MIC50 across20strains | AtleastoneMIC≤16µM; inactiveexcluded; cappedHC50tie broken byMIC50 | HC50 prediction+rank proxy only; official safetywindow unavailable |

Activity means measured MIC≤16µM. MIC assay limit64µM; no inhibition is>64.
HC50 assay limit128µM; nonhemolytic is>128,not128.
Team results average peptide metrics over25randomly sampled candidates fromTop100.
Individual-peptide proxy scores cannot establish team category results.

## Sequence and delivery gates

| Requirement | Status |
| --- | --- |
| 50,000 unique sequences per library | Historical exports verified, bothroutes |
| RankedTop100 drawn fromownlibrary | Verified, bothroutes |
|20standardAA,length8–50 | Verified; actualgeneration8–40 |
|Linearity,freeN/Ctermini,noamidationorothermodification | Design specification; synthesis/QC unverified |
|Librarynoexactofficialreferenceoverlap | Verified byhash-matched historical full-library checker |
|Top100maxreferenceLevenshteinratio≤0.8 | Recomputed for200/200; full50kbothuse same stricterfilter |
|Method,filters,externaltools,rankingdocumented | Included |
|PrivateGitHubinference+weights+organizerreadaccess | Codeincluded; weightsnotreleased; accessnotgranted |
|Publicrepo+permissivecode | MIToriginalcode; reporemainsprivatebyauthorization |
|Fulltrainingdatadisclosure+nonpublicdatarelease | Partial; sourcepermissionsunresolved |
|Pythonversion+uv.lock+noargentrypoint | Implemented; externalR/oracleenvironment remainsreleasegate |
|Twoclean-clonefullgenerationsidentical | Historicalsameenvironmentpassed; portablecold-cloneGPUtestpending |
|seqmecomputationalscreening,DBAASP/dbAMP/APDnearhits | Partialinternalproxies; fullofficialscreeningnotrun |
|SynthesisandpurityidentityQC | Notassessed |
|Officialacceptance/Kagglesubmission | Notperformed |

No numeric overall compliance percentage: the gates are heterogeneous and some are
blocking requirements, not interchangeable points.
