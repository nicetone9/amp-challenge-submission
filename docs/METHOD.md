# Method and ranking

## Models and RL

VQ codec:1280→128 CNN,2 residual blocks,K512,code dimension64,EMA0.99,
commitment0.25. Prior:6 layers,256hidden,8heads. Codec frozen for RL.
DiMA:12 layers,320hidden,16heads;1280latent,tan10,self-conditioning,padding mask;
decoder and reference frozen. Current challenge exports use50 stochastic reverse
steps, not the original baseline plan's500. Do not compare sampling cost as if budgets matched.

Six separate policies per architecture: broad,gram_positive,gram_negative,mdr,
selectivity,joint. Two seeds42/43,20 updates each:24 runs total.
Best checkpoint chosen by validation; step0 can be retained if no improvement.
The independent audit sample does not serve checkpoint selection. The functional
predictor is shared with the reward, so this is not an independent functional oracle.

## Proxy scores

For each predictor, fixed calibration values define a midpoint empirical CDF:
p(x)=(number of values<x + number<=x)/(2N).
ANIA reward=1−p(log-MIC output); HemoPI2 reward=p(predicted HC50 in µM).

- Broad:mean EC,PA,SA rewards.
- Gram+:SA reward.
- Gram−:mean EC,PA.
- MDR proxy:minimum EC,PA,SA reward; **not an MDR isolate predictor**.
- Selectivity proxy:half broad+half HC50 reward; **not HC50/MIC50**.
- Joint:arithmetic mean of the preceding five, range0–1; **not an official score**.

ANIA outputs remain in their native log-MIC scale. No unsupported unit or logarithm
conversion is used to claim MIC≤16µM. HemoPI2 HC50 predictions are not measured HC50.

## Library and Top100

Seed42 samples length from weighted train PDF conditioned to8–40.
Six policies alternate by128-member batches. Reject duplicates,noncanonical sequences,
out-of-range length,maxsingle-AAfraction>0.6,or reference Levenshtein ratio>0.8.
Stop at50,000 accepted,or fail at500,000 raw. No duplication/filling to reach size.

A fixed RNG4202048 chooses2,048 accepted members without replacement.
Greedy Top100 utility=joint proxy−0.15×maximum Indel similarity to previously selected.
First choice has no diversity penalty; ties retain earliest pool index.
No manual substitution. Top100 is not a global score-only sort of the full50k.
Two architectures stay separate; Top200 only concatenates their respective lists.

## Interpretation

Most paired RL confidence intervals cross zero. Raw baseline,raw six-policy mixture,
joint-only,filtered library andTop100 are separate populations.
Higher score afterTop100 filtering is not evidence of a causal RL benefit.
CD50/CD95=CD-HITcluster count/n,global identity50%/95%,bilateral coverage80%.
FD/MMD use independent ProtT5 and finite held-out references,not activity.
