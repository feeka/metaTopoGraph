# Bayesian Error Threshold Derivation

## Goal

Given only the graph output, find the multiplicity $m^*$ where:

> "An edge with multiplicity $\leq m^*$ is more likely a sequencing error than a real genome k-mer."

---

## Step 1 — What is real genome coverage?

`mean_max_branch_ratio` = $C_{real}$

At a branching node, the dominant outgoing edge is the real genome path (high multiplicity), and the other is an error branch (multiplicity ≈ 1). The ratio between them is $C_{real} / 1 = C_{real}$.

**Example:** `mean_max_branch_ratio = 29` → real coverage = 29×.

---

## Step 2 — How many error edges are there?

The graph has two kinds of edges mixed together:
- **Real** edges: multiplicity ≈ $C_{real}$, there are $E_{real}$ of them
- **Error** edges: multiplicity ≈ 1, there are $E_{err}$ of them

The overall mean is just a weighted average:

$$\text{mult\_mean} = \frac{C_{real} \cdot E_{real} + 1 \cdot E_{err}}{E_{real} + E_{err}}$$

Solve for the ratio $E_{err}/E_{real}$:

$$\frac{E_{err}}{E_{real}} = \frac{C_{real} - \text{mult\_mean}}{\text{mult\_mean} - 1}$$

**Example:** $C_{real}=29$, `mult_mean=12.44`:

$$\frac{E_{err}}{E_{real}} = \frac{29 - 12.44}{12.44 - 1} = \frac{16.56}{11.44} = 1.45$$

There are 1.45× more error edges than real edges in the graph.

---

## Step 3 — What is the Poisson rate of error edges?

A specific error k-mer (e.g. the 21-mer that results from the mutation A→C at position 7 of read 1234) appears in a new read only if that read has the **exact same mutation at the exact same position**. That probability is:

$$\frac{p}{3}$$

where $p$ = per-base error rate, and $3$ = number of wrong bases possible at any position (DNA has 4 bases, 1 is correct, 3 are wrong).

So the Poisson rate for a specific error k-mer is:

$$\lambda_e = C_{real} \cdot \frac{p}{3}$$

We don't know $p$ directly, but we can get it from:

$$E_{err} \approx N_{obs} \cdot k \cdot p \quad \Rightarrow \quad p = \frac{E_{err}/E_{real}}{k \cdot C_{real}} = \frac{\text{ratio}}{k \cdot C_{real}}$$

Substituting back:

$$\lambda_e = C_{real} \cdot \frac{p}{3} = \frac{\text{ratio}}{3k}$$

**Example:** ratio=1.45, k=21:

$$\lambda_e = \frac{1.45}{3 \times 21} = 0.023$$

---

## Step 4 — The crossover multiplicity $m^*$

Model each edge's count as Poisson. An edge with multiplicity $m$ is equally likely to be real or error when their Poisson likelihoods match, weighted by how many of each kind exist:

$$E_{real} \cdot \text{Poisson}(m \mid \lambda_r) = E_{err} \cdot \text{Poisson}(m \mid \lambda_e)$$

Take $\ln$ of both sides (the $m!$ cancels):

$$\ln E_{real} + m \ln\lambda_r - \lambda_r = \ln E_{err} + m \ln\lambda_e - \lambda_e$$

Solve for $m$:

$$\boxed{m^* = \frac{\ln(E_{err}/E_{real}) + (\lambda_r - \lambda_e)}{\ln(\lambda_r / \lambda_e)}}$$

**Example:** ratio=1.45, $\lambda_r=29$, $\lambda_e=0.023$:

$$m^* = \frac{\ln(1.45) + (29 - 0.023)}{\ln(29 / 0.023)} = \frac{0.37 + 28.98}{7.14} = \frac{29.35}{7.14} = 4.1$$

→ **Edges with multiplicity ≤ 4 are more likely errors than real genome sequence.**

---

## Canonical graph note

MEGAHIT stores the lexicographically smaller of a k-mer and its reverse complement. A substitution on the forward strand and the complementary substitution on the reverse strand can produce the **same canonical k-mer**, so two independent reads can create the same error edge. This doubles $\lambda_e$:

$$\lambda_e^{\text{canonical}} = \frac{\text{ratio}}{1.5k}$$

The implementation uses `3.0 * kmer_k` (not canonical-corrected) — the correction is small (~0.3 in $m^*$) and conservative (errs toward flagging fewer edges as errors).
