# Error Threshold Derivation

We want a multiplicity cutoff $m^*$ such that any edge with count $\leq m^*$ is more likely a sequencing error than a real genomic k-mer — derived purely from numbers already in the JSON output.

---

## Step 1 — Reading off the real coverage

`mean_max_branch_ratio` is $C_{real}$.

At any branch, one outgoing edge is the true genome path (high multiplicity) and the other is an error (multiplicity ≈ 1). Their ratio is just $C_{real}$, so this feature directly encodes sequencing depth.

**Example:** `mean_max_branch_ratio = 29` means the genome was sequenced at ~29×.

---

## Step 2 — How many errors are in the graph?

Every edge is either a real k-mer (appears ~$C_{real}$ times) or an error k-mer (appears ~1 time). The overall mean $\bar{m}$ (`mult_mean`) is a weighted average of those two populations:

$$\bar{m} = \frac{C_{real} \cdot E_{real} + 1 \cdot E_{err}}{E_{real} + E_{err}}$$

Rearranging gives the error-to-real ratio:

$$r = \frac{E_{err}}{E_{real}} = \frac{C_{real} - \bar{m}}{\bar{m} - 1}$$

**Example:** $C_{real}=29$, $\bar{m}=12.44$:

$$r = \frac{29 - 12.44}{12.44 - 1} = \frac{16.56}{11.44} = 1.45$$

There are 45% more error edges than real genome edges — which makes sense at 29× coverage with ~2% error rate.

---

## Step 3 — Poisson rate of an error k-mer

A specific error k-mer only reappears if another read independently makes the same mistake at the same position. With per-base error rate $p$ and 3 possible wrong bases, that probability is $p/3$ per read. Across $C_{real}$ reads:

$$\lambda_e = C_{real} \cdot \frac{p}{3}$$

We don't observe $p$ directly, but we can infer it: the total number of distinct error k-mers is roughly $E_{err} \approx N_{reads} \cdot k \cdot p$, and real k-mers satisfy $E_{real} \approx N_{reads} \cdot k / C_{real}$, so:

$$p = \frac{r}{k \cdot C_{real}} \quad \Rightarrow \quad \lambda_e = \frac{r}{3k}$$

**Example:** $r=1.45$, $k=21$:

$$\lambda_e = \frac{1.45}{63} \approx 0.023$$

An error k-mer appears on average 0.023 times — almost always just once, occasionally zero.

---

## Step 4 — Where the two populations cross

Model edge counts as Poisson: real edges from $\text{Pois}(C_{real})$, error edges from $\text{Pois}(\lambda_e)$. The crossover $m^*$ is where the posterior probability of being real equals that of being an error:

$$E_{real} \cdot e^{-C_{real}} C_{real}^m = E_{err} \cdot e^{-\lambda_e} \lambda_e^m$$

Take $\ln$ (the $m!$ cancels), then solve for $m$:

$$\boxed{m^* = \frac{\ln(r) + (C_{real} - \lambda_e)}{\ln(C_{real} / \lambda_e)}}$$

**Example:** $r=1.45$, $C_{real}=29$, $\lambda_e=0.023$:

$$m^* = \frac{\ln(1.45) + (29 - 0.023)}{\ln(29 / 0.023)} = \frac{0.37 + 28.98}{7.14} \approx 4.1$$

Any edge with multiplicity $\leq 4$ is more likely noise than signal.

---

## Canonical k-mer note

MEGAHIT stores only the lexicographically smaller of a k-mer and its reverse complement. A substitution on the forward strand and the complementary substitution on the reverse strand can hash to the same canonical k-mer, so two independent reads can increment the same error edge — effectively doubling $\lambda_e$:

$$\lambda_e^{\text{canonical}} = \frac{r}{1.5k}$$

The current implementation uses the non-canonical formula (`3.0 * kmer_k`), which is conservative: it underestimates $\lambda_e$, so $m^*$ lands slightly lower, meaning we flag fewer edges as errors when uncertain. The difference is about 0.3 in $m^*$.
