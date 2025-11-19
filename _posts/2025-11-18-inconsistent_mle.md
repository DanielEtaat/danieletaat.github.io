---
title: 'A Nice Example of an Inconsistent MLE'
date: 2025-11-18
permalink: /posts/2025/11/inconsistent_mle/
tags:
  - statistics
---

<style>
.bbox {
  border: 1px solid black;
  padding: 10px;
  background: #fdfdfd; /* optional: subtly distinguishes proof box */
}

.bbox p {
  margin-bottom: 0;
}
</style>


An important result from classical statistical theory is that the Maximum Likelihood Estimator (MLE) is consistent. In particular you will commonly come across results of the following form:

<div id="theorem" class="bbox" markdown="1">

**Theorem:**
Let $$\{p_\theta : \theta \in \Theta\}$$ be a collection of probability densities. If the collection satisfies condition **________________________**, then the maximum-likelihood estimator $$\hat{\theta}_n$$ is consistent:

$$
\hat{\theta}_n
\;\overset{p_\theta}{\longrightarrow}\;
\theta.
$$

</div>

<br>
You can fill in the blank with your favorite conditions. For example, Wald's proof of consistency would have you assume that $$\Theta$$ is compact and that the score function is dominated (among other things). The result I'm specifically interested in is the following: 

<div id="theorem" class="bbox" markdown="1">

**Theorem:**
Let $$\{p_\theta : \theta \in \Theta\}$$ be a collection of probability densities and $$X_1, \dots, X_n \overset{iid}{\sim} p_{\theta_0}$$ for some $$\theta_0 \in \Theta$$. If the collection satisfies:

* (A1) $$\,$$ $$\Theta$$ is an open interval in $$\mathbb R$$. 

* (A2) $$\,$$ The densities have common support.

* (A3) $$\,$$ The collection is identifiable (no two densities correspond to the same measure).

* (A4) $$\,$$ The likelihood $$\ell_n(\theta) := \sum_{i=1}^n \log p_\theta(X_i)$$ is continuously differentiable in $$\theta$$.

* (A5) $$\,$$ $$\mathbb P_{\theta_0}\bigl(\ell_n' \text{ has a unique root in } \Theta\bigr) \to 1$$ as $$n \to \infty$$.

Then the maximum-likelihood estimator $$\hat{\theta}_n$$ is consistent:

$$
\hat{\theta}_n
\;\overset{p_{\theta_0}}{\longrightarrow}\;
\theta_0.
$$

</div>

<br>
I like this result because the assumptions seem fairly mild. Conditions (A2), (A3), and (A4) are all standard. It would be hard to come up with a consistency proof without using them. Condition (A1) is a bit stronger than usual (most proofs only assume that $$\Theta$$ is open). The most unusual (and strongest) assumption is Condition (A5). This made me wonder if Condition (A1) is really necessary. Since Condition (A5) seems so strong, do we really need to assume that $$\Theta$$ is an open interval? Or can we relax this assumption to the more standard assumption that $$\Theta$$ is only an open subset of $$\mathbb R$$?

It turns out we cannot. My aim in this blog post will be to construct an example where $$\Theta \subset \mathbb R$$ is open and Conditions (A2)–(A5) hold, yet the MLE is not consistent. 

## 1. The Setup

Let’s build the counterexample. We will consider a one-parameter family of normal distributions, but with a twist: near zero the parameter controls the mean, and far from zero it controls the variance. Define the parameter space

$$
\Theta = (-\infty, -2) \cup (-1, 1) \cup (2, \infty)
$$

and for each $$\theta \in \Theta$$ define the distribution of a single observation $$X$$ by

$$
X \mid \theta \sim
\begin{cases}
\mathcal N(\theta, 1), & \theta \in (-1,1), \\
\mathcal N\bigl(0, 1 + 1/\theta\bigr), & \theta \in \Theta \setminus (-1,1).
\end{cases}
$$

Let the true parameter be
$$
\theta_0 = 0,
$$
so that the true distribution is

$$
X_1,\dots,X_n \overset{\text{iid}}{\sim} \mathcal N(0,1).
$$

So under the truth this is just a standard normal sampling model; all the weirdness lives in how we parametrize the family.


## 2. Controlling the number of roots

The key point is that the likelihood in this model is “glued together” from two familiar likelihoods:

* For $$\theta \in (-1,1)$$, the model looks like a normal model with unknown mean and known variance. The MLE in that submodel is the sample mean

$$
\hat\theta_{1,n} = \bar X := \frac{1}{n}\sum_{i=1}^n X_i.
$$

* For $$\theta \in \Theta \setminus (-1,1)$$, the model looks like a normal model with known mean and unknown variance $$\sigma^2(\theta) = 1 + 1/\theta$$. The MLE in that submodel is

$$
\hat\theta_{2,n}
= \left(\frac{1}{n}\sum_{i=1}^n X_i^2 - 1\right)^{-1}
= \bigl(s_n^2 - 1\bigr)^{-1},
$$

where $$s_n^2 := \frac{1}{n}\sum_{i=1}^n X_i^2$$ is the empirical second moment.

In both of these simpler submodels, the MLE is the unique root of the corresponding likelihood equation. So for our glued model, whenever the event 

$$
A_n = \{\hat\theta_{1,n} \in (-1,1)\   \text{ and } \ \hat\theta_{2,n} \in \Theta \setminus (-1,1)\}
$$

occurs, the likelihood equation for the full model has exactly two roots, at $$\hat\theta_{1,n}$$ and $$\hat\theta_{2,n}$$, and the MLE must be one of these two. So if we can show that $$A_n$$ occurs with high probability, then we will have shown that the likelihood equation has two roots with high probability. Condition (A5) requires there to be one root with high probability so it seems we are violating that condition. This will be amended later on. For now we will try to bound the number of roots at two.

To do this define the event

$$
A'_n
=\left\{ |\hat\theta_{1,n}| < n^{-1/4},\ |\hat\theta_{2,n}| > n^{1/4} \right\}.
$$

Under the true model $$X_i \sim \mathcal N(0,1)$$, the CLT implies that

$$
\hat\theta_{1,n} = \bar X = O_p(n^{-1/2}),
\qquad
\hat\theta_{2,n}^{-1} = s_n^2-1 = O_p(n^{-1/2}).
$$

Given these facts, it’s easy to see that
$$
\mathbb P(A'_n) \to 1
$$
as $$n\to\infty$$. Moreover, on the event $$A'_n$$, for large enough $$n$$ we have that $$\hat\theta_{1,n}\in(-1,1)$$ and $$\hat\theta_{2,n}\in\Theta\setminus(-1,1).$$ In other words we have that $$A'_n \subset A_n$$ for large enough $$n$$ and therefore that $$\mathbb P(A_n) \to 1$$. So the likelihood equation has exactly two roots at $$\hat\theta_{1,n}$$ and $$\hat\theta_{2,n}$$ and the MLE is whichever one has the larger log-likelihood with probability tending to $$1$$.


## 3. Step 1: Comparing the Two Log-Likelihoods

We will complete the rest of the construction in three steps:

1. Show that, asymptotically, each of the two candidates wins about half the time:
   $$
   \mathbb P_{\theta_0}\bigl(\ell_n(\hat\theta_{1,n}) < \ell_n(\hat\theta_{2,n})\bigr)\to 0.5.
   $$
2. Restrict the parameter space to a smaller open set $$\Theta_0\subset\Theta$$ that contains only $$\hat\theta_{1,n}$$ with high probability.
3. Show that on the restricted parameter space, the MLE falls outside of $$(-1, 1)$$ with probability not tending to $$0$$. This gives our counterexample.

Let’s go through the first step.

On $$A_n$$ we can write the log-likelihoods at the two candidate roots in simple forms. Doing the usual algebra for the normal likelihood gives

$$
\ell_n(\hat\theta_{1,n}) - \ell_n(0)
= -\frac12 \sum_{i=1}^n\Bigl[(X_i - \bar X)^2 - X_i^2\Bigr]
= \frac{n}{2}\bar X^{2},
$$

and

$$
\ell_n(\hat\theta_{2,n}) - \ell_n(0)
= -\frac{n}{2}\bigl(\log s_n^2 + 1 - s_n^2\bigr).
$$

Let
$$
t_n := s_n^2 - 1.
$$
Then $$s_n^2 = 1 + t_n$$ and

$$
\ell_n(\hat\theta_{2,n}) - \ell_n(0)
= -\frac{n}{2}\Bigl(\log(1+t_n) + 1 - (1+t_n)\Bigr)
= \frac{n}{2}\Bigl[t_n - \log(1+t_n)\Bigr].
$$

Using the Taylor expansion $$\log(1+t) = t - \frac{t^2}{2} + o(t^2)$$ we obtain
$$
t - \log(1+t)
= \frac{t^2}{2} + o(t^2),
$$ 
so

$$
\ell_n(\hat\theta_{2,n}) - \ell_n(0)
= \frac{n}{4}t_n^2 + o_p(1)
= \frac{n}{4}(s_n^2 - 1)^2 + o_p(1).
$$

Under $$\mathcal N(0,1)$$, the central limit theorem gives

$$
\sqrt{n}\bar X \ \xrightarrow{d}\ \mathcal N(0,1),
\qquad
\sqrt{n}(s_n^2 - 1) \ \xrightarrow{d}\ \mathcal N(0,2).
$$

A little moment calculation shows that

$$
\operatorname{Cov}\Bigl(\sqrt{n}\bar X,\ \sqrt{n}(s_n^2 - 1)\Bigr)
= \mathbb E[X^3] = 0,
$$

so the bivariate CLT says that these two quantities are jointly asymptotically normal with zero covariance. Hence they are asymptotically independent.

Now look at the two log-likelihood differences:

$$
\ell_n(\hat\theta_{1,n}) - \ell_n(0)
= \frac{n}{2}\bar X^{2},
\qquad
\ell_n(\hat\theta_{2,n}) - \ell_n(0)
= \frac{n}{4}(s_n^2-1)^2 + o_p(1).
$$

By continuous mapping,

$$
\frac{n}{2}\bar X^{2} \ \xrightarrow{d}\ \tfrac12 \chi^2_1,
\qquad
\frac{n}{4}(s_n^2-1)^2 \ \xrightarrow{d}\ \frac{1}{4}\cdot 2\,\chi^2_1 = \tfrac12 \chi^2_1,
$$

and the two limiting variables are independent. Then let
$$
U,V \stackrel{\text{iid}}{\sim} \tfrac12\chi^2_1,
$$
so that

$$
\bigl(\ell_n(\hat\theta_{1,n}) - \ell_n(0),\,\ell_n(\hat\theta_{2,n}) - \ell_n(0)\bigr)
\ \xrightarrow{d}\ (U,V).
$$

Since the law of $$U-V$$ is symmetric about zero,

$$
\mathbb P(U < V) = \mathbb P(U > V) = \tfrac12.
$$

By convergence of the joint distribution and the continuous mapping theorem we conclude that

$$
\mathbb P_{\theta_0}\bigl(\ell_n(\hat\theta_{1,n}) < \ell_n(\hat\theta_{2,n})\bigr)
\longrightarrow
\frac{1}{2}.
$$

So even in the original parameter space $$\Theta$$, the MLE does not concentrate on a single root of the likelihood equation. Half the time it picks the mean-type root, and half the time it picks the variance-type root, which diverges.


## 4. Step 2: Trimming the Parameter Space

The example so far already shows a pathology: the MLE is "pulled" toward large 

$$|\hat\theta_{2,n}| \sim \sqrt{n}$$ 

about half the time. But we also want an example in which the likelihood equation has a unique root with probability tending to $$1$$, as in Condition (A5).

To do this, we will carve out a smaller open subset $$\Theta_0\subset\Theta$$ that removes the second root $$\hat\theta_{2,n}$$ from the parameter space with high probability. Since we know that $$\hat\theta_{2,n}$$ will appear further and further away from $$0$$ as $$n \to \infty$$ we can achieve our goal by making $$\Theta_0$$ sparser as we move out to $$\infty$$. To do this let $$\{q_i\}_{i=1}^\infty$$ be an enumeration of the rational numbers in $$\Theta \setminus (-1,1)$$ (so each $$q_i$$ lies either below $$-2$$ or above $$2$$). Define

$$
\Theta_0
:= (-1,1)
\cup
\bigcup_{i=1}^\infty \bigl(q_i - 2^{-i}, q_i + 2^{-i}\bigr)\cap\Theta.
$$

Notice the following key properties of $$\Theta_0$$:

* $$\Theta_0$$ is open.

* $$\Theta_0$$ still contains $$(-1,1)$$, so a neighborhood of the true parameter, so it will contain $$\hat\theta_{1,n}$$ on $$A_n$$.

* Its Lebesgue measure $$\lambda(\Theta_0)$$ is bounded, and therefore the tails become thinner and thinner in measure:

  $$
  \lambda\bigl(\Theta_0\setminus(-M,M)\bigr)\ \downarrow\ 0
  \quad\text{as } M\to\infty.
  $$

Since on $$A_n$$, $$\hat\theta_{2,n} = (s_n^2-1)^{-1} \xrightarrow{p} \infty$$ and $$\mathbb P_{\theta_0}$$ is dominated by the Lebesgue measure, it is not hard to see that $$\mathbb P_{\theta_0}(B_n)\to 1$$ where 

$$
B_n := \{\hat\theta_{2,n}\notin \Theta_0\} \cap A_n.
$$

On $$B_n$$ the likelihood equation cannot have a root at $$\hat\theta_{2,n}$$ anymore (it is not in the parameter space). The only remaining root, when it exists, must be the mean-type root $$\hat\theta_{1,n} = \bar X \in(-1,1)$$. Hence on $$\Theta_0$$ the likelihood equation has a unique root with probability tending to $$1$$.


## 5. Step 3: Verifying Condition (A5) and Inconsistency

Up to this point we have:

* Constructed an open parameter space $$\Theta_0 \subset \mathbb R$$.
* Defined a model on $$\Theta_0$$ which satisfies (A2)--(A5). In particular we showed that on $$\Theta_0$$ the likelihood equation has a unique root at $$\hat\theta_{1,n}$$ with probability tending to $$1$$. 

What we will show now is that this root is not the MLE with probability not tending to $$0$$. In fact:

$$
\mathbb P_{\theta_0}\Bigl(\exists\theta\in\Theta_0 \text{ such that }
\ell_n(\theta) > \ell_n(\hat\theta_{1,n})\Bigr)
\longrightarrow
\frac{1}{2}.
$$

Why is this true? On the event $$\{\ell_n(\hat\theta_{2,n}) > \ell_n(\hat\theta_{1,n})\}$$ (which has probability tending to $$0.5$$), the point $$\hat\theta_{2,n}$$ maximizes the likelihood over the full parameter space $$\Theta$$. The log-likelihood $$\ell_n(\theta)$$ is continuous in $$\theta$$, so there is some neighborhood $$U$$ of $$\hat\theta_{2,n}$$ on which

$$
\ell_n(\theta) > \ell_n(\hat\theta_{1,n}) \quad \text{for all } \theta \in U.
$$

The rationals are dense, and for every rational $$q_i$$ we have included a little interval $$\bigl(q_i-2^{-i},q_i+2^{-i}\bigr)$$ inside $$\Theta_0$$. Therefore, any nonempty interval intersects $$\Theta_0$$. In particular we can pick a point
$$
\theta^\star \in U \cap \Theta_0.
$$
By construction,

$$
\ell_n(\theta^\star) > \ell_n(\hat\theta_{1,n}),
$$

and so on this event there exists a parameter in $$\Theta_0$$ with larger likelihood than the unique root $$\hat\theta_{1,n}$$. Since the event $$\{\ell_n(\hat\theta_{2,n}) > \ell_n(\hat\theta_{1,n})\}$$ has asymptotic probability $$0.5$$, the same holds for the event above (actually this implies that the probability of the event above is asymptotically bounded below by $$0.5$$, it's a nice exercise to explain why it is bounded above by $$0.5$$ as well).

From here we can actually conclude that the MLE will be inconsistent. Why? Because if the MLE does not equal $$\hat\theta_{1,n}$$ it must be outside of $$(-1, 1)$$. This is because $$\hat\theta_{1,n}$$ is the maximizer of likelihood on $$(-1, 1)$$. Then with probability not tending to $$0$$ we have that the MLE falls outside of $$(-1, 1)$$ and therefore that it is not consistent. 

This completes the construction. 



## 6. Conclusion

Let’s summarize what this example actually shows.

We started with a perfectly regular-looking one-dimensional model: common support, identifiability, smooth likelihood, and a natural-looking MLE. By gluing together a “mean branch” near zero and a “variance branch” far away, we engineered a situation in which the likelihood equation has two roots: one near the truth ($$\hat\theta_{1,n} = \bar X \to 0$$) and one that diverges ($$\hat\theta_{2,n}\sim \sqrt{n}$$). Step 1 showed that each of these two roots wins the likelihood battle about half the time. 

Then we sparsified the parameter space, passing from $$\Theta$$ to $$\Theta_0$$ by sprinkling small intervals around rationals in the tails. This construction eliminates $$\hat\theta_{2,n}$$ from the parameter space with high probability, while keeping points arbitrarily close to it inside $$\Theta_0$$. As a result, on $$\Theta_0$$:

* the likelihood equation has a unique root (at $$\hat\theta_{1,n}$$) with probability tending to 1, so Condition (A5) holds;
* but with asymptotic probability $$1/2$$ there still exists some $$\theta\in\Theta_0$$ with strictly larger likelihood than $$\hat\theta_{1,n}$$, forcing the MLE to lie outside $$(-1,1)$$ and therefore to be inconsistent.

In other words, on $$\Theta_0$$ Conditions (A2)–(A5) all hold, yet the MLE fails to converge to $$\theta_0$$.

### How strong is (A5) without (A1)?

One way to read this example is as a warning about Condition (A5). At first glance, requiring that the score have a **unique root with high probability** sounds extremely strong. But this construction shows that, on its own, (A5) is not very restrictive if we are allowed to choose an arbitrary open set $$\Theta \subset \mathbb R$$. Given any model where the likelihood equation has multiple roots, we can always “sparsify the tails” in the same spirit as above: remove almost all parameter values in regions where extra roots live, but leave enough tiny islands (like the intervals around rationals) so that there are still parameters with higher likelihood than the surviving root. This lets us enforce (A5) artificially, without doing anything to control the global behavior of the likelihood. Condition (A1), that $$\Theta$$ is an open interval, blocks exactly this trick. 

So the moral is: the geometry of the parameter space matters. Condition (A5) may look strong, but without a topologically simple parameter space like an open interval (A1), it is surprisingly easy to satisfy (A5) and still have an inconsistent MLE.
