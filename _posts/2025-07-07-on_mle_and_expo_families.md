---
title: 'Asymptotic Normality of the MLE for Exponential Families'
date: 2025-07-07
permalink: /posts/2025/07/on_mle_and_expo_families/
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


While studying for an exam on statistical inference I became interested in the following question: **when is the maximum‑likelihood estimator (MLE) for the natural parameter of an exponential family asymptotically efficient?** After crawling through forums and textbooks I became certain that the answer was almost always. But I could not find a straightforward theorem or result to use. Or at least not one with a proof. So I set out to prove such a result for myself. The result of my toils is [Theorem 3](#theorem-3) below. Surprisingly, if you ask ChatGPT for a theorem on this topic it will give you one similar to Theorem 3. However, it cannot provide an accurate reference for it. So this theorem likely exists somewhere on the internet, I was just unable to find it.   

The reasoning I provide in this post will not be the most concise. Instead I hope to provide an intuitive story that aligns more closely with my journey to the proof.

## 1. Exponential families


<div class="bbox" markdown="1">

**Definition 1 (Natural exponential family):**

Let $$(\mathcal X,\mathcal F,\mu)$$ be a σ‑finite measure space. A **natural exponential family (NEF)** with $$s$$ parameters is the collection of probability densities

$$
p_{\eta}(x)
= \exp\bigl\{\eta^{\top}T(x) - A(\eta)\bigr\} \, h(x),
\qquad x \in \mathcal X,\; \eta \in \Xi \subseteq \mathbb R^{s},
$$ 

where

* $$T : \mathcal X \to \mathbb R^{s}$$ is a measurable **sufficient statistic**;
* $$h : \mathcal X \to [0,\infty)$$ is a fixed **carrier**;
* $$A(\eta)$$ is the **log‑partition function** defined as

$$
A(\eta) = \log \int_{\mathcal X} \exp\bigl\{\eta^{\top}T(x)\bigr\} h(x)\, d\mu(x).
$$ 

</div>

### Why Should We Care About Exponential Families?

For those unfamiliar with exponential families, the definition above may seem a bit odd. Why, for example, is the log-partition function $$A(\eta)$$ defined as it is? The answer is normalization. The term $$\exp\{ \eta^\top T(x) \} h(x)$$ may not integrate to 1 on its own, so we define

$$
A(\eta) = \log \int_{\mathcal X} \exp\{ \eta^\top T(x) \} h(x)\, d\mu(x)
$$ 

to ensure that the resulting density integrates to one. In this sense, $$A(\eta)$$ acts as the normalizing constant in the density expression. Another question is, what exactly is $$\Xi$$? We define it as

$$
\Xi := \left\{ \eta \in \mathbb{R}^s : |A(\eta)| < \infty \right\},
$$ 

the set of all parameters for which the integral above converges. This is called the **natural parameter space**—the region where $$p_\eta(x)$$ defines a valid probability density.

But why are we interested in such a form in the first place? Because many of the most commonly used statistical distributions can actually be expressed as exponential families. This includes the Normal, Bernoulli, Poisson, and Gamma distributions, as well as many others.


<div class="bbox" markdown="1">

**Example 1:** Let $$X \sim \text{Bern}(\theta)$$, where $$\theta \in (0,1)$$. The density is:

$$
p_\theta(x) = \theta^x (1 - \theta)^{1 - x}, \quad x \in \{0,1\}.
$$ 

Rewriting this:

$$
p_\theta(x) = \exp\left\{ x \log \left( \frac{\theta}{1 - \theta} \right) + \log(1 - \theta) \right\}
= \exp\left\{ \eta x - A(\eta) \right\},
$$ 

where $$\eta = \log \left( \frac{\theta}{1 - \theta} \right)$$ and $$A(\eta) = \log(1 + e^\eta)$$. This is an exponential family with sufficient statistic $$T(x) = x$$.

</div>

<br>
Given the ubiquity of exponential families, it's not uncommon to encounter the task of estimating the natural parameter $$\eta$$. A common choice is to use the maximum likelihood estimator (MLE). But is this a reasonable choice? If we can show that the MLE is asymptotically efficient for these families, we gain a theoretical justification for the use of the MLE.


## 2. The Proof Idea

Let's try to prove that the MLE in an exponential family is asymptotically efficient. To do so we'll take a "naive" approach. We'll assume many things out of convenience. Pay attention carefully for these hidden assumptions. Our naive approach will be far from rigorous, but it will serve as intuition for a more careful analysis.

Let $$X_1, \dots, X_n \overset{\text{iid}}{\sim} p_{\eta_0}$$, where $$p_{\eta}(x) = \exp\bigl\{ \eta T(x) - A(\eta) \bigr\} h(x)$$ defines a one-parameter natural exponential family (the argument we present can be easily extended to the multi-parameter case). The log-likelihood function is, up to an additive constant,

$$
\ell_n(\eta)
= \sum_{i=1}^n \log p_{\eta}(X_i)
\propto_{\eta} \eta \sum_{i=1}^n T(X_i) - n A(\eta)
= n \left( \eta \overline{T}_n - A(\eta) \right),
$$ 

where $$\overline{T}_n = \frac{1}{n} \sum_{i=1}^n T(X_i)$$ is the empirical mean of the sufficient statistic. To compute the MLE we must find a maximizer of the above equation. Taking the derivative and setting it equal to zero yields the likelihood equation:

$$
\frac{d}{d\eta} \ell_n(\eta)
= n \left( \overline{T}_n - A'(\eta) \right) = 0
\quad \Longrightarrow \quad
A'(\eta) = \overline{T}_n.
$$ 

Let us define the map $$\psi := A'.$$ Assuming $$\psi$$ is invertible, we can write the MLE as

$$
\hat{\eta}_n = \psi^{-1}(\overline{T}_n).
$$ 

Now apply the central limit theorem to $$\overline{T}_n$$:

$$
\sqrt{n} \left( \overline{T}_n - \mathbb{E}_{\eta_0}[T(X)] \right)
\;\xrightarrow{d}\;
\mathcal{N}\left( 0, \operatorname{Var}_{\eta_0}[T(X)] \right).
$$ 

Then, by the delta method, since $$\hat{\eta}_n = \psi^{-1}(\overline{T}_n)$$, we obtain

$$
\sqrt{n} \left( \hat{\eta}_n - \psi^{-1}(\mathbb{E}_{\eta_0}[T(X)]) \right)
\;\xrightarrow{d}\;
\mathcal{N} \left( 0, \left( \psi'(\eta_0) \right)^{-2} \operatorname{Var}_{\eta_0}[T(X)] \right).
$$ 

This already gives us the asymptotic normality of the MLE. Furthermore, we will soon show that

$$
\psi^{-1}(\mathbb{E}_{\eta_0}[T(X)]) = \eta_0
\quad \text{and} \quad
\left( \psi'(\eta_0) \right)^{-2} \operatorname{Var}_{\eta_0}[T(X)] = I^{-1}(\eta_0),
$$ 

where $$I(\eta_0)$$ is the Fisher information at $$\eta_0$$. This implies the MLE is asymptotically efficient, since

$$
\boxed{
\sqrt{n}(\hat{\eta}_n - \eta_0)
\;\xrightarrow{d}\;
\mathcal{N}\left( 0,\, I^{-1}(\eta_0) \right).
}
$$ 

Why isn't this is a rigorous proof? The answer is that there are several hidden assumptions we have used including:

* Any critical point of $$\ell_n(\eta)$$ is a global maximizer.
* $$A(\eta)$$ is twice differentiable.
* $$\psi(\eta) = A'(\eta)$$ is invertible and its inverse is differentiable.
* $$\mathbb{E}_{\eta_0}[T(X)] = \psi(\eta_0)$$ and $$\operatorname{Var}_{\eta_0}[T(X)] = \left( \psi'(\eta_0) \right)^{2} I^{-1}(\eta_0)$$.

These are all subtle points we will need to justify. If you do not see where we used these assumptions I highly encourage you to reread the proof.

## 3. Towards the Main Proof

### Smoothness of $$A(\eta)$$:

Let's begin by addressing the smoothness of $$A(\eta)$$. We must show that it is twice differentiable. In fact we will show a much stronger result:

<div class="bbox" markdown="1">

**Theorem 1:** 
Let $$A(\eta)$$ be the log-partition function of a natural exponential family. Then $$A$$ is infinitely differentiable on the interior $$\Xi^\circ$$ of the natural parameter space.

</div>
<br>
To see why this might be true take a look at the definition of $$A(\eta)$$. It looks suspiciously like the logarithm of an MGF. For MGFs we have the following standard theorem. 
<div class="bbox" markdown="1">

**Theorem 2:** 
Let $$X$$ be a random variable taking values in $$\mathbb R^s$$. If the moment generating function (MGF) of $$X$$ defined as $$M_X(z) := \mathbb{E}[\exp(z^\top X)]$$ is finite in an open neighborhood of $$0$$ then $$M_X(z)$$ is infinitely differentiable in an open neighborhood of $$0$$ and $$X$$ has finite moments of all orders.
</div>
<br>
[A proof of this result can be found in Theorem 1.1 of this reference.](https://researchers.ms.unimelb.edu.au/~xgge@unimelb/Files/Notes/The%20Mathematical%20Theory%20of%20Moment%20Generating%20Functions.pdf)
Unfortunately, $$\exp(A(\eta))$$ does not have the exact form of an MGF. The problem stems from the fact that we are not integrating against a probability measure. However with a little bit of work we can get around this.

<div class="bbox" markdown="1">
*Proof of Theorem 1:*
<br>

Fix any point $$\eta_0 \in \Xi^\circ$$. Let $$\phi(\eta) := \int_{\mathcal{X}} \exp\bigl\{ \eta^\top T(x) \bigr\} h(x) d\mu(x)$$ be the unnormalized integral, so that $$A(\eta) = \log \phi (\eta)$$. We know that for all $$\eta \in \Xi^\circ$$, $$0 < \phi(\eta) < \infty$$ (this follows from the definition of $$\Xi$$). Then we can define a probability measure $$Q_{\eta_0}$$ on $$(\mathcal{X}, \mathcal{F})$$ via

$$
dQ_{\eta_0}(x) = \frac{e^{\eta_0^\top T(x)} h(x)}{\phi(\eta_0)}\, d\mu(x).
$$ 

This is a valid probability measure because the denominator normalizes the density. For $$\eta$$ near $$\eta_0$$, define $$\delta := \eta - \eta_0$$. Then:

$$
\phi(\eta)
= \int e^{(\eta_0 + \delta)^\top T(x)} h(x)\, d\mu(x)
= \phi(\eta_0) \cdot \mathbb{E}_{Q_{\eta_0}}\left[ e^{\delta^\top T(X)} \right]
= \phi(\eta_0) \cdot M_Z(\delta),
$$ 

where $$Z := T(X)$$ and $$M_Z(\delta)$$ is the moment-generating function of $$Z$$ under $$Q_{\eta_0}$$. Now we are in position to apply Theorem 2. Since $$\eta_0 \in \Xi^\circ$$, and $$\Xi^\circ$$ is open, there exists a neighborhood of $$0$$ around $$\delta = 0$$ on which $$M_Z(\delta) = \phi(\eta_0 + \delta) < \infty$$. Then by Theorem 2, $$M_Z$$ is infinitely differentiable at $$\eta_0$$. Since

$$
A(\eta) = \log \phi(\eta) = \log \phi(\eta_0) + \log M_Z(\delta),
$$ 

which is the sum of two $$C^\infty$$ functions in a neighborhood of $$\eta_0$$ we can conclude that $$A$$ is $$C^\infty$$ in this neighborhood. Since $$\eta_0 \in \Xi^\circ$$ was arbitrary, we conclude $$A \in C^\infty(\Xi^\circ).$$ 

</div>

### Moments of $$T(X)$$:
Next, let's address our claims about the moments of $$T(X)$$. To do this, we'll fix some $$\eta_0 \in \Xi^\circ$$ and we'll compute the MGF of $$T(X)$$ under the distribution $$P_{\eta_0}$$:

$$
\begin{aligned}
M_T(z) 
&= \int e^{z^\top T(x)}\, p_{\eta_0}(x)\, d\mu(x) \\
&= \int \exp\left( z^\top T(x) + \eta_0^\top T(x) - A(\eta_0) \right) h(x)\, d\mu(x) \\
&= \int \exp\left( (\eta_0 + z)^\top T(x) \right) h(x)\, d\mu(x)\, \cdot\, \exp\big(- A(\eta_0) \big) \\
&= \exp\big( A(\eta_0 + z) - A(\eta_0) \big).
\end{aligned}
$$

This expression for the moment-generating function,

$$
M_T(z) = \exp\big( A(\eta_0 + z) - A(\eta_0) \big),
$$

is finite in a neighborhood of $$z = 0$$. This follows immediately from the fact that $$\eta_0 + z \in \Xi^\circ$$ for small enough $$z$$, and $$A(\eta_0) < \infty$$ throughout $$\Xi^\circ$$ by definition. Therefore, $$M_T(z) < \infty$$ in some open neighborhood of $$0$$, which guarantees that all moments of $$T(X)$$ exist and are finite under $$P_{\eta_0}$$.

Now, since $$M_T(z) = \exp(A(\eta_0 + z) - A(\eta_0))$$, we can compute the first and second moments by differentiation:

$$
\begin{aligned}
\mathbb{E}_{\eta_0}[T(X)]  
&= \left. \nabla_z \log M_T(z) \right|_{z=0}
= \left. \nabla_z A(\eta_0 + z) \right|_{z=0}
= \nabla_\eta A(\eta) |_{\eta=\eta_0},
\end{aligned}
$$

which shows that the gradient of the log-partition function is the mean of the sufficient statistic. Taking a second derivative gives:

$$
\operatorname{Var}_{\eta_0}[T(X)]
= \left. \nabla_z^2 \log M_T(z) \right|_{z=0}
= \nabla^2_\eta A(\eta)|_{\eta=\eta_0},
$$

Thus, the log-partition function $$A(\eta)$$ simultaneously encodes the first and second moments of the sufficient statistic through its gradient and Hessian:

$$
\nabla_\eta A(\eta) |_{\eta=\eta_0} = \mathbb{E}_{\eta_0}[T(X)], \qquad
\nabla^2_\eta A(\eta)|_{\eta=\eta_0} = \operatorname{Var}_{\eta_0}[T(X)].
$$

<div class="bbox" markdown="1">

**Remark:**
There is an alternative way to prove the smoothness of the log-partition function $$A(\eta)$$ using the MGF derived above. Recall that under $$P_\eta$$, we showed:

$$
M_T(z) := \mathbb{E}_\eta[e^{z^\top T(X)}] = \exp\big(A(\eta + z) - A(\eta)\big),
$$

which is finite in a neighborhood of $$z = 0$$ for all $$\eta \in \Xi^\circ$$. Now fix any $$\eta_0 \in \Xi^\circ$$. Since $$M_T(z)$$ is finite in a neighborhood of zero, Theorem 2 tells us that $$M_T(z)$$ is infinitely differentiable in that neighborhood, and hence so is $$\log M_T(z) = A(\eta_0 + z) - A(\eta_0)$$. Therefore, the map $$z \mapsto A(\eta_0 + z)$$ is smooth in a neighborhood of zero, which implies that $$A(\eta)$$ is infinitely differentiable in a neighborhood of $$\eta_0$$.

</div>

### Minimal Exponential Families:

The next point we’ll address is the assumption that any critical point of the log-likelihood function is a global maximizer. This is crucial for the MLE to be well-behaved, and we’ll now justify it in the setting of exponential families.

Recall from before that the log-likelihood function, up to an additive constant, is given by

$$
\ell_n(\eta) \propto_\eta n \left( \eta^\top \overline{T}_n - A(\eta) \right),
$$

where $$\overline{T}_n$$ is the empirical mean of the sufficient statistic. We already know from the previous section that $$A$$ is twice differentiable, and that

$$
\nabla^2_\eta A(\eta) = \operatorname{Var}_\eta[T(X)].
$$

This matrix is always positive semidefinite since it is a variance/covariance matrix. Hence, $$A(\eta)$$ is convex on $$\Xi^\circ$$. However, to guarantee that the log-likelihood has a unique global maximizer, we need a stronger property: strict convexity of $$A(\eta)$$. In particular, we want:

$$
v^\top \nabla^2_\eta A(\eta)\, v = \operatorname{Var}_\eta\left( v^\top T(X) \right) > 0
\quad \text{for all nonzero } v.
$$

That is, the variance of every nontrivial linear combination of the sufficient statistic must be strictly positive under $$P_\eta$$. Equivalently, for all nonzero vectors $$v \in \mathbb{R}^s$$, the random variable $$v^\top T(X)$$ is **not almost surely constant** under any $$P_\eta$$. An exponential family that satisfies this is called **minimal**.

<div class="bbox" markdown="1">

**Definition 2 (Minimal Exponential Family):**

A natural exponential family with sufficient statistic $$T : \mathcal{X} \to \mathbb{R}^s$$ is said to be **minimal** if there does not exist a nonzero vector $$v \in \mathbb{R}^s$$, a scalar $$c \in \mathbb{R}$$, and a parameter vector $$\eta \in \Xi$$ such that

$$
v^\top T(x) = c \quad \text{almost surely under } P_\eta.
$$

Equivalently, $$T(X)$$ does not satisfy any nontrivial affine constraint almost surely under any $$P_\eta$$.

</div>
<br>
Minimality is not a restrictive constraint in practice. In fact, when a family fails to be minimal, it means that the parameterization is redundant—some components of the sufficient statistic are affine combinations of others and carry no additional information. In such cases, the exponential family can be reparameterized into a lower-dimensional minimal form without loss of generality. So, assuming minimality simply ensures that we're working with the "simplest" or most efficient representation of the family.

Thus, if a natural exponential family is minimal then the log-partition function $$A(\eta)$$ is strictly convex on $$\Xi^\circ$$. This implies that the log-likelihood $$\ell_n(\eta) \propto_\eta n(\eta^\top \overline{T}_n - A(\eta))$$ is strictly concave on $$\Xi^\circ$$ and therefore that any critical point of $$\ell_n(\eta)$$ is a unique global maximizer on $$\Xi^\circ.$$

We can actually strengthen the conclusion slightly. Specifically, any critical point of $$\ell_n(\eta)$$ in the interior $$\Xi^\circ$$ is the unique global maximizer of the log-likelihood over the entire parameter space $$\Xi$$ (not just over $$\Xi^\circ$$). This follows from the fact that the log-partition function $$A(\eta)$$ is strictly convex on all of $$\Xi$$—not just on its interior for a minimal family. A proof of this fact can be found in [Theorem 1 of these lecture notes by Michael Jordan](https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/other-readings/chapter8.pdf).


## 4. The Main Proof

We now have almost all pieces we need to prove the following result.

<div id="theorem-3" class="bbox" markdown="1">

**Theorem 3:**
Let $$\{P_\eta : \eta \in \Xi\}$$ be a minimal exponential family, and suppose $$\eta_0 \in \Xi^\circ$$. Let $$X_1, \dots, X_n \overset{\text{iid}}{\sim} P_{\eta_0}$$, and define the maximum-likelihood estimator $$\hat{\eta}_n$$ as any element of

$$
\hat{\eta}_n \in \operatorname{argmax}_{\eta \in \Xi} \; \ell_n(\eta),
$$

with $$\hat{\eta}_n$$ defined arbitrarily when the argmax is empty. Then, the MLE is asymptotically efficient:

$$
\sqrt{n}(\hat{\eta}_n - \eta_0)
\;\xrightarrow{d}\;
\mathcal{N}\left(0,\, I^{-1}(\eta_0)\right).
$$

</div>
<br>
The last theorem we'll need is the inverse function theorem.

<div class="bbox" markdown="1">

**Theorem 4 (Inverse Function Theorem).**  
Let $$f:\mathbb{R}^s\to\mathbb{R}^s$$ be continuously differentiable on an open neighborhood of a point $$\theta_0$$. If the Jacobian matrix of $$f$$ at $$\theta_0$$ (denoted by $$J_f(\theta_0)$$) is non-singular, then there exists open neighborhoods $$\theta_0 \in \mathcal{O}_1$$ and $$f(\theta_0) \in \mathcal{O}_2$$ such that the restricted function $$f:\mathcal{O}_1\to\mathcal{O}_2$$ is continuously differentiable and invertible and its inverse $$f^{-1}:\mathcal{O}_2\to\mathcal{O}_1$$ is also continuously differentiable.

</div>
<br>
Using this we can prove Theorem 3.

<div class="bbox" markdown="1">

*Proof of Theorem 3.*

Proceed exactly as in our naive proof attempt and derive the likelihood equation  

$$
\nabla_\eta A(\eta)=\overline T_n. 
$$

This time we know that $$\nabla_\eta A(\eta)$$ exists for all $$\eta \in \Xi^\circ$$ and that any solution to the likelihood equation in $$\Xi^\circ$$ is the unique maximiser of the log-likelihood $$\ell_n(\eta)$$. Let $$\psi(\eta) := \nabla_\eta A(\eta).$$ Minimality implies that the Jacobian matrix $$J_\psi(\eta_0) = \nabla^2_\eta A(\eta)|_{\eta=\eta_0} = \operatorname{Var}_{\eta_0}[T(X)]$$ is positive-definite, so $$J_\psi(\eta_0)$$ is non-singular. Then by the inverse function theorem there are open neighborhoods $$\eta_0 \in \mathcal{O}_1 \subset\Xi^\circ$$ and $$\psi(\eta_0) \in \mathcal{O}_2$$ such that  $$\psi:\mathcal{O}_1\to\mathcal{O}_2$$ is invertible and its inverse $$\psi^{-1}:\mathcal{O}_2\to\mathcal{O}_1$$ is also continuously differentiable.
<br>
<br>

Then for $$\overline T_n\in \mathcal{O}_2$$ we have that  

$$
\hat\eta_n = \psi^{-1}(\overline T_n).
$$

From the multivariate CLT we have that

$$
\sqrt{n}\bigl(\overline T_n-\psi(\eta_0)\bigr)\xrightarrow{d}\mathcal{N}\!\bigl(0,\;J_\psi(\eta_0)\bigr).
$$

Since $$\psi^{-1}$$ is differentiable at $$\psi(\eta_0)$$ and we can apply the delta method. Moreover, from multivariate calculus we have that $$J_{\psi^{-1}}(\psi(\eta_0))=(J_\psi(\eta_0))^{-1}$$. Combining these two facts implies that

$$
\sqrt{n}(\psi^{-1}(\overline T_n)-\eta_0)\xrightarrow{d}\mathcal{N}\!\bigl(0,\; (J_\psi(\eta_0))^{-1}\bigr).
$$

Since the Fisher information matrix is $$I(\eta_0) = \mathbb{E}[-\nabla_\eta^2 \log p_{\eta}(x)|_{\eta=\eta_0}] = \mathbb{E}[\nabla_\eta^2 A(\eta)|_{\eta=\eta_0}] = J_\psi(\eta_0)$$ we are nearly done. This would be enough to finish the proof if we always had $$\overline T_n \in \mathcal{O}_2$$. Instead we have something slightly weaker. Since $$\mathcal{O}_2$$ is an open neighborhood of $$\psi(\eta_0) = \mathbb{E}_{\eta_0}[T(X)]$$, the law of large numbers implies that $$\mathbb{P}(\overline T_n\in \mathcal{O}_2)\to1$$. 
<br>
<br>

To use this to finish the proof we let $$\mathbf{1}[A]$$ denote the indicator function of the event $$A$$. Since $$\mathbb{P}(\overline T_n\in \mathcal{O}_2)\to1$$ we have that,

$$
\mathbf{1}[\overline T_n\in \mathcal{O}_2] \xrightarrow{p} 1

\,\,\text{ and }\,\,

\sqrt{n}(\hat{\eta}_n-\eta_0) \mathbf{1}[\overline T_n \not\in \mathcal{O}_2] \xrightarrow{p} 0.
$$

Then the claim follows by applying Slutsky's theorem to the following decomposition:

$$
\sqrt{n}(\hat{\eta}_n-\eta_0) = 
\sqrt{n}(\psi^{-1}(\overline T_n)-\eta_0) \mathbf{1}[\overline T_n\in \mathcal{O}_2] + \sqrt{n}(\hat{\eta}_n-\eta_0) \mathbf{1}[\overline T_n \not\in \mathcal{O}_2].
$$

</div>








