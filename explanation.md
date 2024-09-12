# Fixed-wins proportion sequential test
This is the method suggested [here](https://www.evanmiller.org/sequential-ab-testing.html#cite4). Although it is very well explained there, I will add some notes on the derivation of the formula.

## Gambler's ruin
The method is based on the Gambler's ruin and uses the probability of a gambler getting ruined (wealth becoming $\\$0$) within $N$ rounds after starting with an initial wealth of $\\$d$ and playing a game where in each round they have an equal probability of winning or losing $\\$1$.

To derive the formula for this probability, we need the following:

### 1. Number of paths from $\\$0$ to $\\$d$:

This is the number of ways the gambler can get from $\\$0$ to $\\$d$ in exactly $\mathbb{N}$ rounds, assuming they are allowed to have a negative balance. (This is the opposite of what we want but note that the calculation is symmetric and equal in either direction). The gambler needs $d$ wins for sure and the rest ($n-d$) have to be equally distributed betweeen win and loss (i.e., $(n-d)/2$). Thus, $(n+d)/2$ wins and $(n-d)/2$ losses are needed. The number of ways we can arrange $n$ items is $n!$, but given there are identical items, the permutations are given by:

$$
\begin{align*}
N_{n}(0,d) &= \frac{n!}{
	\left( \frac{n-d}{2} \right)!
	\left( \frac{n+d}{2} \right)!
}
= \binom{n}{\frac{n-d}{2}}
\end{align*}
$$

**Note**: if $n-d$ is not even the binomial coefficient is assumed to be 0 (i.e., there can be no paths when $n-d$ is odd). 

### 2. Probability of not revisiting $\\$0$

This is the probability that starting at $\\$0$, the gambler never comes back to this value given their number of wins ($w=(n+d)/2$) and losses ($l=(n-d)/2$). This probability is given by the [ballot theorem](https://en.wikipedia.org/wiki/Bertrand%27s_ballot_theorem):

$$
\begin{align*}
p_{ahead} &= \frac{w-l}{w+l}
= \frac{d}{n}
\end{align*}
$$

### 3. Combining the above:

Let the number of paths going from $\\$0$ to $\\$d$ without revisiting $\\$0$ be $N_{n}^{0}(0,d)$. Then the probability of getting a path that touches $\\$0$ only once ($\pi$) is given by the ratio of these paths and the total number of paths possible as well as the probability given by the ballot theorem above.

$$
\begin{align*}
\pi &= \frac{N_{n}^{0}(0,d)}{N_{n}(0,d)} = \frac{d}{n} \\
\implies N_{n}^{0}(0,d) &= \frac{d}{n}N_{n}(0,d) \\
&= \frac{d}{n}\binom{n}{\frac{n-d}{2}} \\
\end{align*}
$$

To calculate the probability of taking a path that touches $\$0$ only once, we sum the product of probabilities of the required outcome at each step for each path over all the paths. Note that the probability of the outcomes at each step (win or loss) is 0.5. For a path of $n$ steps, the probability of that particular path is $(1/2)^{n}$. We have $N_{n}^{0}(0,d)$ such paths.

$$
\begin{align*}
\implies p_{n}^{0}(0,d) &= \frac{d}{n}\binom{n}{\frac{n-d}{2}}\left(\frac{1}{2}\right)^{n}
\end{align*}
$$

### 4. Reflectivity:

Note the symmetric nature of the problem means that we can also get the probability we originally wanted (i.e., going from $\\$d$ to $\\$0$).

$$
\begin{align*}
\implies p_{n}^{0}(0,d) &= p_{n}^{0}(d,0) = \frac{d}{n}\binom{n}{\frac{n-d}{2}}\left(\frac{1}{2}\right)^{n}
\end{align*}
$$

### 5. Adding up probabilities:

To get the probability of a gambler getting to $\\$d$ within $\mathbb{N}$ rounds, we can add up the exact probabilities, $p_{n}^{0}(0,d)$:

$$
\begin{align*}
\implies p_{\mathbb{N}}(0,d) &= p_{\mathbb{N}}(d,0) = \sum_{n=1}^{\mathbb{N}}\frac{d}{n}\binom{n}{\frac{n-d}{2}}\left(\frac{1}{2}\right)^{n}
\end{align*}
$$

**Reminder**: if $n-d$ is not even the binomial coefficient is assumed to be 0 (i.e., there can be no paths when $n-d$ is odd).


## Analogy between experiment and gambler's ruin
If we consider the difference between treatment and control successes as wealth ($d=\mathbb{T}-\mathbb{C}$), then the wealth increases by one for each treatment success (think of this as the gambler winning) and decreases by one for each control success (think of this as the adversary winning). Note that we have discarded failurs from both treatment and control and the wealth metric makes it like the zero-sum game that we saw in the gambler's ruin problem. 


## Power analysis

### Null hypothesis
Thus, we know the probability of this difference getting to some value, $d^{\*}$, from 0 in $\mathbb{N}$ rounds (i.e., experiment observations after discarding failures). For a given $\mathbb{N}$, we can find the smallest $d^{\*}$ for which this probability is less than $\alpha$ (type I error). Under the null hypothesis, the distribution of $d$ is centred around zero because both treatment and control have equal probability of success. Therefore, under $H_0$:

$$
\begin{align*}
\sum_{n=1}^{\mathbb{N}}\frac{d^\*}{n}\binom{n}{\frac{n-d^\*}{2}}\left(\frac{1}{2}\right)^{n} < \alpha
\end{align*}
$$

**Side note**: it can be proven that $d^{*} \approx z_{\alpha/2}\sqrt{N}$ (proof to be added).

### Alternative hypothesis
Under the alternate hypothesis, the success probabilities for the treatment and control are different. Let, $\frac{p_t}{p_c}=1+\delta$, where $\delta$ indicates the relative lift. Then, if there are total $K$ observations, the fraction of control successes amongst all successes is:

$$
\begin{align*}
p &= \frac{p_c K}{p_c K + p_t K} = \frac{p_c}{p_c+p_t} \\
&= \frac{1}{1+\frac{p_t}{p_c}} = \frac{1}{2+\delta}
\end{align*}
$$

And the fraction of treatment successes amongst all successes is:

$$
\begin{align*}
q &= 1-p = \frac{1+\delta}{2+\delta}
\end{align*}
$$

With different outcome probabilities, the value of $p_{\mathbb{N}}(0,d)$ changes (win and loss probabilities are now included separately):

$$
\begin{align*}
p_{\mathbb{N}}(0,d) &= \sum_{n=1}^{\mathbb{N}}\frac{d}{n}\binom{n}{\frac{n-d}{2}}
p^{(n-d)/2}q^{(n+d)/2}\\
&= \sum_{n=1}^{\mathbb{N}}\frac{d}{n}\binom{n}{\frac{n-d}{2}}
\frac{1}{2+\delta}^{(n-d)/2}\frac{1+\delta}{2+\delta}^{(n+d)/2}
\end{align*}
$$

And we would like this probability to be higher than power ($1-\beta$):

$$
\begin{align*}
\sum_{n=1}^{\mathbb{N}}\frac{d^\*}{n}\binom{n}{\frac{n-d^\*}{2}}
\frac{1}{2+\delta}^{(n-d^\*)/2}\frac{1+\delta}{2+\delta}^{(n+d^\*)/2} &< 1-\beta
\end{align*}
$$

### Calculation
Finally, solving the following to calculate the smallest $\mathbb{N}$ and $d^*$ which will satisfy our experiment requirements will give us a one-sided sequential test:

$$
\begin{align*}
\sum_{n=1}^{\mathbb{N}}\frac{d^\*}{n}\binom{n}{\frac{n-d^\*}{2}}\left(\frac{1}{2}\right)^{n} &< \alpha\\
\sum_{n=1}^{\mathbb{N}}\frac{d^\*}{n}\binom{n}{\frac{n-d^\*}{2}}
\frac{1}{2+\delta}^{(n-d^\*)/2}\frac{1+\delta}{2+\delta}^{(n+d^\*)/2} &< 1-\beta
\end{align*}
$$

## Application
After each observation in the experiment, we calculate $\mathbb{T}-\mathbb{C}$ and $\mathbb{T}+\mathbb{C}$. Our experiment continues until either:

1. $\mathbb{T}-\mathbb{C}=d^\*$: then declare treatment the winner (reject $H_0$)
2. $\mathbb{T}+\mathbb{C}=\mathbb{N}$: then we are unable to reject $H_0$ 
