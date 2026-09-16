# Draft replacement for equation (3.8)

Replace the displayed equation (3.8) and the sentence introducing it with:

> Since $T_n$ and $T_{n-1}^{(-i)}$ are two-sample $U$-statistics with the same
> kernel $h$, evaluated on $(n_1,n_2)$ and on $(n_1-1,n_2)$ or $(n_1,n_2-1)$
> observations respectively, both are unbiased for $D$. Hence
> $$E[V_i] \;=\; n\,D-(n-1)\,D \;=\; D, \qquad i=1,\dots ,n, \tag{3.8}$$
> for every pseudo-value, irrespective of $n_1$, $n_2$ and of which of the two
> samples $W_i$ came from. The constraint in (3.7) is therefore
> $\sum_{i=1}^{n}p_i\,(V_i-\theta)=0$.

LaTeX:

```latex
Since $T_n$ and $T_{n-1}^{(-i)}$ are two-sample $U$-statistics with the same
kernel $h$, evaluated on $(n_1,n_2)$ and on $(n_1-1,n_2)$ or $(n_1,n_2-1)$
observations respectively, both are unbiased for $D$. Hence
\begin{equation}\label{eq:EV}
  E[V_i] \;=\; n\,D-(n-1)\,D \;=\; D, \qquad i=1,\dots ,n,
\end{equation}
for every pseudo-value, irrespective of $n_1$, $n_2$ and of which of the two
samples $W_i$ came from. The constraint in \eqref{eq:EL} is therefore
$\sum_{i=1}^{n}p_i\,(V_i-\theta)=0$.
```

Two knock-on edits:

1. Anywhere the text refers to the "expected pseudo-values" being
   group-specific, or to constants depending on `n_1` and `n_2`, that wording
   goes with the old (3.8).
2. Section 3.1's variance estimator (3.2) uses the *group-specific*
   pseudo-values `V_{i,0} = n_1 U - (n_1-1) U^{(-i)}` and
   `V_{0,j} = n_2 U - (n_2-1) U^{(-j)}`. That is a different and equally
   standard construction from the pooled-`n` pseudo-values of (3.5), and the
   code keeps the two separate. It is worth one sentence saying so, because a
   reader who assumes (3.2) and (3.5) use the same `V` will be confused by the
   different multipliers.

Tables 5 and 6 must be regenerated; the values are in
`Uncensored Simulations/expected_outputs/`.
