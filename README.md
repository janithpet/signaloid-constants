# signaloid-constants

Implements a method for calculating the constant $e$ using the [Signaloid UxHw API](https://docs.signaloid.io/docs/hardware-api/).


## The theory

### Calcuating $e$ using the Monte Carlo method

The following is a summary of the method present in the paper [Estimating the value of $e$ using the Monte Carlo method by K. G. Russel](https://www.jstor.org/stable/2685243?seq=1#page_scan_tab_contents).

The method is as follows. Let $E$ be a random variable defined as the minimum number $n \in \mathcal{N}$ such that $\sum_{i=1}^n u_i > 1$, where $u_i \sim \mathcal{U}(0, 1)$. $\mathcal{U}(0, 1)$ denotes the uniform distribution on the interval $[0, 1]$. Then, $e$ can be approximated by the expected value $\mathbb{E}[E]$.

## Implementation
Implementations of the raw Monte Carlo estimate and the Signaloid UxHw method can be found in `src/main.c`.

### Raw Monte Carlo method

The code for the raw Monte Carlo method is as follows:
```C
double simple_e(void) {
	double n_monte_carlo = 50000;
	double total_count   = 0;

	for (int i = 0; i < n_monte_carlo; i++) {
		double rolling_sum = 0;
		double count       = 0;

		while (rolling_sum < 1) {
			count += 1;
			rolling_sum += sample_uniform(0, 1);
		}

		total_count += count;
	}

	return total_count / n_monte_carlo;
}
```

We carry out `n_monte_carlo` repetitions of
```C
double rolling_sum = 0;
double count       = 0;

while (rolling_sum < 1) {
	count += 1;
	rolling_sum += sample_uniform(0, 1);
}

total_count += count;
```
where we keep a `count` of the number of samples it took for the `rolling_sum` to exceed 1. During each iteration of the `while` loop, we sample from a uniform using `sample_uniform(0, 1)` and add it to the `rolling_sum`.

The `total_count` keeps a rolling sum of the `count` after each Monte Carlo iteration, and is divided by `n_monte_carlo` at the end to give the average number of samples it took for the `rolling_sum` to exceed 1. This is the expected value of $e$.

## Signaloid UxHw method
The code for this method is below:

```C
double e(void) {
	double rolling_sum = 0;
	double n           = 10;
	double count       = 0;

	for (int i = 0; i < n; i++) {
		double p = UxHwDoubleProbabilityGT(rolling_sum, 1.0);
		rolling_sum += UxHwDoubleUniformDist(0, 1);

		count += UxHwDoubleMixture(1, 0, 1 - p);
	}

	return UxHwDoubleNthMoment(count, 1);
}
```

The key difference here is that we don't need an outer Monte Carlo loop. Instead we deal with the distributions directly. That is, we initialize the `rolling_sum` and `count` to be zero. We then loop for `n` times the following code:
```C
double p 	= UxHwDoubleProbabilityGT(rolling_sum, 1.0);
count 		+= UxHwDoubleMixture(1, 0, 1 - p);

rolling_sum += UxHwDoubleUniformDist(0, 1);
```

The first line calculates the probability that the `rolling_sum` is greate than 1. We can think of this as the proportion of the current distribution of the `rolling_sum` that is greate than 1.

We then add to the `count` a mixture between 1 and 0, where the probability of adding 1 is equal to $1-p$. We can think of this as an if-condition, which adds 1 to the `count` if the `rolling_sum` is *less* than 1, and 0 otherwise. This method simulates the process of keeping count of the number of iterations until the `rolling_sum` exceeds 1.

At the end of the loop, we add a uniform distribution to the `rolling_sum`.

Finally, we return the first moment of the `count`, which is the expected value of `count` using
```C
return UxHwDoubleNthMoment(count, 1);
```

## Notes
We can see that the Signaloid UxHw implementation is much closer to the mathematical description of the process. The only trick that we had to use was the use of a mixture `UxHwDoubleMixture` to simulate an `if-condition`.


