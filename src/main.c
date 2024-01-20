#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <uxhw.h>

/* \brief Samples from a uniform distribution.
 *
 * \param min The lower bound of the uniform distribution.
 * \param max The upper bound of the uniform distribution.
 * \return A sample from U(min, max).
 */
double sample_uniform(double min, double max) {
	double s = (double)rand() / RAND_MAX;
	return (s * (max - min)) + min;
}

/* \brief Estimates `e` using the Monte Carlo method.
 *
 * The Monte Carlo method to estimate `e` is to sample from a Uniform
 * distrubtion repeatedly and maintain a rolling sum until the sum is greater
 * than 1, keeping track of the number of iterations taken to become greater
 * than 1.
 *
 * Repeat this, and the expected number of iterations taken to become greater
 * than 1 approaches `e` as the number of repetitions approaces infinity.
 *
 * \return An estimate of `e`.
 */
double simple_e(void) {
	double n_monte_carlo = 50000; // number of Monte Carlo repetitions.
	double total_count   = 0;     // total count for taking average at end.

	for (int i = 0; i < n_monte_carlo; i++) {
		double rolling_sum = 0;
		double count       = 0; // double because we need to take an average.

		while (rolling_sum < 1) {
			count += 1;
			rolling_sum += sample_uniform(0, 1);
		}

		total_count += count;
	}

	return total_count / n_monte_carlo;
}

/*
 * \brief Estimates `e` using a Laplace implementation of the Monte Carlo
 * method.
 *
 * The Monte Carlo method to estimate `e` is to sample from a Uniform
 * distrubtion repeatedly and maintain a rolling sum until the sum is greater
 * than 1, keeping track of the number of iterations taken to become greater
 * than 1.
 *
 * Repeat this, and the expected number of iterations taken to become greater
 * than 1 approaches `e` as the number of repetitions approaces infinity.
 *
 * This function uses the UxHw API to implement the process of summing samples
 * from a Uniform random variable, and eliminates the outer loop (see the
 * simple_e function).
 *
 * \return An estimate of `e`.
 */
double e(void) {
	double rolling_sum = 0;
	double n           = 10; // Converges around 8.
	double count       = 0;

	for (int i = 0; i < n; i++) {
		double p = UxHwDoubleProbabilityGT(rolling_sum, 1.0); // calculates the proportion of rolling_sum that is > 1.
		count += UxHwDoubleMixture(1, 0, 1 - p);              // add to count when rolling_sum < 1.

		rolling_sum += UxHwDoubleUniformDist(0, 1); // updates rolling sum by adding another uniform random variable.
	}

	return UxHwDoubleNthMoment(count, 1);
}

int main(void) {
	srand(time(NULL)); // set random seed.

	printf("Laplace e:\t%f\n", e());
	printf("Simple e:\t%f\0", simple_e());
}
