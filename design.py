from typing import Tuple

import numpy as np
from scipy.stats import norm


class Sequential:
	def __init__(self, alpha: float = 0.05, beta: float = 0.2, delta: float = 0.1):
		self.alpha = alpha
		self.beta = beta
		self.delta = delta

	def fixed_wins_proportion(self) -> Tuple[int, int]:
		N_min = 2
		if self.delta < 0.01:
			raise NotImplementedError(
				"""
				Effect sizes <= 1%% are currently not supported.
				They can be implemented by searching within a smaller
				space of N.
				"""
			)
			N_min = 250_000
			N_max = 500_000
		elif self.delta <= 0.02:
			N_min = 200_000
			N_max = 300_000
		elif self.delta <= 0.05:
			N_min = 10_000
			N_max = 70_000
		elif self.delta <= 0.1:
			N_min = 2_500
			N_max = 10_000
		else:
			N_max = 3_000

		log_i_cumsum = np.log(np.arange(1, N_max + 1)).cumsum()
		log_p_h1 = np.log(1 / (2 + self.delta))
		log_q_h1 = np.log((1 + self.delta) / (2 + self.delta))
		log_p_h0 = np.log(1 / 2)

		d_star_range = np.unique(
			np.ceil(
				np.abs(norm.ppf(self.alpha / 2))
				* np.sqrt(np.arange(N_min, N_max + 1, 2))
				/ 2
			)
			* 2
		).astype(int)
		d_star_values = np.tile(d_star_range, (N_max // 2, 1)).T
		n_values = np.tile(np.arange(2, N_max + 1, 2), (len(d_star_range), 1))

		n_d_star_2_minus = (n_values - d_star_values) // 2
		n_d_star_2_minus_max = np.maximum(n_d_star_2_minus - 1, 0)

		n_d_star_2_plus = (n_values + d_star_values) // 2
		n_d_star_2_plus = n_d_star_2_plus - 1

		h1_terms = np.exp(
			np.log(d_star_values / n_values)
			+ log_i_cumsum[n_values - 1]
			- log_i_cumsum[n_d_star_2_minus_max]
			- log_i_cumsum[n_d_star_2_plus]
			+ (n_d_star_2_minus) * log_p_h1
			+ (n_d_star_2_plus + 1) * log_q_h1
		)
		h1 = np.cumsum(h1_terms, axis=1)

		h0_terms = np.exp(
			np.log(d_star_values / n_values)
			+ log_i_cumsum[n_values - 1]
			- log_i_cumsum[n_d_star_2_minus_max]
			- log_i_cumsum[n_d_star_2_plus]
			+ n_values * log_p_h0
		)
		h0 = np.cumsum(h0_terms, axis=1)

		d_star_idx, N_idx = np.unravel_index(
			np.argmax((h0 < self.alpha) & (h1 > (1 - self.beta))), h1.shape
		)
		N = n_values[0, :][N_idx]
		d_star = d_star_range[d_star_idx]

		print(f"""N:{N}, d*: {d_star}""")
		return N, d_star
