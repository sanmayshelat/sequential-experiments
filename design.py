from typing import Tuple

import numpy as np


class Sequential:
	def __init__(self, alpha: float = 0.05, beta: float = 0.2, delta: float = 0.1):
		self.alpha = alpha
		self.beta = beta
		self.delta = delta

	def fixed_sample_proportion(self) -> Tuple[int, int]:
		if self.delta <= 0.02:
			raise NotImplementedError(
				"""
				Effect sizes <= 2%% are currently not supported.
				They can be implemented by searching within a smaller
				space of N and d_star.
				"""
			)
			N_max = 300_000
			d_star_max = 1500
		elif self.delta <= 0.05:
			N_max = 100_000
			d_star_max = 500
		elif self.delta <= 0.1:
			N_max = 12_000
			d_star_max = 250
		else:
			N_max = 3_000
			d_star_max = 120

		log_i_cumsum = np.log(np.arange(1, N_max + 1)).cumsum()
		log_p_h1 = np.log(1 / (2 + self.delta))
		log_q_h1 = np.log((1 + self.delta) / (2 + self.delta))
		log_p_h0 = np.log(1 / 2)

		d_star_range = np.arange(2, d_star_max, 2)
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

		d_star, N = np.unravel_index(
			np.argmax((h0 < self.alpha) & (h1 > (1 - self.beta))), h1.shape
		)
		N = (N + 1) * 2
		d_star = (d_star + 1) * 2

		print(f"""N:{N}, d*: {d_star}""")
		return N, d_star
