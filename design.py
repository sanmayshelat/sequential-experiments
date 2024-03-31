from typing import Tuple

import numpy as np


class Sequential:
	def __init__(self, alpha: float = 0.05, beta: float = 0.2, delta: float = 0.1):
		self.alpha = alpha
		self.beta = beta
		self.delta = delta

	def fixed_sample_proportion(self) -> Tuple[int, int]:
		# search space
		N_max = 5000
		d_star_max = 200

		log_i_cumsum = np.log(range(1, N_max)).cumsum()
		log_p_h1 = np.log(1 / (2 + self.delta))
		log_q_h1 = np.log((1 + self.delta) / (2 + self.delta))
		log_p_h0 = np.log(1 / 2)

		flag = 0

		for N in range(N_max):
			for d_star in range(2, d_star_max, 2):
				h1 = sum(
					np.exp(
						np.log(d_star / n)
						+ log_i_cumsum[n - 1]
						- log_i_cumsum[max(int((n - d_star) / 2) - 1, 0)]
						- log_i_cumsum[int((n + d_star) / 2) - 1]
						+ ((n - d_star) / 2) * log_p_h1
						+ ((n + d_star) / 2) * log_q_h1
					)
					for n in range(2, N + 1, 2)
				)
				if h1 <= 1 - self.beta:
					break

				h0 = sum(
					np.exp(
						np.log(d_star / n)
						+ log_i_cumsum[n - 1]
						- log_i_cumsum[max(int((n - d_star) / 2) - 1, 0)]
						- log_i_cumsum[int((n + d_star) / 2) - 1]
						+ n * log_p_h0
					)
					for n in range(2, N + 1, 2)
				)
				if (h1 > 1 - self.beta) & (h0 < self.alpha):
					flag = 1
					break

			if flag == 1:
				break
		print(
			f"""
            P(H0):{h0}
            P(H1):{h1}
            N:{N}, d*: {d_star}
        """
		)
		return N, d_star
