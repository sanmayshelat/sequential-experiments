from typing import Annotated, Optional, Sequence, Union

import numpy as np
import pandas as pd
from pydantic import AfterValidator, BaseModel, model_validator


def check_prob_vector(prob):
	prob = np.asarray(prob)
	return (all(prob >= 0)) & (sum(prob) == 1)


class ExperimentImpact(BaseModel):
	perc_change: Sequence[float]
	expt_traffic: float = 1
	expt_prop: Annotated[Sequence[float], AfterValidator(check_prob_vector)] = [
		0.5,
		0.5,
	]

	@model_validator(mode='after')
	def check_same_size(self) -> 'ExperimentImpact':
		if len(self.perc_change) != len(self.expt_prop):
			raise ValueError('perc_change length should be same as expt_prop')
		return self

	@model_validator(mode='after')
	def convert_to_numpy(self) -> 'ExperimentImpact':
		self.perc_change = np.asarray(self.perc_change)
		self.expt_prop = np.asarray(self.expt_prop)
		return self


class DgpProportion:
	def __init__(
		self,
		prob_success: Union[Sequence[float], float],
		cluster_prob: Optional[Sequence[float]] = [1],
		seed: Optional[int] = None,
	):
		if isinstance(prob_success, (float, int)):
			prob_success = [prob_success]
		self.prob_success = np.asarray(prob_success)
		self.rng = np.random.default_rng()

		self.cluster_prob = np.asarray(cluster_prob)
		if not check_prob_vector(self.cluster_prob):
			raise ValueError('cluster_prob is not a probability vector')

		if not (len(self.cluster_prob) == len(prob_success)):
			raise ValueError(
				'cluster_prob and prob_success should be of the same length'
			)

	def dgp(
		self, n_units: int, expt_impact: Optional[ExperimentImpact] = None
	) -> pd.DataFrame:
		unit_clusters = self.rng.choice(
			a=len(self.cluster_prob), size=n_units, p=self.cluster_prob
		)

		unit_success_probs = self.prob_success[unit_clusters]

		if expt_impact is not None:
			unit_treatments = self.rng.choice(
				a=len(expt_impact.expt_prop), size=n_units, p=expt_impact.expt_prop
			)

			unit_success_probs *= 1 + expt_impact.perc_change[unit_treatments] / 100

		outcomes = self.rng.binomial(n=1, p=unit_success_probs)

		df = pd.DataFrame(
			{
				'unit_id': range(n_units),
				'cluster_id': unit_clusters,
				'success_prob': unit_success_probs,
				'outcome': outcomes,
			}
		)

		if expt_impact is not None:
			df = df.assign(treatmeant=unit_treatments)

		return df
