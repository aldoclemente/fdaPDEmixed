#ifndef __REGRESSION_DATA_IMP_H__
#define __REGRESSION_DATA_IMP_H__

// costructor with WeightMatrix
template <typename MatrixType>
RegressionData<MatrixType>::RegressionData(Real *locations, UInt n_locations, UInt ndim, VectorXr &observations, UInt order, MatrixXr &covariates,
										   VectorXr &WeightsMatrix, std::vector<UInt> &bc_indices, std::vector<Real> &bc_values, MatrixXi &incidenceMatrix, bool arealDataAvg, UInt search) : locations_(locations, n_locations, ndim), observations_(observations), arealDataAvg_(arealDataAvg), WeightsMatrix_(WeightsMatrix),
																																															  order_(order), bc_values_(bc_values), bc_indices_(bc_indices), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
																																															  flag_SpaceTime_(false), search_(search)
{
	nRegions_ = incidenceMatrix_.rows();
	if (locations_.nrows() == 0 && nRegions_ == 0)
	{
		locations_by_nodes_ = true;
		for (UInt i = 0; i < observations_.size(); ++i)
			observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

template <typename MatrixType>
RegressionData<MatrixType>::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP Rcovariates,
										   SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch) : locations_(Rlocations)
{
	flag_SpaceTime_ = false;

	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);

	order_ = INTEGER(Rorder)[0];
	search_ = INTEGER(Rsearch)[0];

	UInt length_indexes = Rf_length(RBCIndices);
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) + length_indexes);
	bc_values_.assign(REAL(RBCValues), REAL(RBCValues) + Rf_length(RBCIndices));

	arealDataAvg_ = INTEGER(RarealDataAvg)[0];
}

// mixed
template <typename MatrixType>
RegressionData<MatrixType>::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP RnumUnits, SEXP RRandomEffect, SEXP Rorder, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, 
										   SEXP Rsearch, SEXP RFLAG_ITERATIVE, SEXP Rthreshold, SEXP Rmax_num_iteration, SEXP Rthreshold_residual, SEXP verbose, SEXP anderson_memory) : locations_(Rlocations)
{
	flag_SpaceTime_ = false;
	flag_Mixed_ = true;

	num_units_ = INTEGER(RnumUnits)[0];
	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservationsTime(Robservations);
	if constexpr (std::is_same<SpMat, MatrixType>::value)
	{
		setCovariates(Rcovariates, RRandomEffect);
	}
	else // fixed effects case
	{
		setCovariates(Rcovariates);
	}

	order_ = INTEGER(Rorder)[0];
	search_ = INTEGER(Rsearch)[0];
	flag_iterative_ = INTEGER(RFLAG_ITERATIVE)[0];
	verbose_ = INTEGER(verbose)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ = REAL(Rthreshold)[0];
	threshold_residual = REAL(Rthreshold_residual)[0];
	UInt length_indexes = Rf_length(RBCIndices);
	anderson_memory_ = INTEGER(anderson_memory)[0];
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) + length_indexes);

	bc_values_.assign(REAL(RBCValues), REAL(RBCValues) + Rf_length(RBCIndices));

	arealDataAvg_ = INTEGER(RarealDataAvg)[0];
}

template <typename MatrixType>
RegressionData<MatrixType>::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Rflag_iterative, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Ric, SEXP Rsearch) : locations_(Rlocations)
{
	flag_SpaceTime_ = true;

	setTimeLocations(Rtime_locations);
	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservationsTime(Robservations);
	setCovariates(Rcovariates);

	order_ = INTEGER(Rorder)[0];
	search_ = INTEGER(Rsearch)[0];
	flag_mass_ = INTEGER(Rflag_mass)[0];
	flag_parabolic_ = INTEGER(Rflag_parabolic)[0];
	flag_iterative_ = INTEGER(Rflag_iterative)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ = REAL(Rthreshold)[0];

	UInt length_indexes = Rf_length(RBCIndices);
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) + length_indexes);
	bc_values_.assign(REAL(RBCValues), REAL(RBCValues) + Rf_length(RBCIndices));

	arealDataAvg_ = INTEGER(RarealDataAvg)[0];

	UInt length_ic = Rf_length(Ric);
	ic_.resize(length_ic);
	for (UInt i = 0; i < length_ic; ++i)
	{
		ic_(i) = REAL(Ric)[i];
	}
}

RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
											   SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
											   SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch) : RegressionData<>(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates,
																														   RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch),
																										  K_(RK), beta_(Rbeta), c_(REAL(Rc)[0]) {}

RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder,
											   SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
											   SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Rflag_iterative, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Ric, SEXP Rsearch) : RegressionData<>(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, Rcovariates,
																																																											RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration, Rthreshold, Ric, Rsearch),
																																																						   K_(RK), beta_(Rbeta), c_(REAL(Rc)[0]) {}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
																	   SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
																	   SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch) : RegressionData<>(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch),
																																  K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{
	;
}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder,
																	   SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg,
																	   SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Rflag_iterative, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Ric, SEXP Rsearch) : RegressionData<>(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, Rcovariates,
																																																						 RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration, Rthreshold, Ric, Rsearch),
																																																		K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{
	;
}

template <typename MatrixType>
void RegressionData<MatrixType>::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if (locations_.nrows() == 0 && nRegions_ == 0)
	{
		locations_by_nodes_ = true;
		for (auto i = 0; i < n_obs_; ++i)
		{
			if (!ISNA(REAL(Robservations)[i]))
			{
				observations_[count] = REAL(Robservations)[i];
				count++;
				observations_indices_.push_back(i);
			}
		}
		observations_.conservativeResize(count, Eigen::NoChange);
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		locations_by_nodes_ = false;
		for (auto i = 0; i < n_obs_; ++i)
		{
			observations_[i] = REAL(Robservations)[i];
		}
	}

	// std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	// for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

template <typename MatrixType>
void RegressionData<MatrixType>::setObservationsTime(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);
	if (isSpaceTime())
		observations_na_.resize(time_locations_.size());
	else if (isMixed())
		observations_na_.resize(num_units_);

	UInt count = 0;
	locations_by_nodes_ = (locations_.nrows() == 0 && nRegions_ == 0);

	for (int i = 0; i < n_obs_; ++i)
	{
		if (!ISNA(REAL(Robservations)[i]))
		{
			observations_(i) = REAL(Robservations)[i];
			observations_indices_.push_back(i);
		}
		else
		{
			observations_(i) = 0.0;
			observations_na_[i / getNumberofSpaceObservations()].push_back(i % getNumberofSpaceObservations());
		}
	}

	bool empty = true;
	for (int i = 0; i < observations_na_.size(); ++i)
	{
		if (!observations_na_[i].empty())
		{
			empty = false;
			break;
		}
	}
	if (empty)
	{
		observations_na_.clear();
		observations_na_.shrink_to_fit();
	}
}

template <typename MatrixType>
void RegressionData<MatrixType>::setCovariates(SEXP Rcovariates)
{
	int n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	int p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];
	UInt k = 0;

	covariates_.resize(n_, p_);

	for (size_t j = 0; j < p_; ++j)
	{
		for (size_t i = 0; i < n_; ++i)
		{
			covariates_(i, j) = REAL(Rcovariates)[i + n_ * j];
		}
	}

	for (size_t j = 0; j < covariates_.cols(); ++j)
	{
		for (size_t k = 0; k < observations_na_.size(); ++k)
		{
			for (size_t i = 0; i < observations_na_[k].size(); i++)
			{
				covariates_(observations_na_[k][i] + getNumberofSpaceObservations() * k, j) = 0;
			}
		}
	}

	if (covariates_.cols() != 1)
	{
		MatrixXr tmp(covariates_.cols(), covariates_.cols());
		tmp.selfadjointView<Eigen::Upper>().rankUpdate(covariates_.transpose());
		dec_mat_type WTW_;
		WTW_.compute(tmp.selfadjointView<Eigen::Upper>());
		tmp.setIdentity();
		WTW_inv = WTW_.solve(tmp);
	}
	else
	{
		WTW_inv.resize(1, 1);
		WTW_inv(0, 0) = 1. / (covariates_.transpose() * covariates_).coeff(0, 0);
	}
}

template <>
inline void RegressionData<SpMat>::setCovariates(SEXP Rcovariates, SEXP RRandomEffect)
{
	int N_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	int q_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];
	int p_ = Rf_length(RRandomEffect);
	int nlocations = N_ / num_units_;

	covariates_.resize(N_, num_units_ * p_ + q_ - p_);
	std::vector<coeff> tripletAll;
	tripletAll.reserve(q_ * N_);
	std::vector<size_t> complement;
	for (size_t i = 0, index = 0; i != q_; ++i) // I exploit the fact that the random effect vector has already been sorted
	{
		if (i != INTEGER(RRandomEffect)[index])
			complement.push_back(i);
		else if (index != p_ - 1)
			++index;
	}

	for (size_t j = 0; j < complement.size(); ++j)
	{
		for (size_t i = 0; i < N_; ++i)
		{
			tripletAll.push_back(coeff(i, j, REAL(Rcovariates)[i + N_ * complement[j]]));
		}
	}

	for (size_t j = 0; j < p_; j++)
	{
		for (size_t i = 0; i < N_; i++)
		{
			tripletAll.push_back(coeff(i, q_ - p_ + j + p_ * (i / nlocations), REAL(Rcovariates)[i + N_ * (INTEGER(RRandomEffect)[j])]));
		}
	}
	covariates_.setFromTriplets(tripletAll.begin(), tripletAll.end());

	for (size_t j = 0; j < covariates_.cols(); ++j)
	{
		for (size_t k = 0; k < observations_na_.size(); ++k)
		{
			for (size_t i = 0; i < observations_na_[k].size(); i++)
			{
				if (covariates_.coeff(observations_na_[k][i] + getNumberofSpaceObservations() * k, j) != 0)
					covariates_.coeffRef(observations_na_[k][i] + getNumberofSpaceObservations() * k, j) = 0;
			}
		}
	}
	covariates_.prune(0.);
	covariates_.makeCompressed();
	if (verbose_)
		Rprintf("Non zeros of design matrix X: %i\n", covariates_.nonZeros());
	SpMat tmp(covariates_.cols(), covariates_.cols());
	tmp.selfadjointView<Eigen::Upper>().rankUpdate(covariates_.transpose());
	dec_mat_type WTW_;
	WTW_.compute(tmp.selfadjointView<Eigen::Upper>());
	if (verbose_)
		Rprintf("Non zeros of L: %i\n", static_cast<SpMat>(WTW_.matrixL()).nonZeros());
	if (verbose_)
		Rprintf("Non zeros of U: %i\n", static_cast<SpMat>(WTW_.matrixU()).nonZeros());
	if (verbose_)
		Rprintf("Size of P (should be like the dimension of WTW ? = %i): %i\n", static_cast<SpMat>(WTW_.matrixU()).cols(), WTW_.permutationP().size());
	tmp.setIdentity();
	WTW_inv = WTW_.solve(tmp);
}

template <typename MatrixType>
void RegressionData<MatrixType>::setTimeLocations(SEXP Rtime_locations)
{
	UInt n_time_loc_ = Rf_length(Rtime_locations);
	time_locations_.resize(n_time_loc_);

	for (auto i = 0; i < n_time_loc_; ++i)
	{
		time_locations_[i] = REAL(Rtime_locations)[i];
	}
}

template <typename MatrixType>
void RegressionData<MatrixType>::setBaryLocations(SEXP RbaryLocations)
{
	// RECIEVE BARYCENTER INFORMATION FROM R
	if (TYPEOF(RbaryLocations) != 0)
	{ // TYPEOF(RbaryLocations) == 0 means SEXPTYPE is NILSXP (Description is NULL)
		UInt *id_ = INTEGER(VECTOR_ELT(RbaryLocations, 1));
		Real *bary_ = REAL(VECTOR_ELT(RbaryLocations, 2));

		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[0];
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[1]; // barycenter column dimension

		barycenters_.resize(n_, p_);
		element_ids_.resize(n_);

		if (n_ > 0)
		{
			for (auto i = 0; i < n_; ++i)
			{
				for (auto j = 0; j < p_; ++j)
				{
					barycenters_(i, j) = bary_[i + n_ * j];
				}
				element_ids_(i) = id_[i];
			}
		}
		locations_by_barycenter_ = true;
	}
	else
	{
		locations_by_barycenter_ = false;
	}
}

template <typename MatrixType>
void RegressionData<MatrixType>::setIncidenceMatrix(SEXP RincidenceMatrix)
{
	nRegions_ = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt p = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1];

	incidenceMatrix_.resize(nRegions_, p);

	for (auto i = 0; i < nRegions_; ++i)
	{
		for (auto j = 0; j < p; ++j)
		{
			incidenceMatrix_(i, j) = INTEGER(RincidenceMatrix)[i + nRegions_ * j];
		}
	}
}

template <typename MatrixType>
void RegressionData<MatrixType>::printObservations(std::ostream &out) const
{

	for (auto i = 0; i < observations_.size(); i++)
	{
		out << i << "\t" << observations_(i) << std::endl;
	}
}

template <typename MatrixType>
void RegressionData<MatrixType>::printCovariates(std::ostream &out) const
{

	for (auto i = 0; i < covariates_.rows(); i++)
	{
		for (auto j = 0; j < covariates_.cols(); j++)
		{
			out << covariates_(i, j) << "\t";
		}
		out << std::endl;
	}
}

template <typename MatrixType>
void RegressionData<MatrixType>::printLocations(std::ostream &out) const
{
	if (locations_.ncols() == 2)
		for (UInt i = 0; i < locations_.nrows(); i++)
			out << getLocations<2>(i) << std::endl;
	else
		for (UInt i = 0; i < locations_.nrows(); i++)
			out << getLocations<3>(i) << std::endl;
}

template <typename MatrixType>
void RegressionData<MatrixType>::printIncidenceMatrix(std::ostream &out) const
{
	for (auto i = 0; i < incidenceMatrix_.rows(); i++)
	{
		for (auto j = 0; j < incidenceMatrix_.cols(); j++)
		{
			out << incidenceMatrix_(i, j) << "\t";
		}
		out << std::endl;
	}
}

// -- GAM CONSTRUCTORS --
// Laplace
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
														SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
														SEXP Rmax_num_iteration, SEXP Rthreshold) : RegressionData<>(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ = REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	this->isGAM = true;
}

// PDE
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
														SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
														SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold) : RegressionDataElliptic(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc,
																																													Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ = REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	this->isGAM = true;
}

// PDE SpaceVarying
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
														SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
														SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold) : RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru,
																																																Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ = REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	this->isGAM = true;
}

#endif