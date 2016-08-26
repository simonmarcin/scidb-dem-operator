#include "Firdem2.h"

/*
 * Addapted version of algorithm by:
 * 
 * Plowman, Joseph, Charles Kankelborg, and Petrus Martens. 
 * "Fast differential emission measure inversion of solar coronal data." 
 * The Astrophysical Journal 771.1 (2013): 2.
 * 
 */

using namespace std;
using namespace boost::numeric::ublas;

Firdem2::Firdem2()
{

	// Init vectors and matrices
	tr_norm.resize(nchan);
	datavec.resize(nchan);
	datavec_diff.resize(nchan);
	data2vec.resize(nchan);
	b_inv.resize(nchan, nchan);
	b_inv_tmp.resize(nchan, nchan);
	sigs2_diag.resize(nchan, nchan);

	// Values from iterate
	dem_test.resize(40);
	dempos_test_iter.resize(40);
	datapos_extrap.resize(6);
	chi2totarr.resize(6);
	std::fill(chi2totarr.data().begin(), chi2totarr.data().end(), (1.0 / 6));
	deltatot_iter.resize(1, 6);
	dem_noextrap.resize(40);
	lastdem.resize(40);
	coffs_iter.resize(6);
	deltadem.resize(40);
	tmp_mat_nchan_1.resize(1, 6);

	data_out.resize(6);
	errsin.resize(6);
	datapos.resize(6);
	dem_initial.resize(40);
	dem.resize(40);
	dem_out.resize(40);
	coffs.resize(40);
	dempos_test.resize(40);
	datapos_test.resize(6);
	deltatot.resize(6);
	mincond.resize(nconfs);
	dempos.resize(40);
	dempos_out.resize(40);
	deltadata.resize(6);
	wp.resize(6, 6);
	ainv_arr.resize(6, 6);

	// Values for temp results. needed to speed up prod
	tmp_vec_nchan.resize(nchan);
	tmp_vec_nchan2.resize(nchan);
	tmp_mat_nchan.resize(6, 6);

}


Firdem2::~Firdem2()
{
}


/* Matrix inversion routine.
Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T>(A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}



/* Regulazie input data. */
void Firdem2::firdem_regularize_data(boost::numeric::ublas::vector<int>& datain)
{
	// init variables (before quick tests)
	total = 0;

	// Check to see if the data values are so small that zero is an acceptable solution
	for (int n = 0; n < aStruc.nchan; n++)
	{
		total += pow(pow((datain(n) / errsin(n)), 2), 2);
	}

	if (total < chi2_reg_current) {
		for (int i = 0; i < aStruc.nchan; i++) {
			data_out(i) = 0.0;
		}
		return;
	}

	// Essentially no regularization required
	if (chi2_reg_current / aStruc.nchan < 0.0001) {
		for (int i = 0; i < aStruc.nchan; i++) {
			data_out(i) = datain(i);
		}
		return;
	}

	// init variables (after quick tests)
	bisect_start = false;
	chi2 = 0.0;
	chi2_high = 0.0;
	alpha_low = 0.0;
	alpha_high = 0.0;
	alpha = 0.0;
	tmp = 0.0;

	std::fill(datavec_diff.data().begin(), datavec_diff.data().end(), 0);
	std::fill(data2vec.data().begin(), data2vec.data().end(), 0);
	std::fill(b_inv.data().begin(), b_inv.data().end(), 0);
	std::fill(b_inv_tmp.data().begin(), b_inv_tmp.data().end(), 0);
	std::fill(sigs2_diag.data().begin(), sigs2_diag.data().end(), 0);

	// a_array is computed with normalized temperature response functions so that it has unit
	// diagonals and off - diagonal elements less than one.The data vector must be divided by
	// the same normalization factor to match.
	// a_inv_scaled = a_inv/(tr_norm#tr_norm) is in preprocessing
	noalias(tr_norm) = element_prod(aStruc.tr_norm, aStruc.exptimes);
	for (int n = 0; n < aStruc.nchan; n++) {
		sigs2_diag(n, n) = pow((1.0 / errsin(n)), 2);
	}

	// put dynamic part of datavec to datavec_diff inside the loop
	noalias(datavec) = element_div(datain, tr_norm);

	// loop of firdem_regularize_data
	for (int k = 0; k < niter_max_reg; k++) {

		// Intermediate vector appearing multiple times in the calculation
		axpy_prod(aStruc.a_inv, datavec, tmp_vec_nchan, true);
		noalias(data2vec) = element_div(tmp_vec_nchan, tr_norm);

		// Multiplicative factor ensuring chi squared between the regularized data vectors
		// for this step vs.the previous step is chi20
		if (k == 0) {
			for (int i = 0; i < aStruc.nchan; i++) {
				alpha += pow(errsin(i), 2)*pow(data2vec(i), 2);
			}
			alpha = sqrt(9.0*chi2_reg_current / alpha);
		}

		// The difference between the regularized data vector for this step and the previous setp
		for (int i = 0; i < aStruc.nchan; i++) {
			for (int j = 0; j < aStruc.nchan; j++) {
				b_inv_tmp(i, j) = alpha*aStruc.a_inv_scaled(i, j) + sigs2_diag(i, j);
			}
		}
		InvertMatrix(b_inv_tmp, b_inv);

		// datavec(*) = datavec-(alpha*b_inv#data2vec)/tr_norm
		for (int i = 0; i < aStruc.nchan; i++) {
			tmp = 0.0;
			for (int j = 0; j < aStruc.nchan; j++) {
				tmp += b_inv(i, j)*data2vec(j);
			}
			datavec_diff(i) = datavec(i) - (alpha*tmp) / tr_norm(i);
		}

		// Finish if chi squared between datavec and datavec0 reaches chi2_end
		chi2 = 0.0;
		for (int i = 0; i < aStruc.nchan; i++) {
			chi2 += pow(((datain(i) - datavec_diff(i)*tr_norm(i)) / errsin(i)), 2);
		}
		if (k == 0) { chi2_high = chi2; }

		// Update aplha
		if (bisect_start) {
			if (chi2 < chi2_reg_current) { alpha_low = alpha; }
			if (chi2 > chi2_reg_current) { alpha_high = alpha; }
			alpha = alpha_low + 0.5*(alpha_high - alpha_low);
		}
		if (chi2 < chi2_reg_current && !bisect_start) { alpha = 5 * alpha; }
		if (chi2 > chi2_reg_current && !bisect_start) {
			bisect_start = true;
			chi2_high = chi2;
			alpha_high = alpha;
			alpha = alpha_low + 0.5*(alpha_high - alpha_low);
		}

		// End loop if:
		if (abs(chi2 - chi2_reg_current) / chi2_reg_current < chi2_tol) { break; }

	}

	// Calculate the result      
	noalias(data_out) = element_prod(datavec_diff, tr_norm);


}




/* Core of iterative process. */
void Firdem2::firdem_iterate(int& iter_mincond) {

	// Set up for iteration
	dem = dem_initial;
	std::fill(dempos.begin(), dempos.end(), 0);

	for (int i = 0; i < 40; i++) {
		if (dem_initial(i) > 0) { dempos(i) = dem_initial(i); }
	}
	dempos_out = dempos;
	//ToDo: if(n_elements(datain) eq 0) then datain = reform(a2_array#dem)

	// Data values for DEM with negatives zeroed
	axpy_prod(aStruc.a2_array, dempos, datapos, true);
	noalias(deltadata) = datapos - data_out;

	// Not needed:
	//if(n_elements(mincond) eq 0) then mincond = 0.02 

	std::fill(wp.data().begin(), wp.data().end(), 0);
	for (int i = 0; i < 6; i++) {
		if (aStruc.wvec(i) >= mincond(iter_mincond)*aStruc.max_wvec) {
			wp(i, i) = aStruc.wpvec(i);
		}
	}
	//ainv_arr = prod((aStruc.v), ( (boost::numeric::ublas::matrix<double>)prod(wp, aStruc.u)));
	axpy_prod(wp, aStruc.u2, tmp_mat_nchan, true);
	axpy_prod(aStruc.v2, tmp_mat_nchan, ainv_arr, true);

	// 	coffs = prod((ainv_arr), (element_div(datain, aStruc.normfac)));
	noalias(tmp_vec_nchan) = element_div(data_out, aStruc.normfac);
	axpy_prod(ainv_arr, tmp_vec_nchan, coffs_iter, true);
	int nbad_recent = 0;

	// When extrapolating, DEMs are initially assigned to these.If they improve chi squared,
	// they are copied into the main DEM and data arrays
	std::fill(dem_test.data().begin(), dem_test.data().end(), 0);
	std::fill(dempos_test_iter.data().begin(), dempos_test_iter.data().end(), 0);
	std::fill(datapos_extrap.data().begin(), datapos_extrap.data().end(), 0);
	std::fill(chi2totarr.data().begin(), chi2totarr.data().end(), (1.0 / 6));
	std::fill(deltatot_iter.data().begin(), deltatot_iter.data().end(), 0);

	// Compute inital chi squared
	datapos = element_div(datapos, errsin);
	row(deltatot_iter, 0) = (boost::numeric::ublas::vector<double>)element_prod(datapos, datapos);

	//double chi2min = (prod((element_prod(deltatot, deltatot)), chi2totarr))(0);
	noalias(tmp_mat_nchan_1) = element_prod(deltatot_iter, deltatot_iter);
	boost::numeric::ublas::vector<double> tmp_vec_1(1);
	axpy_prod(tmp_mat_nchan_1, chi2totarr, tmp_vec_1, true);
	double chi2min = tmp_vec_1(0);

	double chi2;
	double chi2_test;
	double chi2_noextrap;
	double lastchi2 = 0;
	int itcounter;

	// Begin iteration
	for (itcounter = 0; itcounter < niter_max - 1; itcounter++) {

		// Stop if we've gone over nbad_recent points since last chi squared improvemen
		if (nbad_recent > nbad_iter_max) { break; }

		// Compute DEM correction coefficients in the instrument response basis.
		// These will restore the DEM to agreement with the data, but generally
		// introduce some negative emission.
		//coffs = prod((ainv_arr), (element_div(deltadata, aStruc.normfac)));
		noalias(tmp_vec_nchan) = element_div(deltadata, aStruc.normfac);
		axpy_prod(ainv_arr, tmp_vec_nchan, coffs_iter, true);

		// Convert these coefficients to basis2
		//deltadem = prod(aStruc.basis22, coffs);
		axpy_prod(aStruc.basis22, coffs_iter, deltadem, true);
		dem_initial = dempos - deltadem;
		for (int i = 0; i < 40; i++) {
			if (dem_initial(i) < 0) { dempos(i) = 0; }
			else { dempos(i) = dem_initial(i); }
		}

		//Compute difference between data for positive DEM and initial data
		//datapos = prod(aStruc.a2_array, dempos);
		axpy_prod(aStruc.a2_array, dempos, datapos, true);
		deltadata = datapos - data_out;

		// Compute chi squared for this step
		//row(deltatot, 0) = element_prod(element_div(deltadata, sigmas), element_div(deltadata, sigmas));
		noalias(tmp_vec_nchan) = element_div(deltadata, errsin);
		tmp_vec_nchan = element_prod(tmp_vec_nchan, tmp_vec_nchan);
		row(deltatot_iter, 0) = tmp_vec_nchan;

		//chi2 = (prod(deltatot, chi2totarr))(0);
		axpy_prod(deltatot_iter, chi2totarr, tmp_vec_1, true);
		chi2 = tmp_vec_1(0);

		dem_noextrap = dem_initial;
		chi2_noextrap = chi2;

		if (itcounter >= aStruc.extrap_start + 1) {

			// If current chi squared is better than the previous one, extrapolate. 
			// Otherwise, do nothing
			if (chi2 < lastchi2) {

				//Simple linear extrapolation scheme: Change is computed between
				// the most recent DEM and the previous DEM(before extrapolation).
				// This is multiplied by extrap_fac*itcount and added to the most recent DEM
				noalias(deltadem) = ((dem_initial - lastdem)* aStruc.extrap_fac)*itcounter;
				dem_test = dem_initial + deltadem;
				for (int i = 0; i < 40; i++) {
					if (dem_test(i) < 0) { dempos_test_iter(i) = 0; }
					else { dempos_test_iter(i) = dem_test(i); }
				}

				// Compute chi squared resulting from the extrapolation
				//datapos_extrap = prod(aStruc.a2_array, dempos_test);
				axpy_prod(aStruc.a2_array, dempos_test_iter, datapos_extrap, true);
				row(deltatot_iter, 0) = element_prod(element_div((datapos_extrap - data_out), errsin),
					element_div((datapos_extrap - data_out), errsin));
				//chi2_test = (prod(deltatot, chi2totarr))(0);
				axpy_prod(deltatot_iter, chi2totarr, tmp_vec_1, true);
				chi2_test = tmp_vec_1(0);

				// If the resulting chi squared is better than the previous one,
				// update the DEM with the extrapolation
				if (chi2_test < chi2) {
					chi2 = chi2_test;
					datapos = datapos_extrap;
					deltadata = datapos - data_out;
					dem_initial = dem_test;
					dempos = dempos_test_iter;
				}
			}
		}

		// If chi squared got better, update output chi squared and dems
		if (chi2 < chi2min) {
			nbad_recent = 0;
			chi2min = chi2;
			dem = dem_initial;
			dempos_out = dempos;
			// Finish if chi squared is good enough
			if (chi2 < chi2_iter_thold) { break; }

		}

		lastdem = dem_noextrap;
		lastchi2 = chi2_noextrap;
		nbad_recent++;

	}

	// Finished iterating. dem contains the DEM corresponding to the best chi squared
	// (with negatives zeroed) found.Although chi squared is calculated with negatives
	// zeroed, the returned DEM retains them for good measure(the user will generally
	// want to zero the negatives before using the DEMs; a trivial exercise)
	//datapos = prod(aStruc.a2_array, dempos_out);
	axpy_prod(aStruc.a2_array, dempos_out, datapos, true);
	row(deltatot_iter, 0) = element_prod(element_div((datapos - data_out), errsin),
		element_div((datapos - data_out), errsin));
	iters = itcounter;
	//chi2_final = deltatot#chi2totarr
}



void Firdem2::calculateDem(boost::numeric::ublas::vector<int>& datain,
	boost::numeric::ublas::vector<double>& coffs, int& itcount) {

	//Calculate standard AIA_Errors
	for (int i = 0; i < 6; i++) {
		errsin(i) = (int)sqrt(std::max(datain(i)*aStruc.dnpp(i), 0.0) + (aStruc.rdn(i)*aStruc.rdn(i)));
	}
	for (int i = 0; i < 6; i++) {
		if (datain(i) < 0) { datain(i) = 0; }
	}


	// inside each pixel (i,j), do:
	for (int k = 0; k < aStruc.nconfs; k++) {

		// calculate chi2 helber variables
		chi2_reg_current = aStruc.chi2_reg_ends(k);
		chi2_iter_thold = std::max((chi2_target - aStruc.chi2_reg_ends(0)), (chi2_target - aStruc.chi2_reg_ends(k))) / nchan;

		// Regularize the data so that its chi squared(relative to the original data)
		// is equal to chi2_reg_current :
		firdem_regularize_data(datain);

		//Compute the first pass DEM corresponding to the regularized data
		//dem_initial = prod(aStruct.basis22, (boost::numeric::ublas::vector<double>)prod(aStruct.a_inv, element_div(data_out, aStruct.normfac)));
		noalias(tmp_vec_nchan) = element_div(data_out, aStruc.normfac);
		axpy_prod(aStruc.a_inv, tmp_vec_nchan, tmp_vec_nchan2, true);
		axpy_prod(aStruc.basis22, tmp_vec_nchan2, dem_initial, true);

		// Check to see if the first pass DEM needs iteration to remove negative emission
		for (int i = 0; i < 40; i++) {
			if (dem_initial(i) > 0) {
				dempos_test(i) = dem_initial(i);
			}
			else {
				dempos_test(i) = 0;
			}
		}

		//datapos_test = prod(aStruct.a2_array, dempos_test);
		axpy_prod(aStruc.a2_array, dempos_test, datapos_test, true);
		noalias(deltatot) = element_div((datapos_test - datain), errsin);
		deltatot = element_prod(deltatot, deltatot);
		chi2_current = 0.0;
		for (int i = 0; i < nchan; i++) {
			chi2_current += deltatot(i) / nchan;
		}
		

		if (chi2_current <= aStruc.chi2thold) {
			chi2 = chi2_current;
			coffs = dempos_test;
			break;
		}
		else {

			//Attempt to iterate away negative emission in the first pass DEM
			firdem_iterate(k);
			itcount += iters;

			//Compute chi squared relative to the original data
			deltatot = element_prod(element_div((datapos - datain), errsin),
				element_div((datapos - datain), errsin));

			chi2_current = 0;
			for (int i = 0; i < nchan; i++) {
				chi2_current += deltatot(i) / nchan;
			}
			if (k == 0) { chi2 = chi2_current; }

			// Update the DEM if the new chi squared is better than the previous best
			if (k == 0 || chi2_current < chi2) {
				chi2 = chi2_current;
				coffs = dem;
			}
			if (chi2 < aStruc.chi2thold) { break; }

		}
	}

	//chi2s(i, j) = chi2
	//its(i, j) = itcount

}



