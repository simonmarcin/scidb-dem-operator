#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include "PreProcesing.h"

class Firdem2
{
public:
	Firdem2();
	void calculateDem(boost::numeric::ublas::vector<int>& datain,
		boost::numeric::ublas::vector<double>& coffs, int& itcount);
	~Firdem2();

private:


	void firdem_regularize_data(boost::numeric::ublas::vector<int>& datain);
	void firdem_iterate(int& iter_mincond);

	// Values from preprocessing.
	double chi2_target = 12.5916;
	int nchan = 6;
	int nconfs = 3;
	double chi2_current = 0.0;
	double chi2 = 0.0;
	double chi2_reg_current;
	double chi2_iter_thold;
	int niter_max = 1000;
	int nbad_iter_max = 250;
	int iters = 0;
	int itcount = 0;

	// Values used in regulaize function
	int niter_max_reg = 50;
	bool bisect_start = false;
	double chi2_regularize = 0.0;
	double chi2_tol = 0.05;
	double chi2_high = 0.0;
	double alpha_low = 0.0;
	double alpha_high = 0.0;
	double alpha = 0.0;
	double tmp = 0.0;
	double total = 0;
	boost::numeric::ublas::vector<double> tr_norm;
	boost::numeric::ublas::vector<double> datavec;
	boost::numeric::ublas::vector<double> datavec_diff;
	boost::numeric::ublas::vector<double> data2vec;
	boost::numeric::ublas::matrix<double> b_inv;
	boost::numeric::ublas::matrix<double> b_inv_tmp;
	boost::numeric::ublas::matrix<double> sigs2_diag;

	// Values used in iterate function
	boost::numeric::ublas::vector<double> dem_test;
	boost::numeric::ublas::vector<double> dempos_test_iter;
	boost::numeric::ublas::vector<double> datapos_extrap;
	boost::numeric::ublas::vector<double> chi2totarr;
	boost::numeric::ublas::matrix<double> deltatot_iter;
	boost::numeric::ublas::vector<double> dem_noextrap;
	boost::numeric::ublas::vector<double> lastdem;
	boost::numeric::ublas::vector<double> coffs_iter;
	boost::numeric::ublas::vector<double> deltadem;
	boost::numeric::ublas::matrix<double> tmp_mat_nchan_1;

	boost::numeric::ublas::vector<int> errsin;
	boost::numeric::ublas::vector<double> data_out;
	boost::numeric::ublas::vector<double> datapos;
	boost::numeric::ublas::vector<double> dem_initial;
	boost::numeric::ublas::vector<double> dem;
	boost::numeric::ublas::vector<double> dem_out;
	boost::numeric::ublas::vector<double> coffs;
	boost::numeric::ublas::vector<double> dempos_test;
	boost::numeric::ublas::vector<double> datapos_test;
	boost::numeric::ublas::vector<double> deltatot;
	boost::numeric::ublas::vector<double> mincond;
	boost::numeric::ublas::vector<double> dempos;
	boost::numeric::ublas::vector<double> dempos_out;
	boost::numeric::ublas::vector<double> deltadata;

	boost::numeric::ublas::matrix<double> wp;
	boost::numeric::ublas::matrix<double> ainv_arr;

	// Values for temp results. needed to speed up prod
	boost::numeric::ublas::vector<double> tmp_vec_nchan;
	boost::numeric::ublas::vector<double> tmp_vec_nchan2;
	boost::numeric::ublas::matrix<double> tmp_mat_nchan;

	int tot_iterations = 0;

	// Init a preprocessing Object
	PreProcesing aStruc;

};

