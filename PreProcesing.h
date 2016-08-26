#pragma once

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

class PreProcesing {

public:

	// Define values
	boost::numeric::ublas::vector<double> chi2_reg_ends;
	boost::numeric::ublas::vector<double> tr_norm;
	boost::numeric::ublas::vector<double> exptimes;
	boost::numeric::ublas::vector<double> normfac;
	boost::numeric::ublas::vector<double> wvec;
	boost::numeric::ublas::vector<double> wpvec;
	boost::numeric::ublas::vector<double> dnpp;
	boost::numeric::ublas::vector<double> rdn;
	boost::numeric::ublas::matrix<double> a_inv;
	boost::numeric::ublas::matrix<double> basis22;
	boost::numeric::ublas::matrix<double> a2_array;
	boost::numeric::ublas::matrix<double> v;
	boost::numeric::ublas::matrix<double> v2;
	boost::numeric::ublas::matrix<double> u;
	boost::numeric::ublas::matrix<double> u2;
	boost::numeric::ublas::matrix<double> a_inv_scaled;
	double chi2thold;
	double extrap_fac;
	double max_wvec;
	int nconfs;
	int nchan;
	int hour;
	int extrap_start;

	PreProcesing ();


};
