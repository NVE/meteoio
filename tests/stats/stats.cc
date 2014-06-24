#include <cstdlib>
#include <time.h>
#include <meteoio/MeteoIO.h>

using namespace std;
using namespace mio;

const double rand_range = 1000.;

void cr_fixed_vectors(vector<double> &X, vector<double> &Y) {
	X.clear(); X.resize(10);
	Y.clear(); Y.resize(10);

	X[0] = -499.882; Y[0] = 421.516;
	X[1] = -114.998; Y[1] = 184.937;
	X[2] = -103.033; Y[2] = 241.216;
	X[3] = -86.3743; Y[3] = -140.725;
	X[4] = 57.4737; Y[4] = -377.073;
	X[5] = 127.179; Y[5] = -376.651;
	X[6] = 250.781; Y[6] = 24.0065;
	X[7] = -999.; Y[7] = 496.798;
	X[8] = 393.442; Y[8] = -999.;
	X[9] = 408.396; Y[9] = 105.45;
}

void cr_rand_vectors(vector<double> &X, vector<double> &Y) {
	const size_t N = 10;
	srand( time(NULL) );
	X.clear(); X.resize(N);
	Y.clear(); Y.resize(N);

	for(size_t ii=0; ii<N; ++ii) {
		X[ii] = rand()/(double)RAND_MAX*rand_range - rand_range/2.;
		Y[ii] = rand()/(double)RAND_MAX*rand_range - rand_range/2.;
	}
}

void print_vector(const vector<double>& X) {
	for(size_t ii=0; ii<X.size(); ++ii)
		std::cout << "X[" << ii << "] = " << X[ii] << "\n";
}

void print_vectors(const vector<double>& X, const vector<double>& Y) {
	for(size_t ii=0; ii<X.size(); ++ii)
		std::cout << "X[" << ii << "] = " << X[ii] << "; Y[" << ii << "] = " << Y[ii] << ";\n";
}

////////////////////////////////// start real tests
bool check_sort(const vector<double>& x, const vector<double>& y) {
	vector<double> X(x), Y(y);
	Interpol1D::sort(X, Y);
	bool status = true;

	for(size_t ii=1; ii<X.size(); ++ii) {
		const double value = X[ii];
		if(value < X[ii-1]) status=false;

		const size_t index = find(x.begin(), x.end(), value) - x.begin();
		if(y[index] != Y[ii]) status=false;
	}

	if(status)
		std::cout << "Sorting of vectors: success\n";
	else
		std::cout << "Sorting of vectors: failed\n";
	return status;
}

bool check_bin(const vector<double>& x, const vector<double>& y) {
	vector<double> X(x), Y(y);
	Interpol1D::equalBin(3, X, Y);
	const bool status1 = IOUtils::checkEpsilonEquality(X[0],-272.812,1e-3) &&
	                     IOUtils::checkEpsilonEquality(X[1],181.327,1e-3) &&
	                     IOUtils::checkEpsilonEquality(X[2],635.466,1e-3) &&
	                     IOUtils::checkEpsilonEquality(Y[0],176.736,1e-3) &&
	                     IOUtils::checkEpsilonEquality(Y[1],-432.179,1e-3) &&
	                     IOUtils::checkEpsilonEquality(Y[2],602.248,1e-3) ;

	/*vector<double> X2(x), Y2(y);
	Interpol1D::equalCountBin(2, X2, Y2);
	const bool status2 = IOUtils::checkEpsilonEquality(X2[0],-201.072,1e-3) &&
	                     IOUtils::checkEpsilonEquality(X2[1],210.957,1e-3) &&
	                     IOUtils::checkEpsilonEquality(Y2[0],176.736,1e-3) &&
	                     IOUtils::checkEpsilonEquality(Y2[1],-156.067,1e-3) ;*/
	const bool status2 = true;

	const bool status = status1 && status2;
	if(status)
		std::cout << "Binning: success\n";
	else
		std::cout << "Binning: failed\n";
	return status;
}

bool check_quantiles(const vector<double>& X) {
	static const double arr[] = {.1, .2, .4, .5, .75, .95};
	static const double results[] = {-191.9748, -107.819, -57.6047, 57.4737, 250.781, 402.4144};
	const vector<double> quartiles(arr, arr + sizeof(arr) / sizeof(arr[0]));
	const vector<double> quantiles = Interpol1D::quantiles(X, quartiles);

	bool status = true;
	for(size_t ii=0; ii<quantiles.size(); ++ii) {
		if(!IOUtils::checkEpsilonEquality(quantiles[ii], results[ii], 1e-6)) {
			std::cerr << setprecision(18) << "Quantile[" << ii << "] should be " << results[ii] << ", computed " << quantiles[ii] << " instead\n";
			status=false;
		}
	}

	if(status)
		std::cout << "Quantiles: success\n";
	else
		std::cout << "Quantiles: failed\n";
	return status;
}

bool check_basics(const vector<double>& X, const vector<double>& Y) {
	//min, max
	const double min = Interpol1D::min_element(Y);
	const double min_ref = -377.073;
	const double max = Interpol1D::max_element(Y);
	const double max_ref = 496.798;
	const bool min_status = IOUtils::checkEpsilonEquality(min, min_ref, 1e-4);
	const bool max_status = IOUtils::checkEpsilonEquality(max, max_ref, 1e-4);
	if(min_status && max_status)
		std::cout << "Min/max: success\n";
	else {
		std::cout << setprecision(12) << "Min should be " << min_ref << ", computed " << min << "\n";
		std::cout << setprecision(12) << "Max should be " << max_ref << ", computed " << max << "\n";
	}

	//median
	const double median = Interpol1D::getMedian(X);
	const double median_ref = 57.4737;
	const bool median_status = IOUtils::checkEpsilonEquality(median, median_ref, 1e-6);
	if(median_status)
		std::cout << "Median: success\n";
	else
		std::cout << setprecision(12) << "Median should be " << median_ref << ", computed " << median << " instead\n";

	//MAD
	const double mad = Interpol1D::getMedianAverageDeviation(X);
	const double mad_ref = 172.4717;
	const bool mad_status = IOUtils::checkEpsilonEquality(mad, mad_ref, 1e-6);
	if(mad_status)
		std::cout << "MAD: success\n";
	else
		std::cout << setprecision(12) << "MAD should be " << mad_ref << ", computed " << mad << " instead\n";

	//variance
	const double variance = Interpol1D::variance(X);
	const double variance_ref = 83038.1246275;
	const bool variance_status = IOUtils::checkEpsilonEquality(variance, variance_ref, 1e-6);
	if(variance_status)
		std::cout << "Variance: success\n";
	else
		std::cout << setprecision(12) << "Variance should be " << variance_ref << ", computed " << variance << " instead\n";

	//stddev
	const double stddev = Interpol1D::std_dev(X);
	const double stddev_ref = 288.163364478;
	const bool stddev_status = IOUtils::checkEpsilonEquality(stddev, stddev_ref, 1e-6);
	if(stddev_status)
		std::cout << "Stddev: success\n";
	else
		std::cout << setprecision(12) << "Stddev should be " << stddev_ref << ", computed " << stddev << " instead\n";

	//weighted mean
	const double d1 = 288.1643545;
	const double d2 = 384.1562055;
	const double w1 = 0.232326, w2 = 0.68125;
	const double m1_ref = 310.465757275, m2_ref = 353.558802994;
	const double m1 = Interpol1D::weightedMean(d1, d2, w1);
	const double m2 = Interpol1D::weightedMean(d1, d2, w2);
	double wmean_status = true;
	if(IOUtils::checkEpsilonEquality(m1, m1_ref, 1e-6) && IOUtils::checkEpsilonEquality(m2, m2_ref, 1e-6))
		std::cout << "Weighted mean: success\n";
	else {
		std::cout << setprecision(12) << "Weighted means should be " << m1_ref << " and " << m2_ref;
		std::cout <<  ", computed " << m1 << " and " << m2 << " instead\n";
		wmean_status = false;
	}

	//vector weighted mean
	const vector<double> weights(X.size(), 1./(X.size()));
	const double vector_mean = Interpol1D::weightedMean(X, weights);
	const double mean = Interpol1D::arithmeticMean(X);
	const double vector_mean_status = IOUtils::checkEpsilonEquality(vector_mean, mean, 1e-6);
	if(vector_mean_status)
		std::cout << "Vector weighted mean: success\n";
	else
		std::cout << setprecision(12) << "Vector mean should be " << mean << ", conputed " << vector_mean << " instead\n";

	return (min_status || max_status || mad_status || variance_status || median_status || stddev_status || wmean_status || vector_mean_status);
}

bool check_covariance(const vector<double>& x, const vector<double>& y) {
	const double cov = Interpol1D::covariance(x, y);
	const double cov_ref = -35272.1266148;

	const bool status = IOUtils::checkEpsilonEquality(cov, cov_ref, 1e-6);
	if(status)
		std::cout << "Covariance: success\n";
	else
		std::cout << setprecision(12) << "Covariance should be " << cov_ref << ", conputed " << cov << " instead\n";

	return status;
}

bool check_derivative(const vector<double>& x, const vector<double>& y) {
	static const double results[] = {-999., -0.6146761, -0.454329, -11.377355, -3.8521071, -1.104764, 2.0748285, 3.241513, 0.516724, -999.};
	vector<double> X(x), Y(y);
	Interpol1D::sort(X, Y);
	const vector<double> der = Interpol1D::derivative(X, Y);

	bool status = true;
	for(size_t ii=0; ii<der.size(); ++ii) {
		if(!IOUtils::checkEpsilonEquality(der[ii], results[ii], 1e-6)) {
			std::cerr << setprecision(12) << "Derivative[" << ii << "] should be " << results[ii] << ", conputed " << der[ii] << " instead\n";
			status = false;
		}
	}

	if(status)
		std::cout << "Derivative: success\n";
	else
		std::cout << "Derivative: failed\n";
	return status;
}

int main() {
	vector<double> x,y;
	//cr_rand_vectors(x, y);
	cr_fixed_vectors(x, y);

	const bool sort_status = check_sort(x,y);
	const bool bin_status = check_bin(x,y);
	const bool basics_status = check_basics(x,y);
	const bool covariance_status = check_covariance(x,y);
	const bool der_status = check_derivative(x,y);
	const bool quantiles_status = check_quantiles(x);

	if(!basics_status || !sort_status || !bin_status || !quantiles_status || !covariance_status || !der_status)
		throw IOException("Statistical functions error!", AT);


	return 0;
}
