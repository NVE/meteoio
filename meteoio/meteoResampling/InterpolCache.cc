        // Copy constructor
        InterpolARIMA(const InterpolARIMA &other)
            : gap_loc(other.gap_loc), N_gap(other.N_gap), time(other.time), pred_forward(other.pred_forward),
              pred_backward(other.pred_backward), data(other.data), xreg_vec_f(other.xreg_vec_f), xreg_vec_b(other.xreg_vec_b), 
                data_forward(other.data_forward), data_backward(other.data_backward), new_xreg_vec_f(other.new_xreg_vec_f), 
                new_xreg_vec_b(other.new_xreg_vec_b), N_data_forward(other.N_data_forward), N_data_backward(other.N_data_backward),
                max_p(other.max_p), max_d(other.max_d), max_q(other.max_q), start_p(other.start_p), start_q(other.start_q),
                max_P(other.max_P), max_D(other.max_D), max_Q(other.max_Q), start_P(other.start_P), start_Q(other.start_Q),
                r(other.r), s(other.s), method(other.method), opt_method(other.opt_method), stepwise(other.stepwise),
                approximation(other.approximation), num_models(other.num_models), seasonal(other.seasonal), stationary(other.stationary) 
        {    
            std::cout << "copy constructor called" << std::endl;
            auto_arima_forward = auto_arima_copy(other.auto_arima_forward);
            auto_arima_backward = auto_arima_copy(other.auto_arima_backward);
            xreg_f = (xreg_vec_f.empty()) ? NULL : &xreg_vec_f[0]; 
            xreg_b = (xreg_vec_b.empty()) ? NULL : &xreg_vec_b[0];
            new_xreg_f = (new_xreg_vec_f.empty()) ? NULL : &new_xreg_vec_f[0];
            new_xreg_b = (new_xreg_vec_b.empty()) ? NULL : &new_xreg_vec_b[0];
        }

        // Copy assignment operator
        InterpolARIMA& operator=(const InterpolARIMA &other)
        {
            std::cout << "copy assignment operator called" << std::endl;
            if (this != &other) // protect against invalid self-assignment
            {
                auto_arima_forward = auto_arima_copy(other.auto_arima_forward);
                auto_arima_backward = auto_arima_copy(other.auto_arima_backward);

                // 3: copy all the other fields from the other object
                gap_loc = other.gap_loc;
                N_gap = other.N_gap;
                time = other.time;
                pred_forward = other.pred_forward;
                pred_backward = other.pred_backward;
                data = other.data;
                xreg_vec_f = other.xreg_vec_f;
                xreg_vec_b = other.xreg_vec_b;
                data_forward = other.data_forward;
                data_backward = other.data_backward;
                new_xreg_vec_f = other.new_xreg_vec_f;
                new_xreg_vec_b = other.new_xreg_vec_b;
                N_data_forward = other.N_data_forward;
                N_data_backward = other.N_data_backward;
                max_p = other.max_p;
                max_d = other.max_d;
                max_q = other.max_q;
                start_p = other.start_p;
                start_q = other.start_q;
                max_P = other.max_P;
                max_D = other.max_D;
                max_Q = other.max_Q;
                start_P = other.start_P;
                start_Q = other.start_Q;
                r = other.r;
                s = other.s;
                method = other.method;
                opt_method = other.opt_method;
                stepwise = other.stepwise;
                approximation = other.approximation;
                num_models = other.num_models;
                seasonal = other.seasonal;
                stationary = other.stationary;

                // 4: handle the pointers to the vectors
                xreg_f = (xreg_vec_f.empty()) ? NULL : &xreg_vec_f[0]; 
                xreg_b = (xreg_vec_b.empty()) ? NULL : &xreg_vec_b[0];
                new_xreg_f = (new_xreg_vec_f.empty()) ? NULL : &new_xreg_vec_f[0];
                new_xreg_b = (new_xreg_vec_b.empty()) ? NULL : &new_xreg_vec_b[0];
            }
            // by convention, always return *this
            return *this;
        }

        // Destructor
        ~InterpolARIMA() {
            std::cout << "destructor called" << std::endl;
            auto_arima_free(auto_arima_forward);
            auto_arima_free(auto_arima_backward);
            delete xreg_f;
            delete xreg_b;
            delete new_xreg_f;
            delete new_xreg_b;
        }


        	// if im at the end of the data i need to have a specific gap, so that it is not just a single point and not a new gap is computed every time a later point is called


		} else if (position == ResamplingAlgorithms::end && resampling_date > gap.endDate && resampling_date >= gap.startDate && gap.startDate == vecM[vecM.size()-1].date) {
			std::cout << "at date: " << resampling_date.toString(Date::ISO) << " there is a gap at the end of the data" << std::endl;
			// we have already computed a gap after the end of data, now we can use it
			int n_steps = static_cast<int>((resampling_date - gap.startDate).getJulian(true) * gap.sampling_rate);
			std::vector<double> predictions = cache_end_arima.predict(n_steps);
			std::vector<Date> pred_dates(n_steps);
			for (int i = 0; i < n_steps; i++) {
				pred_dates[i] = gap.startDate + i / gap.sampling_rate;
			}
			auto it = findDate(pred_dates, resampling_date);
			if (it != pred_dates.end()) {
				size_t idx = std::distance(pred_dates.begin(), it);
				md(paramindex) = predictions[idx];
				return;
			} else {
				// linearly interpolate the data
				size_t idx = std::distance(pred_dates.begin(), std::lower_bound(pred_dates.begin(), pred_dates.end(), resampling_date));
				md(paramindex) = interpolVecAt(predictions, pred_dates, idx, resampling_date);
				return;
			}
		}
	}



auto_arima_object auto_arima_copy(auto_arima_object original) {
	// Allocate a new auto_arima_set on the heap
	int pqdmax[] = {original->pmax, original->dmax, original->qmax};
	int PQDmax[] = {original->Pmax, original->Dmax, original->Qmax};
	auto_arima_object copy = auto_arima_init(pqdmax, PQDmax, original->s, original->r, original->N);

	// Copy the values from the original to the copy
	copy->N = original->N;
	copy->Nused = original->Nused;
	copy->method = original->method;
	copy->optmethod = original->optmethod;
	copy->pmax = original->pmax;
	copy->dmax = original->dmax;
	copy->qmax = original->qmax;
	copy->Pmax = original->Pmax;
	copy->Dmax = original->Dmax;
	copy->Qmax = original->Qmax;
	copy->p = original->p;
	copy->d = original->d;
	copy->q = original->q;
	copy->s = original->s;
	copy->P = original->P;
	copy->D = original->D;
	copy->Q = original->Q;
	copy->r = original->r;
	copy->M = original->M;
	copy->ncoeff = original->ncoeff;
	copy->lvcov = original->lvcov;
	copy->mean = original->mean;
	copy->var = original->var;
	copy->loglik = original->loglik;
	copy->ic = original->ic;
	copy->retval = original->retval;
	copy->start = original->start;
	copy->imean = original->imean;
	copy->idrift = original->idrift;
	copy->stationary = original->stationary;
	copy->seasonal = original->seasonal;
	copy->Order_max = original->Order_max;
	copy->p_start = original->p_start;
	copy->q_start = original->q_start;
	copy->P_start = original->P_start;
	copy->Q_start = original->Q_start;
	copy->stepwise = original->stepwise;
	copy->num_models = original->num_models;
	copy->approximation = original->approximation;
	copy->verbose = original->verbose;
	copy->alpha_test = original->alpha_test;
	copy->alpha_seas = original->alpha_seas;
	copy->lambda = original->lambda;
	copy->sigma2 = original->sigma2;
	copy->aic = original->aic;
	copy->bic = original->bic;
	copy->aicc = original->aicc;

	memcpy(copy->information_criteria, original->information_criteria, 10 * sizeof(int));
	memcpy(copy->test, original->test, 10 * sizeof(int));
	memcpy(copy->type, original->type, 10 * sizeof(int));
	memcpy(copy->seas, original->seas, 10 * sizeof(int));
	int res_size;
	if (copy->P == 0 && copy->D == 0 && copy->Q==0) {
		res_size = copy->N-copy->d;
	} else {
		res_size = copy->N-copy->D*copy->s-copy->d;
	}
	memcpy(copy->res, original->res, res_size * sizeof(double));

	// If there are any members that are pointers (like phi, theta, etc.), you'll need to allocate new memory for them and copy the values
	copy->phi = (double*)malloc(original->p * sizeof(double));
	memcpy(copy->phi, original->phi, original->p * sizeof(double));
	copy->theta = (double*)malloc(original->q * sizeof(double));
	memcpy(copy->theta, original->theta, original->q * sizeof(double));
	copy->PHI = (double*)malloc(original->P * sizeof(double));
	memcpy(copy->PHI, original->PHI, original->P * sizeof(double));
	copy->THETA = (double*)malloc(original->Q * sizeof(double));
	memcpy(copy->THETA, original->THETA, original->Q * sizeof(double));
	copy->exog = (double*)malloc(original->M * sizeof(double));
	memcpy(copy->exog, original->exog, original->N*original->r * sizeof(double));
	copy->vcov = (double*)malloc(copy->lvcov * sizeof(double));

	int param_size;
	if (copy->r == 0) {
		if (copy->P == 0 && copy->D == 0 && copy->Q==0) {
			param_size = copy->p + copy->q + res_size + copy->lvcov;
		} else {
			param_size = copy->p + copy->q + copy->P + copy->Q + res_size + copy->lvcov;
		}
	} else {
		param_size = copy->p + copy->q + copy->P + copy->Q + copy->M + res_size + copy->lvcov; 
	}
	memcpy(copy->params, original->params, param_size * sizeof(double));

	return copy;
}

auto_arima_object auto_arima_copy(auto_arima_object original);



