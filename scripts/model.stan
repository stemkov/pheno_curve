functions{
	
	// linear interpolation from https://discourse.mc-stan.org/t/linear-interpolation-and-searchsorted-in-stan/13318/5
	// modified slightly to declare all vars at top of block
	real linear_interpolation_v(real x, vector x_pred, vector y_pred){
    	int K = rows(x_pred);
    	vector[K] deltas = x - x_pred;
		real ans;
		real t;
		real w;
		int i;
		real x1;
		real x2;
		real y1;
		real y2;
		if(x<x_pred[1] || x>x_pred[K]) reject("x is outside of the x_pred grid!");
		if(rows(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
		//this is which.min()
		i = sort_indices_asc(fabs(deltas))[1];
		if(deltas[i]<=0) i -= 1;
		ans = y_pred[i];
		x1 = x_pred[i];
		x2 = x_pred[i + 1];
		y1 = y_pred[i];
		y2 = y_pred[i + 1];
		ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
		t = (x-x1)/(x2-x1);
		w = 1-t;
		//w = 1/(1 + exp(1/(1-t) - 1/t));
		//w = 1 - 3*pow(t,2) + 2*pow(t,3);
		ans = w*y1 + (1-w)*y2;
		return ans;
  }
	
  real verhulst(int t_given, real B0, real e, real f, real m, real onset){
		
		int t_length = 150; // maximum length of simulation - make sure it's greater than biggest onset+t_given
		real t_out = 1.0*t_given - onset; // 1.0 is int-to-real trick
		real P_out;
		vector [t_length] P_vec;
		vector [t_length] t_vec;
		for (i in 1:(t_length)) t_vec[i] = i-1; // -1 because time at index 1 should actually be 0
		
		
		if (onset >= t_given){
			return 0.0;
		}
		else{
		
			// shift stuff
			// int t_length = t_max - shift; // not allowed bc real and int don't mix
			
			
			// initial conditions
		  real P[t_length];
		  real B[t_length];
		  P[1] = 0.0;
		  B[1] = B0;
		  
		  // simulation
		  for(t in 2:(t_length)){
		    B[t] = B[t-1] - e*B[t-1] + f*B[t-1]^2;
		    P[t] = P[t-1] + e*B[t-1] - f*B[t-1]^2 - m*P[t-1];
		  }
			
			P_vec = to_vector(P);
			//t_vec = to_vector(1:(t_length));
			
			// interpolated value between timesteps
			P_out = linear_interpolation_v(t_out, t_vec, P_vec);
			
			
		  // return P[t_max];
		  return P_out;
    }
  }

}

data {
	int n;
	int n_ts;
	int n_year;
	int ts[n];
	int year[n];
	int time[n];
	
	real P[n];
	real flowers[n];
}

parameters {
	vector <lower=100, upper=200>[n_year] onset;
	vector <lower=1>[n_ts] B0;
	// pars for reparametrization (Matt's trick)
	vector [n_year] e_tilde;
	vector [n_year] f_tilde;
	vector [n_year] m_tilde;
	// hyperparameters
	// vector<lower=1>[n_ts] B0_mu; // do I need [n_ts] here?
	// vector<lower=0>[n_ts] B0_sigma;
	real <lower=1> B0_mu;
	real <lower=0> B0_sigma;
	real <lower=100, upper=200> onset_mu;
	real <lower=0, upper=100> onset_sigma;

	real e_tilde_mu;
	real e_tilde_sigma;
	real f_tilde_mu;
	real f_tilde_sigma;
	real m_tilde_mu;
	real m_tilde_sigma;

	// covariates
	real a_flr;

	real <lower=0> sigma;
}

transformed parameters {
	// reparameterization
	vector <lower=0>[n_year]  e;
	vector <lower=0>[n_year]  f;
	vector <lower=0>[n_year]  m;

	e = 0.2 + e_tilde*0.1;
	f = 0.004 + f_tilde*0.001;
	m = 0.2 + m_tilde*0.1;
}

model{
  // priors
	// B0 ~ normal(40,10); // for no pooling
	B0_mu ~ normal(40,10); // hyper mean
	B0_sigma ~ lognormal(0,1); // hyper deviation
	B0 ~ normal(B0_mu, B0_sigma); // for partial pooling

	// onset ~ normal(160, 40);
	onset_mu ~ normal(160, 10);
	onset_sigma ~ lognormal(1,1);
	onset ~ normal(onset_mu, onset_sigma);

	// hyperparameters for rates
	e_tilde_mu ~ normal(0,1);
	e_tilde_sigma ~ lognormal(0,1);
	f_tilde_mu ~ normal(0,1);
	f_tilde_sigma ~ lognormal(0,1);
	m_tilde_mu ~ normal(0,1);
	m_tilde_sigma ~ lognormal(0,1);

	e_tilde ~ normal(e_tilde_mu,e_tilde_sigma); // implies mean 0.2, sigma 0.1
	f_tilde ~ normal(f_tilde_mu,f_tilde_sigma); // implies mean 0.004, sigma 0.001
	m_tilde ~ normal(m_tilde_mu,m_tilde_sigma); // implies mean 0.1, sigma 0.1

	a_flr ~ normal(0,5);

	sigma ~ lognormal(0,1);
  

	for (i in 1:n){
    P[i] ~ normal(verhulst(time[i], B0[ts[i]], e[year[i]], f[year[i]], m[year[i]], onset[year[i]]) + a_flr*flowers[i], sigma);
  }

  
}

