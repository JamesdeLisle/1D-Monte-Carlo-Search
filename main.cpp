#include <eigen3/Eigen/Dense>
#include <complex>
#include <random>
#include <ctime>
#include <iostream>
#define PI 3.14159

//function that generates random doubles from within a range
double getRandom( double min, double max ) 
{
	std::uniform_real_distribution<double> dist(min, max);
	std::mt19937 rng;
	rng.seed(std::random_device{}());
	return dist(rng);
}

//function that computes the symmetry equation
//C(H(-p)^*C^{\dagger}+H(p)=\alpha(p) at different
//points in the parameter space
std::vector<double> alphaGenerate( double mu_a, double mu_b, double t_1, double t_2, double p ) //function that computes the 
{
	double t_0, del_t, mu_0, del_mu; 
	mu_0 = ( mu_a + mu_b )/2.0;
	del_mu = ( mu_a - mu_b )/2.0;
	t_0 = ( t_1 + t_2 )/2.0;
	del_t = ( t_1 - t_2 )/2.0;
	Eigen::Matrix2cd Ham;
	Eigen::Matrix2cd HamMp;
	Eigen::Matrix2cd sigZ;
	Eigen::Matrix2cd alpha;
	std::complex<double> I(0.0,1.0);
	std::vector<double> alpha_elements(4);

	sigZ(0,0) = 1.0;
	sigZ(0,1) = 0.0;
	sigZ(1,0) = 0.0;
	sigZ(1,1) = -1.0;

	Ham(0,0) = mu_0 + del_mu;
	Ham(0,1) = ( t_0 + del_t ) + ( t_0 - del_t ) * ( cos( p ) + I * cos( p ) );
	Ham(1,0) = ( t_0 + del_t ) + ( t_0 - del_t ) * ( cos( p ) - I * cos( p ) );
	Ham(1,1) = ( mu_0 - del_mu );

	HamMp(0,0) = mu_0 + del_mu;
	HamMp(0,1) = ( t_0 + del_t ) + ( t_0 - del_t ) * ( cos( -p ) + I * cos( -p ) );
	HamMp(1,0) = ( t_0 + del_t ) + ( t_0 - del_t ) * ( cos( -p ) - I * cos( -p ) );
	HamMp(1,1) = ( mu_0 - del_mu );
	
	alpha = sigZ * Ham.conjugate() * sigZ + HamMp;
	

	alpha_elements[0] = std::abs( alpha(0,0) );
	alpha_elements[1] = std::abs( alpha(0,1) );
	alpha_elements[2] = std::abs( alpha(1,0) );
	alpha_elements[3] = std::abs( alpha(1,1) );
	//std::cout << alpha_elements[3] << std::endl;
	return alpha_elements;
}

int main()
{
	
	std::vector<double> parameters(4), parametersTemp(4);
	double p;
	double mu_a_min = -4.0, mu_a_max = 4.0;
	double mu_b_min = -4.0, mu_b_max = 4.0;
	double t_1_min = -4.0, t_1_max = 4.0;
	double t_2_min = -4.0, t_2_max = 4.0;
	double ball_radius = 1.0e-3;

	//generate initial point in param space	
	parameters[0] = getRandom( mu_a_min, mu_a_max );
	parameters[1] = getRandom( mu_b_min, mu_b_max );
	parameters[2] = getRandom( t_1_min, t_1_max );
	parameters[3] = getRandom( t_2_min, t_2_max );
	p = getRandom( -PI, PI ); 
	
	int stopFlag = 1, holdFlag00 = 1, holdFlag01 = 1, holdFlag11 = 1;
	//compute inital alpha
	std::vector<double> alpha = alphaGenerate( parameters[0], parameters[1], parameters[2], parameters[3], p ), alphaTemp;
	
	double choice_prob;

	while ( stopFlag == 1 )
	{
		if ( holdFlag00 == 1 )
		{
			parametersTemp[0] = getRandom( parameters[0] - ball_radius, parameters[0] + ball_radius );
		}
		if ( holdFlag11 == 1 )
		{
			parametersTemp[1] = getRandom( parameters[1] - ball_radius, parameters[1] + ball_radius );
		}
		if ( holdFlag01 == 1 )
		{
			parametersTemp[2] = getRandom( parameters[2] - ball_radius, parameters[2] + ball_radius );
			parametersTemp[3] = getRandom( parameters[3] - ball_radius, parameters[3] + ball_radius );
		}
		
		alphaTemp = alphaGenerate( parametersTemp[0], parametersTemp[1], parametersTemp[2], parametersTemp[3], p );
		if ( alphaTemp[0] == 0.0 ) { holdFlag00 == 0; }
		if ( alphaTemp[1] == 0.0 ) { holdFlag01 == 0; }
		if ( alphaTemp[3] == 0.0 ) { holdFlag11 == 0; }
		
		if ( holdFlag00 == 1 )
		{
			if ( alphaTemp[0] <= alpha[0] )
			{
				parameters[0] = parametersTemp[0];
				alpha[0] = alphaTemp[0];
			}
			else
			{
				choice_prob = getRandom( 1.0, 2.0 );
				if ( choice_prob < 1.00001 )
				{
					parameters[0] = parametersTemp[0];
					alpha[0] = alphaTemp[0];
				}
			}
		}

		if ( holdFlag01 == 1 )
		{
			if ( alphaTemp[1] <= alpha[1] )
			{
				parameters[2] = parametersTemp[2];
				parameters[3] = parametersTemp[3];
				alpha[1] = alphaTemp[1];
				alpha[2] = alphaTemp[2];
			}
			else
			{
				choice_prob = getRandom( 1.0, 2.0 );
				if ( choice_prob < 1.00001 )
				{
					parameters[2] = parametersTemp[2];
					parameters[3] = parametersTemp[3];
					alpha[1] = alphaTemp[1];
					alpha[2] = alphaTemp[2];

				}
			}
		}
		
		if ( holdFlag11 == 1 )
		{
			if ( alphaTemp[3] <= alpha[3] )
			{
				parameters[1] = parametersTemp[1];
				alpha[3] = alphaTemp[3];
			}
			else
			{
				choice_prob = getRandom( 1.0, 2.0 );
				if ( choice_prob < 1.00001 )
				{
					parameters[1] = parametersTemp[1];
					alpha[3] = alphaTemp[3];
				}
			}
		}

		//check if alpha is close to zero
		if (( alpha[0] < 1.0e-4 ) && ( alpha[1] < 1.0e-4 ) && ( alpha[2] < 1.0e-4 ) && ( alpha[3] < 1.0e-4 ))
		{
			stopFlag = 0;
		}
		 
		
	}
	std::cout << "mu_a: " << parameters[0] << " " << "mu_b: " << parameters[1] << " " << "t_1: " << parameters[2] << " " << "t_2: " << parameters[3] << std::endl; 
		
	return 0;	
}




