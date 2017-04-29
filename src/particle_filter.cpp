/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>

#include "particle_filter.h"
using namespace std;
//#define DEBUG 1

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles=500;
	
	random_device rd;
	default_random_engine gen(rd());
	normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double>dist_y(y,std[1]);
    normal_distribution<double>dist_psi(theta,std[2]);
    
	//double weight_initial=1.0/double(num_particles);
	double weight_initial=1.0;
	for (int i=0;i<num_particles;i++)
	{
		Particle particle_temp;
		particle_temp.id=i;
		particle_temp.x=dist_x(gen);
		particle_temp.y=dist_y(gen);
		particle_temp.theta=dist_psi(gen);
		particle_temp.weight=weight_initial;
		
		particles.push_back(particle_temp);
		weights.push_back(weight_initial);
	}
			
	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//add random Gaussian noise to control input
	random_device rd;
	default_random_engine gen(rd());

	//set the siagma of velocity control input noise 
	double sigma_velocity=0.0001;
	//set the siagma of yaw_rate control input noise 
	double sigma_yaw_rate=0.0001;
	normal_distribution<double> dist_velocity(velocity, sigma_velocity);
    normal_distribution<double>dist_yaw_rate(yaw_rate,sigma_yaw_rate);	
	double velocity_noised=dist_velocity(gen);
	double yaw_rate_noised=dist_yaw_rate(gen);
	
	//add measurements to each particle

	for (int i=0;i<num_particles;i++)
	{
		if(fabs(yaw_rate_noised)!=0.0)
		
		{
			particles[i].x+=velocity_noised/yaw_rate_noised*(sin(particles[i].theta +delta_t*yaw_rate_noised)-sin(particles[i].theta));
			particles[i].y+=velocity_noised/yaw_rate_noised*(cos(particles[i].theta)-cos(particles[i].theta + delta_t*yaw_rate_noised));
			particles[i].theta+= delta_t*yaw_rate_noised;
		}
		else   //go straight
		{
			particles[i].x+=velocity_noised*cos(particles[i].theta);
			particles[i].y+=velocity_noised*sin(particles[i].theta);
			//particles[i].theta= particles[i].theta;
		}
	}

/*
		// add noise to prediction by using std_pos  instead
		
		for (int i=0;i<num_particles;i++)
	{
		if(fabs(yaw_rate)>0.000001)
		
		{
			particles[i].x+=velocity/yaw_rate*(sin(particles[i].theta +delta_t*yaw_rate)-sin(particles[i].theta));
			particles[i].y+=velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta + delta_t*yaw_rate));
			particles[i].theta+= delta_t*yaw_rate;
		}
		else   //go straight
		{
			particles[i].x+=velocity*cos(particles[i].theta);
			particles[i].y+=velocity*sin(particles[i].theta);
			//particles[i].theta= particles[i].theta;
		}
		
		std::normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		std::normal_distribution<double> dist_theta(particles[i].theta,std_pos[2]);

		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}
*/
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(int i=0;i<observations.size();i++)
	{
		double min_dist=100.0;
		double dist;
		int id_matched;
		for(int j=0;j<predicted.size();j++)
		{
			dist=sqrt(pow((observations[i].x-predicted[j].x),2)+pow((observations[i].y-predicted[j].y),2));
			if(dist<min_dist)
			{
				min_dist=dist;
				id_matched=predicted[j].id;
				//id_matched=j;
			}	
		}
		observations[i].id=id_matched;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	
	// 0) update weights for every particle
    double w_sum=0.0;

	for(int i=0;i<particles.size();i++)
	{
		std::vector<LandmarkObs> observations_global;
		//std::vector<LandmarkObs> valid_observations;
        std::vector<LandmarkObs> landmarks;
		LandmarkObs tf;
		double dist;
		double w;

		// 1) for every particle ,tranform the local coordonees measurements to global coordonee measurements

		for(int j=0;j<observations.size();j++)
		{
			dist=sqrt(pow((observations[j].x),2)+pow((observations[j].y),2));
			if(dist<sensor_range)
			{
				//tf.id=i;
				tf.x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x ;
				tf.y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
				observations_global.push_back(tf);
				//valid_observations.push_back(observations[j]);
			}
		}

        for (int j =0; j < map_landmarks.landmark_list.size(); j++) {
			dist=sqrt(pow((map_landmarks.landmark_list[j].x_f-particles[i].x),2)+pow((map_landmarks.landmark_list[j].y_f-particles[i].y),2));
			if(dist<sensor_range)
			{
				tf.id = map_landmarks.landmark_list[j].id_i;
				tf.x= map_landmarks.landmark_list[j].x_f;
				tf.y= map_landmarks.landmark_list[j].y_f;
				landmarks.push_back(tf);
			}
        }
		// 2) for every particle,make the measurement data association to the landmarks
		dataAssociation(landmarks,observations_global);
		
		// 3)for every particle ,calculate the weights 
		w=1.0;
		for(int k=0;k<observations_global.size();k++)
		{ 
			//for better readable 
			int id=observations_global[k].id-1;
			double mu_x=map_landmarks.landmark_list[id].x_f;
			double mu_y=map_landmarks.landmark_list[id].y_f;
			double x=observations_global[k].x;
			double y=observations_global[k].y;
			
			w*=exp(-0.5*(pow((x-mu_x),2)/pow(std_landmark[0],2)+pow((y-mu_y),2)/pow(std_landmark[1],2)))/sqrt(2.0 * M_PI * std_landmark[0] * std_landmark[1]);
		}
	w_sum+=w;
	 weights[i] = w;
     particles[i].weight = w;	
	}
	/*
	// 4) normalized weights
	for(int i=0;i<particles.size();i++)
	{
		weights[i] /= w_sum;
		particles[i].weight /= w_sum;
	}
	*/
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	random_device rd;
	default_random_engine gen(rd());
	std::vector<Particle> resampled_particles;
    std::discrete_distribution<int> resample_id(weights.begin(), weights.end());

    for (int i = 0; i < particles.size(); i++) {
        int id = resample_id(gen);
        resampled_particles.push_back(particles[id]);
    }

    particles = resampled_particles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
