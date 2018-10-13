/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first
    // position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method
    // (and others in this file).

    num_particles = 3;
    // TODO: Set standard deviations for x, y, and theta
    default_random_engine gen;

    // line creates a normal (Gaussian) distribution
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 0; i < num_particles; ++i) {
        Particle pa;
        pa.id = i;
        pa.x = dist_x(gen);
        pa.y = dist_y(gen);
        pa.theta = dist_theta(gen);
        pa.weight = 0;
        particles.push_back(pa);

        // Print your samples to the terminal.
        cout << "Sample " << i + 1 << " " << pa.x << " " << pa.y << " "
             << pa.theta << endl;
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and
    // std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    for (int i = 0; i < num_particles; i++) {
        Particle& pa = particles[i];
        if (fabs(yaw_rate) < 0.00001) {
            pa.x += velocity * delta_t * cos(pa.theta);
            pa.y += velocity * delta_t * sin(pa.theta);
        } else {
            pa.x += +velocity / yaw_rate *
                    (sin(pa.theta + yaw_rate * delta_t) - sin(pa.theta));
            pa.y += +velocity / yaw_rate *
                    (cos(pa.theta) - cos(pa.theta + yaw_rate * delta_t));
            pa.theta += yaw_rate * delta_t;
        }
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
                                     std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed
    // measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will
    // probably find it useful to
    //   implement this method and use it as a helper during the updateWeights
    //   phase.
    for (int i = 0; i < observations.size(); i++) {
        LandmarkObs& o = observations[i];

        int association_id = -1;
        double min_dist = numeric_limits<double>::max();

        for (int j = 0; j < predicted.size(); j++) {
            LandmarkObs p = predicted[j];

            // cal the dist between predicted and observations
            double cur_dist =
                sqrt((o.x - p.x) * (o.x - p.x) + (o.y - p.y) * (o.y - p.y));
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                association_id = p.id;
            }
        }
        o.id = association_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs>& observations,
                                   const Map& map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian
    // distribution. You can read
    //   more about this distribution here:
    //   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your
    // particles are located
    //   according to the MAP'S coordinate system. You will need to transform
    //   between the two systems. Keep in mind that this transformation requires
    //   both rotation AND translation (but no scaling). The following is a good
    //   resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to
    //   implement (look at equation 3.33
    //   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional
    // to their weight. NOTE: You may find std::discrete_distribution helpful
    // here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
}

Particle ParticleFilter::SetAssociations(Particle& particle,
                                         const std::vector<int>& associations,
                                         const std::vector<double>& sense_x,
                                         const std::vector<double>& sense_y) {
    // particle: the particle to assign each listed association, and
    // association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed
    // association sense_x: the associations x mapping already converted to
    // world coordinates sense_y: the associations y mapping already converted
    // to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
