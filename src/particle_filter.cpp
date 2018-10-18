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

    num_particles = 100;
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
        pa.weight = 1.0;
        particles.push_back(pa);
        weights.push_back(pa.weight);
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

    default_random_engine gen;

    for (int i = 0; i < num_particles; i++) {
        auto pa = particles[i];
        if (fabs(yaw_rate) < 0.0001) {
            pa.x += velocity * delta_t * cos(pa.theta);
            pa.y += velocity * delta_t * sin(pa.theta);
        } else {
            pa.x += velocity / yaw_rate *
                    (sin(pa.theta + yaw_rate * delta_t) - sin(pa.theta));
            pa.y += velocity / yaw_rate *
                    (cos(pa.theta) - cos(pa.theta + yaw_rate * delta_t));
            pa.theta += yaw_rate * delta_t;
        }

        normal_distribution<double> dist_x(pa.x, std_pos[0]);
        normal_distribution<double> dist_y(pa.y, std_pos[1]);
        normal_distribution<double> dist_theta(pa.theta, std_pos[2]);

        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
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

    // for normalizing weights
    double weight_normalizer = 0.0;

    for (int i = 0; i < num_particles; i++) {
        double px = particles[i].x;
        double py = particles[i].y;
        double ptheta = particles[i].theta;

        vector<LandmarkObs> transformed_observations;

        // Observations coordinates Transformation from vehicle to map
        for (int j = 0; j < observations.size(); j++) {
            int oid = observations[j].id;
            double ox = px + (cos(ptheta) * observations[j].x) -
                        (sin(ptheta) * observations[j].y);
            double oy = py + (sin(ptheta) * observations[j].x) +
                        (cos(ptheta) * observations[j].y);
            transformed_observations.push_back(LandmarkObs{oid, ox, oy});
        }

        // Map ROI filter
        vector<LandmarkObs> predicted_landmarks;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            int mid = map_landmarks.landmark_list[j].id_i;
            double mx = map_landmarks.landmark_list[j].x_f;
            double my = map_landmarks.landmark_list[j].y_f;

            if (fabs(mx - px) <= sensor_range &&
                fabs(my - py) <= sensor_range) {
                predicted_landmarks.push_back(LandmarkObs{mid, mx, my});
            }
        }

        // Data association
        dataAssociation(predicted_landmarks, transformed_observations);

        // Calculate the weights with Multivariate Gaussian distribution.
        double sigma_x = std_landmark[0];
        double sigma_y = std_landmark[1];
        double sigma_x_2 = pow(sigma_x, 2);
        double sigma_y_2 = pow(sigma_y, 2);
        double normalizer = (1.0 / (2.0 * M_PI * sigma_x * sigma_y));

        /*Calculate the weight of particle based on the multivariate Gaussian
         * probability function*/
        for (int k = 0; k < transformed_observations.size(); k++) {
            double trans_obs_x = transformed_observations[k].x;
            double trans_obs_y = transformed_observations[k].y;
            double trans_obs_id = transformed_observations[k].id;
            double multi_prob = 1.0;

            for (int l = 0; l < predicted_landmarks.size(); l++) {
                double pred_landmark_x = predicted_landmarks[l].x;
                double pred_landmark_y = predicted_landmarks[l].y;
                double pred_landmark_id = predicted_landmarks[l].id;

                if (trans_obs_id == pred_landmark_id) {
                    multi_prob =
                        normalizer *
                        exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2) /
                                     (2.0 * sigma_x_2)) +
                                    (pow((trans_obs_y - pred_landmark_y), 2) /
                                     (2.0 * sigma_y_2))));
                    particles[i].weight *= multi_prob;
                }
            }
        }
        weight_normalizer += particles[i].weight;
    }

    // Normalize the weights
    for (int i = 0; i < particles.size(); i++) {
        particles[i].weight /= weight_normalizer;
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional
    // to their weight. NOTE: You may find std::discrete_distribution helpful
    // here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> resampled_particles;
    default_random_engine gen;

    // generate random starting index for resampling wheel
    uniform_int_distribution<int> p_index(0, num_particles - 1);
    auto index = p_index(gen);

    double beta = 0.0;

    // get max weight
    double max_weight = 2.0 * *max_element(weights.begin(), weights.end());

    // spin the resample wheel!
    for (int i = 0; i < num_particles; i++) {
        uniform_real_distribution<double> random_weight(0.0, max_weight);
        beta += random_weight(gen);

        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        resampled_particles.push_back(particles[index]);
    }

    particles = resampled_particles;
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
