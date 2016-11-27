package org.cloudbus.cloudsim;

/**
 * Particle class implements the particle for building a swarm.
 *
 * TODO: How exactly is the fitness calculated? Does each particle calculate 
 * execution times of the entire power set of cloudlets on only a single Vm or 
 * on every Vm? 
 *
 * TODO: Also, does each particle only calculate execution times of a 
 * partition of the power set of cloudlets?
 *
 * TODO: How does the particle express the personal current and best 
 * solutions (in terms of pairings)?
 */

import java.util.*;

class Particle {
    /* List of cloudlets. */
    List<Cloudlet> cloudletList;

    /* D-dimensional space for positions. */
    ArrayList<int[]> position;

    /* D-dimensional space for velocities. */
    ArrayList<double[]> velocity;

    /* Personal best solution. */
    ArrayList<int[]> bestPosition;

    /* Fitness value. */
    double fitness;

    /* Best fitness value. */
    double bestFitness;

    /* Uniform random number from swarm. */
    int r;

    /* Random number one. */
    int r1;

    /* Random number two. */
    int r2;

    /* Random Object. */
    Random random;

    /**
     * Particle constructor. 
     *
     * @param cloudletList List of cloudlets to be considered.
     * @param n Number of Vms available for scheduling.
     */
    public Particle(ArrayList<int[]> position, double fitness, 
            ArrayList<double[]> velocity, ArrayList<int[]> bestPosition, 
            double bestFitness) {
        this.position = position;
        this.fitness = fitness;
        this.velocity = velocity;
        this.bestPosition = bestPosition;
        this.bestFitness = bestFitness;
    }
}
