/*
 * Title:        CloudSim Toolkit
 * Description:  CloudSim (Cloud Simulation) Toolkit for Modeling and Simulation of Clouds
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2009-2012, The University of Melbourne, Australia
 */

package org.cloudbus.cloudsim;

/**
 * BinaryPSO implements the popular population-based metaheuristic algorithm 
 * for scheduling cloudlets on VMs.
 * 
 * @author Matt Levan 
 * @since CloudSim Toolkit 2.0
 */
public class BinaryPSO {
    /* List of Vms for submission to cloud resources. */
    List<Vm> vmList;
    
    /* List of cloudlets for submission to cloud resources. */
    List<Cloudlet> cloudletList;

    /* Particle swarm for evaluating resources and finding a Vm 
     * solution for each cloudlet. */
    List<Particle> swarm;

    /* Global best solution. */
    double g;

    /* Inertial constant. */
    double w;

    /* Cognitive constant. */
    double c1;

    /* Social constant. */
    double c2;

    /* Uniform random number. */
    double r;

    /* Random number one. */
    double r1;

    /* Random number two. */
    double r2;

    /* Random Object. */
    Random random;


    /*
     * Constructor accepts a list of Vms and cloudlets.
     *
     */
    public BinaryPSO(List<Vm> vmList, List<Cloudlet> cloudletList) {
        this.vmList = vmList;
        this.coudletList = cloudletList;
        initSwarm(100); // Initialize swarm with 100 particles.
        g = Double.POSITIVE_INFINITY; // Initialize global best.
        c1, c2 = 1.49445; // Initialize cognitive and social constants.
        w = 2.0; // Set the inertial weight.
        random = new Random();
        r = random.nextDouble(); // Initialize random number.
    }


    /**
     * The fitness function calculates the execution times of all possible 
     * cloudlet combinations on every cloud resource and then returns the 
     * highest execution time as the fitness value of the particle, where 
     * exc(Vm vm, List<Cloudlet> cloudletList) is the execution time of 
     * running the cloudletList on vm.
     *
     * @pre $none
     * @post $none
     */
    protected double calculateFitnessValue(Vm vm, List<Cloudlet>) {
         
    }

    /** Initialize the swarm with n particles. */
    protected void initSwarm(int n) {
        this.swarm = new List<Particle>;

        for (int i = 0; i < n; i++) {
            Particle p = new Particle(cloudletList);
            swarm.add(p);
        }
    }

    /** Evaluate a given solution. */
    protected void evaluateSolution(double solution) {
        if (solution < g) {
            g = solution;
        }
    }

    /** Run the particle swarm optimization algorithm. 
     *  
     *  @pre $none
     *  @post $none
     *  @return vmIds List of Vm ids to match cloudlets with. */
    protected List<Integer> run {
        
    }
}


/**
 * Particle class implements the particle for building a swarm.
 */
class Particle {
    /* List of cloudlets. */
    List<Cloudlet> cloudletList;

    /* Velocity. */
    public double v;

    /* Personal best solution. */
    public double pBest;

    /* Personal current solution. */
    public double pCurrent;

    /* Constructor. */
    public Particle(List<Cloudlet> cloudletList) {
        this.cloudletList = cloudletList;
        this.v = 0.0;
        this.pCurrent = 0.0;
        this.pBest = 0.0;
    }

    /** Calculate velocity. */
    public void calcVelocity(double w, double c1, double c2, 
                             double r1, double r2, double g) {
        v = w*v+c1*r1*(pBest-pCurrent)+c2*r2*(g-pCurrent);
    }

    /** Calculate personal current solution. */
    public void calcpCurrent() {
        pCurrent = pCurrent+v;
    }

    /** Sigmoid function. */
    public double s(double pCurrent) {
        return 1/(1+Math.exp(-1*pCurrent));
    }

    /** Calculate position. */
    public double calcPosition(double s) {
        if (s < 
    }
}
