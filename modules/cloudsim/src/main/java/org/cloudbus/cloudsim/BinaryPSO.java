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
 * @author Amlan Chatterjee, Crosby Lanham, Matt Levan, Mishael Zerrudo
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

    /* Random Object. */
    Random random;

    /* Number of iterations. */
    int numIterations;

    /* Number of particles. */
    int numParticles;

    /* Inertia calculation technique choice. */
    int inertiaTechnique;

    /* Fixed inertia weight. */
    int fixedIntertiaWeight = 0.5;


    /**
     * Constructor accepts a list of Vms and cloudlets.
     *
     * @param vmList List of Vms Objects for consideration.
     * @param cloudletList List of cloudlet jobs to be run.
     * @param numIterations Number of iterations.
     * @param numParticles Number of particles.
     * @param inertiaTechnique Integer to specify which technique to use for 
     * calculating the inertia:
     * <ol start="0">
     *  <li>Fixed Inertia Weight (FIW)</li>
     *  <li>Random Inertia Weight (RIW)</li>
     *  <li>Linearly Decreasing Inertia Weight (LDIW)</li>
     *  <li>Combination</li>
     * </ol>
     */
    public BinaryPSO(List<Vm> vmList, List<Cloudlet> cloudletList, 
                     int numIterations, int numParticles, int inertiaTechnique) {
        this.vmList = vmList;
        this.coudletList = cloudletList;
        this.numIterations = numIterations;
        this.numParticles = numParticles;
        this.inertiaTechnique = inertiaTechnique;
        initSwarm(numParticles); // Initialize swarm with 100 particles.
        g = Double.POSITIVE_INFINITY; // Initialize global best.
        c1 = 1.49445; // Initialize cognitive constant.
        c2 = 1.49445; // Initialize social constants.
        w = 2.0; // Set the inertial weight.
        random = new Random();
        r = random.nextDouble(); // Initialize random number.
    }

    /**
     * The inertia weight w can be calculated with a variety of techniques:
     *
     * @todo:   Implement RIW, LDIW, and Combination techniques.
     * @param   numIterations Number of iterations.
     * @param   numParticles Number of individual particles.
     * @return  Sets the global inertia value w.
     */
    protected void calculateInertia(int numIterations, int numParticles
                                    int inertiaTechnique) {
        switch (inertiaTechnique) {
            /* Fixed Inertia Weight (FIW). */
            case 0: 
                w = fixedInertiaWeight;
                break;
            /* Random Inertia Weight (RIW). */
            case 1: 
                break;
            /* Linearly Reducing Inertia Weight (LDIW). */
            case 2: 
                break;
            /* Combination. */
            case 3: 
                break;
            default:
                break;
        }
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
    protected double calculateFitnessValue(Vm vm, List<Cloudlet> cloudletList) {
         
    }

    /** Initialize the swarm with n particles. */
    protected void initSwarm(int n) {
        swarm = new List<Particle>();

        for (int i = 0; i < n; i++) {
            Particle p = new Particle(cloudletList);
            swarm.add(p);
        }
    }

    /**
     * Evaluate a given solution. 
     *
     * @param solution Don't know yet.
     */
    protected void evaluateSolution(double solution) {
        if (solution < g) {
            g = solution;
        }
    }

    /** 
     * Run the particle swarm optimization algorithm. 
     *  
     *  @pre $none
     *  @post $none
     *  @return vmIds List of Vm ids to match cloudlets with. */
    protected List<Integer> run() {
        List<Integer> vmIds = new List<Integer>();

        for (int i = 0; i < numIterations; i++) {
            for (int j = 0; j < numParticles; j++) {

            }
        }
    }

    /** 
     * Set fixed inertia weight (FIW).
     *
     * @param fixedInertiaWeight The desired value.
     */
    protected void setFixedInertiaWeight(int fixedInertiaWeight) {
        this.fixedInertiaWeight = fixedInertiaWeight;
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

    /* Random number one. */
    double r1;

    /* Random number two. */
    double r2;

    /* Random Object. */
    Random random;

    /* Constructor. */
    public Particle(List<Cloudlet> cloudletList) {
        cloudletList = cloudletList;
        v = 0.0;
        pCurrent = 0.0;
        pBest = 0.0;
        random = new Random();
    }

    /** Calculate velocity. */
    public void calcVelocity(double w, double c1, double c2, double g) {
        r1 = random.nextDouble();
        r2 = random.nextDouble();
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
    public double calcPosition() {
        if (s(pCurrent) <= r)
            pCurrent = 0;
        else
            pCurrent = 1;
    }

    /** 
     * Generate power set of cloudletList. 
     *
     * @todo: Write this code.
     */
}
