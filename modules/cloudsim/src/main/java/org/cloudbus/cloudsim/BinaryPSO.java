/*
 * Title:        CloudSim Toolkit
 * Description:  CloudSim (Cloud Simulation) Toolkit for Modeling and Simulation of Clouds
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2009-2012, The University of Melbourne, Australia
 */

package org.cloudbus.cloudsim;

import java.util.*;

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

    /* Power set of list of cloudlets. */
    List<List<Cloudlet>> ps;

    /* Particle swarm for evaluating resources and finding a Vm 
     * solution for each cloudlet. */
    List<Particle> swarm;

    /* List of all possible solutions. */
    List<ArrayList<ArrayList<Integer>>> solutions;

    /* Global best solution. */
    ArrayList<double[]> g;

    /* Inertial constant. */
    double w;

    /* Cognitive constant. */
    final double c1;

    /* Social constant. */
    final double c2;

    /* Uniform random number. */
    final double r;

    /* Random Object. */
    Random random;

    /* Number of iterations. */
    final int numIterations;

    /* Number of particles. */
    final int numParticles;

    /* Inertia calculation technique choice. */
    final int inertiaTechnique;

    /* Fixed inertia weight. */
    final double fixedIntertiaWeight = 0.5;

    /* Solution counter for use in the calcSolutions method. */
    int iteration = 0;


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
        this.ps = powerSet(cloudletList);
        this.solutions = new ArrayList<ArrayList<ArrayList<Integer>>>();
        initSwarm(); // Initialize swarm with 100 particles.
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
     * Fitness for a single PSO particle. 
     *
     * The fitness function calculates the execution times of all possible 
     * cloudlet combinations on every cloud resource and then returns the 
     * highest execution time as the fitness value of the particle, where 
     * exc(Vm vm, List<Cloudlet> cloudletList) is the execution time of 
     * running the cloudletList on vm.
     */
    protected double calculateFitnessValue() {
        List<Double> executionTimes = new ArrayList<Double>();
        
        for (List<Cloudlet> subset : ps) {
            for (Vm vm : vmList) {
                executionTimes.add(calcExecutionTime(vm, subset));
            }
        }

        return Collections.max(executionTimes);
    }

    /**
     * Calculate execution time of the power set of cloudlets on a single Vm.
     *
     * @param vm Vm to calculate execution time on.
     * @return Execution time of all cloudlets on vm.
     */
    private double calcExecutionTime(Vm vm, List<Cloudlet> subset) {
        double total = 0.0;

        for (Cloudlet cloudlet : subset) {
            total += cloudlet.getCloudletLength() / 
                vm.getHost().getTotalAllocatedMipsForVm(vm);
        }
        
        return total;
    }

    /** 
     * Initialize the swarm with numParticles particles. 
     * 
     * @TODO: Write code to initialize first position, velocity matrices,
     * and fitness and use those values to instantiate all particles in 
     * the swarm.
     */
    protected void initSwarm() {
        swarm = new ArrayList<Particle>();

        for (int i = 0; i < numParticles; i++) {
            swarm.add(new Particle(position, fitness, velocity, 
                    bestPosition, bestFitness));
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

    /**
     * Generates the search space by calculating all possible solutions as 
     * a function of the cloudletList and vmList.
     */
    protected void genSearchSpace() {
        /* Generate a list of cloudletId integers to pass to calcSolutions. */
        ArrayList<Integer> cloudletIds = new ArrayList<Integer>();
        for (Cloudlet cloudlet : cloudletList) {
            cloudletIds.add(cloudlet.getCloudletId());
        }

        /* Instantiate a new list that stores calculations for the 
         * recursive calcSolutions() to use. */ 
        ArrayList<ArrayList<Integer>> pastCalculations = 
            new ArrayList<ArrayList<Integer>>();

        /* Run the recursive helper function. */
        calcSolutions(cloudletIds, pastCalculations, iteration);
    }

    /**
     * Calculates all possible solutions as a function of the cloudletList and 
     * vmList (helper function for genSearchSpace())).
     * 
     * @param remainingList List of remaining cloudlets.
     * @param pastCalculations List of previously calculated solutions.
     * @param iteration Iteration counter.
     */
    protected void calcSolutions(ArrayList<Integer> remainingList, 
            ArrayList<ArrayList<Integer>> pastCalculations, int iteration) {
        List<List<Integer>> ps = powerSet(remainingList);

        if (iteration >= cloudletList.size()-2) {
            for (List<Integer> x : ps) {
                ArrayList<ArrayList<Integer>> t = 
                    (ArrayList<ArrayList<Integer>>) pastCalculations.clone();
                ArrayList<Integer> remaining = new 
                    ArrayList<Integer>(remainingList);
                remaining.removeall(x);
                t.add((ArrayList<Integer>) x);
                t.add(remaining);
                solutions.add(t);
            }

            /* Reset the iteration counter. */
            iteration = 0;
        }
        else {
            for (List<Integer> x : ps) {
                ArrayList<ArrayList<Integer>> t = 
                    (ArrayList<ArrayList<Integer>>) pastCalculations.clone();
                ArrayList<Integer> remaining = new 
                    ArraList<Integer>(remainingList);
                t.add((ArrayList<Integer>) x);
                remaining.removeAll(x);
                calcSolutions((ArrayList<Integer>) remaining, t, iteration+1);
            }
        }
    }

    /** 
     * Generate power set. 
     *
     * @param list Collection<T>.
     * @return Power set of the list.
     */
	public static <T> List<List<T>> powerSet(Collection<T> list) {
	    List<List<T>> ps = new ArrayList<List<T>>();
	    ps.add(new ArrayList<T>());   // Add the empty set.
	 
	    // For every item in the original list...
	    for (T item : list) {
	        List<List<T>> newPs = new ArrayList<List<T>>();
	 
	        for (List<T> subset : ps) {
	            // Copy all of the current powerset's subsets.
	            newPs.add(subset);
	 
	            // Plus the subsets appended with the current item.
	            List<T> newSubset = new ArrayList<T>(subset);
	            newSubset.add(item);
	            newPs.add(newSubset);
	        }
	 
	        ps = newPs;
	    }
	    return ps;
	}
}
