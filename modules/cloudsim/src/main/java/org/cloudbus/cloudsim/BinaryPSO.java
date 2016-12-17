/*
 * Title:        CloudSim Toolkit
 * Description:  CloudSim (Cloud Simulation) Toolkit for Modeling and 
 * Simulation of Clouds
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
    /** Uniform rand number. */
    private final double r;
    /** Number of iterations. */
    private final int numIterations;
    /** Number of particles. */
    private final int numParticles;
    /** Inertia calculation technique choice. */
    private final int inertiaTechnique;
    /**
     * Cognitive constant.
     */
    protected double c1;
    /**
     * Social constant.
     */
    protected double c2;
    /**
     * Fixed inertia weight.
     */
    protected double fixedInertiaWeight = 0.5;
    /** List of Vms for submission to cloud resources. */
    private List<Vm> vmList;
    /** List of cloudlets for submission to cloud resources. */
    private List<Cloudlet> cloudletList;
    /**
     * Particle swarm for evaluating resources and finding a Vm
     * solution for each cloudlet.
     */
    private List<Particle> swarm;
    /** List of all possible solutions. */
    private List<ArrayList<ArrayList<Integer>>> solutions;
    /** Global best solution. */
    private ArrayList<int[]> globalBest;
    /** Average fitness matrix. */
    private ArrayList<double[]> averageFitnesses;
    /** Global best fitness. */
    private double globalBestFitness;
    /** Run times. */
    private ArrayList<double[]> runTime;
    /** Number of Vms. */
    private int n;
    /** Number of cloudlets. */
    private int m;
    /** Random Object. */
    private Random rand;
    /** Solution counter for use in the calcSolutions method. */
    private int iteration = 0;


    /**
     * Constructor accepts a list of Vms and cloudlets.
     *
     * @param vmList           List of Vms Objects for consideration.
     * @param cloudletList     List of cloudlet jobs to be run.
     * @param numIterations    Number of iterations.
     * @param numParticles     Number of particles.
     * @param inertiaTechnique Integer to specify which technique to use for
     *                         calculating the inertia:
     *                         <ol start="0">
     *                         <li>Fixed Inertia Weight (FIW)</li>
     *                         <li>Random Inertia Weight (RIW)</li>
     *                         <li>Linearly Decreasing Inertia Weight (LDIW)</li>
     *                         <li>Combination</li>
     *                         </ol>
     */
    public BinaryPSO(List<Vm> vmList, List<Cloudlet> cloudletList,
                     int numIterations, int numParticles, int inertiaTechnique) {
        this.vmList = vmList;
        this.cloudletList = cloudletList;
        this.numIterations = numIterations;
        this.numParticles = numParticles;
        this.inertiaTechnique = inertiaTechnique;
        this.solutions = new ArrayList<ArrayList<ArrayList<Integer>>>();
        this.n = vmList.size();
        this.m = cloudletList.size();
        this.globalBestFitness = Double.MAX_VALUE; // Small values are better.
        this.globalBest = new ArrayList<int[]>();
        c1 = 1.49445; // Initialize cognitive constant.
        c2 = 1.49445; // Initialize social constants.
        rand = new Random();
        r = rand.nextDouble(); // Initialize rand number.
    }


    /**
     * Fitness for a single PSO particle.
     * <p>
     * The fitness function calculates the execution times of all possible
     * cloudlet combinations on every cloud resource and then returns the
     * highest execution time as the fitness value of the particle, where
     * exc(Vm vm, List<Cloudlet> cloudletList) is the execution time of
     * running the cloudletList on vm.
     *
     * @param runTime   The runtime matrix.
     * @param positions The positions matrix.
     * @return The double fitness value.
     */
    private double calcFitness(ArrayList<double[]> runTime,
                               ArrayList<int[]> positions) {
        double[] sum = new double[n];

        for (int i = 0; i < n; i++) {
            double[] time = runTime.get(i);
            int[] pos = positions.get(i);

            for (int j = 0; j < pos.length; j++) {
                if (pos[j] == 1)
                    sum[i] += time[j];
            }
        }

        /* Find the longest execution time. */
        double result = 0;

        for (int i = 0; i < n; i++) {
            if (result < sum[i])
                result = sum[i];
        }

        return result;
    }

    /**
     * The inertia weight w can be calculated with a variety of techniques.
     *
     * @param particleNumber   Particle number.
     * @param averageFitnesses Average particle fitnesses matrix.
     * @return The global inertia value w.
     */
    private double calcInertia(int particleNumber,
                               ArrayList<double[]> averageFitnesses) {
        int k = 5;
        double w = 0.0;

        double w_max = 0.9;
        double w_min = 0.1;

        double t_max = numIterations;
        double t = iteration;

        if (t % k == 0 && t != 0) {
            double p = 0.0; // Annealing probability.

            double currentFitness =
                    averageFitnesses.get(particleNumber)[iteration];
            double previousFitness =
                    averageFitnesses.get(particleNumber)[iteration - k];

            if (previousFitness <= currentFitness) {
                p = 1;
            } else {
                /* Annealing temperature. */
                double coolingTemp_Tt = 0.0;

                Particle currParticle = swarm.get(particleNumber);
                double bestFitness = currParticle.bestFitness;

                double particleFitnessAverage = 0;

                int counter = 0;
                for (int i = 0; i < iteration; i++) {
                    if (averageFitnesses.get(particleNumber)[i] > 0) {
                        particleFitnessAverage +=
                                averageFitnesses.get(particleNumber)[i];
                        counter++;
                    }
                }

                particleFitnessAverage = particleFitnessAverage / counter;
                coolingTemp_Tt = particleFitnessAverage / bestFitness - 1;
                p = Math.exp(-(previousFitness - currentFitness) / coolingTemp_Tt);
            }

            int random = rand.nextInt(2);

            /* New inertia weight. */
            if (p >= random) {
                w = random / 2 + 1;
            } else {
                w = random / 2;
            }
        } else {
            /* New inertia weight using LDIW. */
            double w_fraction = (w_max - w_min) * (t_max - t) / t_max;
            w = w_max - w_fraction;
        }

        switch (inertiaTechnique) {
            /* Fixed Inertia Weight (FIW). */
            case 0:
                w = fixedInertiaWeight;
                break;
            /* Random Inertia Weight (RIW). */
            case 1:
                w = randomInertiaWeight(particleNumber, averageFitnesses);
                break;
            /* Linearly Reducing Inertia Weight (LDIW). */
            case 2:
                w = linearlyDecreasingInertiaWeight(particleNumber,
                        averageFitnesses);
                break;
            /* Combination. */
            case 3:
                w = combinationInertiaWeight(particleNumber, averageFitnesses);
                break;
        }

        return w;
    }

    /**
     * Calculate the initial positions matrix.
     *
     * @return Initial positions, ensured for complete jobs assignment.
     */
    private ArrayList<int[]> calcInitPositions() {
        ArrayList<int[]> initPositions = new ArrayList<int[]>();
        int[] assignedTasks = new int[m];

        /* Iterate through the VMs. */
        for (int i = 0; i < n; i++) {
            int[] randomPositions = new int[m];

            /* Iterate through the cloudlets. */
            for (int j = 0; j < m; j++) {
                if (assignedTasks[j] == 0) {
                    randomPositions[j] = rand.nextInt(2);

                    if (randomPositions[j] == 1) {
                        assignedTasks[j] = 1;
                    }
                } else {
                    randomPositions[j] = 0;
                }
            }

            initPositions.add(randomPositions);
        }

        /* Ensure that all jobs have been assigned to at least one VM. */
        ArrayList<int[]> newPositions = ensureCompleteAssignment(initPositions,
                assignedTasks);

        return newPositions;
    }

    /**
     * Calculate the initial velocities matrix.
     *
     * @return Initial velocities, ensured for complete jobs assignment.
     */
    private ArrayList<double[]> calcInitVelocities() {
        ArrayList<double[]> initVelocities = new ArrayList<double[]>();
        int[] assignedTasks = new int[m];

        /* Iterate through the VMs. */
        for (int i = 0; i < n; i++) {
            double[] randomPositions = new double[m];

            /* Iterate through the cloudlets. */
            for (int j = 0; j < m; j++) {
                if (assignedTasks[j] == 0) {
                    randomPositions[j] = rand.nextInt(2);

                    if (randomPositions[j] == 1)
                        assignedTasks[j] = 1;
                } else {
                    randomPositions[j] = 0;
                }
            }

            initVelocities.add(randomPositions);
        }

        return initVelocities;
    }

    /**
     * Calculate new positions.
     *
     * @param w Inertia weight for the particle.
     * @param p Particle.
     */
    private void calcNewPositions(double w, Particle p) {
        ArrayList<int[]> newPositionsMatrix = new ArrayList<int[]>();
        int[] assignedTasksArrayInPositionsMatrix = new int[m];

        /* For each VM. */
        for (int i = 0; i < p.velocity.size(); i++) {
            double[] vmVelocities = p.velocity.get(i);

            int[] newPosition = new int[m];

            /* For each cloudlet. */
            for (int j = 0; j < vmVelocities.length - 1; j++) {
                p.r = rand.nextInt(2);

                if (assignedTasksArrayInPositionsMatrix[j] == 0) {
                    double s = s(vmVelocities[j]);

                    if (s > p.r) {
                        newPosition[j] = 1;
                    } else {
                        newPosition[j] = 0;
                    }

                    if (newPosition[j] == 1) {
                        assignedTasksArrayInPositionsMatrix[j] = 1;
                    }
                } else {
                    newPosition[j] = 0;
                }
            }

            /* Add new velocities to the ArrayList. */
            newPositionsMatrix.add(newPosition);
        }

        /* Update the particle's position. */
        newPositionsMatrix = ensureCompleteAssignment(newPositionsMatrix,
                assignedTasksArrayInPositionsMatrix);
        newPositionsMatrix = rebalancePSO(newPositionsMatrix, runTime);
        p.position = newPositionsMatrix;
    }

    /**
     * Calculate new velocities.
     *
     * @param w Inertia weight for the particle.
     * @param p Particle.
     */
    private void calcNewVelocities(double w, Particle p) {
        ArrayList<double[]> newVelocitiesMatrix = new ArrayList<double[]>();
        int[] assignedTasksArrayInVelocityMatrix = new int[m];

        for (int i = 0; i < p.velocity.size(); i++) {
            double[] vmVelocities = p.velocity.get(i);
            int[] vmBestPositions = p.bestPosition.get(i);
            int[] vmPositions = p.position.get(i);
            int[] vmGlobalBestPositions = globalBest.get(i);

            double[] newVelocities = new double[m];

            for (int j = 0; j < vmVelocities.length - 1; j++) {
                p.r1 = rand.nextInt(2);
                p.r2 = rand.nextInt(2);

                if (assignedTasksArrayInVelocityMatrix[j] == 0) {
                    /* Velocity vector. */
                    newVelocities[j] = w * vmVelocities[j + 1] + c1 * p.r1 *
                            (vmBestPositions[j] - vmPositions[j] + c2 * p.r2 *
                                    (vmGlobalBestPositions[j] - vmPositions[j]));
                    /* Formula 4. */
                    if (newVelocities[j] < 0) {
                        newVelocities[j] = 0;
                    } else if (newVelocities[j] > 1) {
                        newVelocities[j] = 1;
                    }
                } else {
                    newVelocities[j] = 0;
                }
            }

            /* Add new velocities to the ArrayList. */
            newVelocitiesMatrix.add(newVelocities);
        }

        /* Update the particle's velocity. */
        p.velocity = newVelocitiesMatrix;
    }

    /**
     * Calculate the runtime, which is the execution time of each cloudlet on
     * each VM, separately.
     */
    private ArrayList<double[]> calcRunTime() {
        ArrayList<double[]> runTime = new ArrayList<double[]>();

        /* Iterate through all the VMs. */
        for (int i = 0; i < n; i++) {
            Vm vm = vmList.get(i);

            /* Execution times of running cloudlet j on VM i. */
            double[] times = new double[m];

            /* Iterate through all the cloudlets. */
            for (int j = 0; j < m; j++) {
                Cloudlet cloudlet = cloudletList.get(j);
                times[j] = (double) cloudlet.getCloudletLength() /
                        vm.getHost().getTotalAllocatedMipsForVm(vm);
            }

            runTime.add(times);
        }

        return runTime;
    }

    /**
     * Calculates the inertia weight using the combination technique.
     *
     * @param particleNumber   Particle number.
     * @param averageFitnesses Average particle fitnesses matrix.
     * @return The global inertia value w.
     */
    private double combinationInertiaWeight(int particleNumber,
                               ArrayList<double[]> averageFitnesses) {
        int k = 5;
        double w;

        double w_max = 0.9;
        double w_min = 0.1;

        double t_max = numIterations;
        double t = iteration;

        if (t % k == 0 && t != 0) {
            double p; // Annealing probability.

            double currentFitness =
                    averageFitnesses.get(particleNumber)[iteration];
            double previousFitness =
                    averageFitnesses.get(particleNumber)[iteration - k];

            if (previousFitness <= currentFitness) {
                p = 1;
            } else {
                /* Annealing temperature. */
                double coolingTemp_Tt = 0.0;

                Particle currParticle = swarm.get(particleNumber);
                double bestFitness = currParticle.bestFitness;

                double particleFitnessAverage = 0;

                int counter = 0;
                for (int i = 0; i < iteration; i++) {
                    if (averageFitnesses.get(particleNumber)[i] > 0) {
                        particleFitnessAverage +=
                                averageFitnesses.get(particleNumber)[i];
                        counter++;
                    }
                }

                particleFitnessAverage = particleFitnessAverage / counter;
                coolingTemp_Tt = particleFitnessAverage / bestFitness - 1;
                p = Math.exp(-(previousFitness - currentFitness) / coolingTemp_Tt);
            }

            int random = rand.nextInt(2);

            /* New inertia weight. */
            if (p >= random) {
                w = random / 2 + 1;
            } else {
                w = random / 2;
            }
        } else {
            /* New inertia weight using LDIW. */
            double w_fraction = (w_max - w_min) * (t_max - t) / t_max;
            w = w_max - w_fraction;
        }

        return w;
    }

    /**
     * Ensure that all jobs are assigned to at least one VM in the given
     * positions matrix.
     *
     * @param initPositions The positions matrix to be checked for completeness.
     * @param assignedTasks The list of already assigned tasks.
     * @return New positions matrix.
     */
    private ArrayList<int[]> ensureCompleteAssignment(
            ArrayList<int[]> initPositions, int[] assignedTasks) {
        ArrayList<int[]> newPositions = initPositions;

        /* Check if task is not yet assigned. */
        for (int i = 0; i < assignedTasks.length; i++) {
            if (assignedTasks[i] == 0) {
                /* Pick a rand VM index. */
                int x = rand.nextInt(n);

                int[] positions = newPositions.get(x);
                positions[i] = 1;

                newPositions.set(x, positions);
            }
        }

        return newPositions;
    }

    /**
     * Evaluate a given solution and update the particle and global bests if
     * necessary.
     *
     * @param newFitness      The current particle's fitness.
     * @param currentPosition The current particle's position.
     * @param currentParticle The current particle.
     */
    private void evaluateSolution(double newFitness,
                                  ArrayList<int[]> currentPosition,
                                  Particle currentParticle) {
        /* If the newFitness is better than the current particle's fitness... */
        if (newFitness < currentParticle.fitness) {
            /* Update the particle's fitness. */
            currentParticle.fitness = newFitness;

            /* Update the particle's position. */
            currentParticle.position = currentPosition;
        }

        /* If the newFitness is better than the globalBestFitness... */
        if (newFitness < globalBestFitness) {
            globalBestFitness = newFitness;
            globalBest = currentPosition;
        }
    }

    /**
     * Gets the cloudlet positions and returns them to DatacenterBroker.java.
     *
     * @param globalBest The global best positions matrix.
     * @return The VM positions for each cloudlet to be
     * allocated by the broker.
     */
    private ArrayList<Integer> getCloudletPositions(ArrayList<int[]> globalBest) {
        Integer[] cloudletPositions = new Integer[m];

        /* For each VM. */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (globalBest.get(i)[j] == 1) {
                    cloudletPositions[j] = i;
                }
            }
        }

        // System.out.println(Arrays.toString(cloudletPositions));

        return new ArrayList<Integer>(Arrays.asList(cloudletPositions));
    }

    /**
     * Instantiate average fitness matrix.
     */
    private void initAverageFitnesses() {
        averageFitnesses = new ArrayList<double[]>();

        for (int i = 0; i < numParticles; i++) {
            averageFitnesses.add(new double[numIterations]);
        }
    }

    /**
     * Initialize the swarm with numParticles particles.
     *
     * @TODO: Write code to initialize first position, velocity matrices,
     * and fitness and use those values to instantiate all particles in
     * the swarm.
     */
    private void initSwarm() {
        swarm = new ArrayList<Particle>();
        runTime = calcRunTime();

        /* Initialize each particle in the swarm. */
        for (int i = 0; i < numParticles; i++) {
            /* Calculate the initial positions matrix. */
            ArrayList<int[]> positions = calcInitPositions();

            /* Calculate the initial velocities matrix. */
            ArrayList<double[]> velocities = calcInitVelocities();

            /* Calculate the initial fitness value. */
            double fitness = calcFitness(runTime, positions);

            swarm.add(new Particle(positions, fitness, velocities, positions,
                    fitness));

            /* Initialize global best fitness value. */
            if (swarm.get(i).fitness < globalBestFitness) {
                globalBestFitness = swarm.get(i).fitness;
                globalBest = swarm.get(i).position;
            }
        }

        initAverageFitnesses();
    }

    /**
     * The inertia weight w can be calculated with a variety of techniques:
     * <p>
     * TODO:    Implement RIW, LDIW, and combination techniques.
     *
     * @param particleNumber   Particle number.
     * @param averageFitnesses Average particle fitnesses matrix.
     * @return The global inertia value w.
     */
    private double linearlyDecreasingInertiaWeight(int particleNumber,
                               ArrayList<double[]> averageFitnesses) {
        double w_max = 0.9;
        double w_min = 0.1;

        double t_max = numIterations;
        double t = iteration;

        /* New inertia weight using LDIW. */
        double w_fraction = (w_max - w_min) * (t_max - t) / t_max;
        double w = w_max - w_fraction;

        return w;
    }

    /**
     * Prints a 2-dimensional matrix.
     *
     * @param matrix Matrix to print.
     */
    private void printMatrix(ArrayList<int[]> matrix) {
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix.get(i).length; j++) {
                System.out.print(matrix.get(i)[j]);
            }

            System.out.println();
        }
    }

    /**
     * Calculates a random inertia weight.
     *
     * @param particleNumber   Particle number.
     * @param averageFitnesses Average particle fitnesses matrix.
     * @return The global inertia value w.
     */
    private double randomInertiaWeight(int particleNumber,
                               ArrayList<double[]> averageFitnesses) {
        int k = 5;
        double w;

        double p; // Annealing probability.

        double currentFitness =
                averageFitnesses.get(particleNumber)[iteration];
        double previousFitness =
                averageFitnesses.get(particleNumber)[iteration - k];

        if (previousFitness <= currentFitness) {
            p = 1;
        } else {
            /* Annealing temperature. */
            double coolingTemp_Tt;

            Particle currParticle = swarm.get(particleNumber);
            double bestFitness = currParticle.bestFitness;

            double particleFitnessAverage = 0;

            int counter = 0;
            for (int i = 0; i < iteration; i++) {
                if (averageFitnesses.get(particleNumber)[i] > 0) {
                    particleFitnessAverage +=
                            averageFitnesses.get(particleNumber)[i];
                    counter++;
                }
            }

            particleFitnessAverage = particleFitnessAverage / counter;
            coolingTemp_Tt = particleFitnessAverage / bestFitness - 1;
            p = Math.exp(-(previousFitness - currentFitness) / coolingTemp_Tt);
        }

        int random = rand.nextInt(2);

        /* New inertia weight. */
        if (p >= random) {
            w = random / 2.0 + 1;
        } else {
            w = random / 2.0;
        }

        return w;
    }

    /**
     * Rebalances the solution found by PSO for better solutions.
     *
     * @param newPositionsMatrix Positions matrix.
     * @param runTime            The run time matrix.
     * @return Rebalanced positions matrix.
     */
    private ArrayList<int[]> rebalancePSO(ArrayList<int[]> newPositionsMatrix,
                                          ArrayList<double[]> runTime) {
        boolean done = false;
        int counter = 0;

        while (!done) {
            double[] sum = new double[n];

            for (int i = 0; i < n; i++) {
                double[] time = runTime.get(i);
                int[] position = newPositionsMatrix.get(i);

                for (int j = 0; j < position.length; j++) {
                    if (position[j] == 1) {
                        sum[i] += time[j];
                    }
                }
            }

            int heaviestLoad = 0;
            int lightestLoad = 0;

            for (int i = 1; i < n; i++) {
                if (sum[heaviestLoad] < sum[i]) {
                    heaviestLoad = i;
                }
                if (sum[lightestLoad] < sum[i]) {
                    lightestLoad = i;
                }
            }

            int[] heaviestPosition = newPositionsMatrix.get(heaviestLoad);
            int[] lightestPosition = newPositionsMatrix.get(lightestLoad);

            for (int i = 0; i < heaviestPosition.length; i++) {
                int cloudletNumber = 0;

                if (heaviestPosition[i] == 1) {
                    cloudletNumber = 1;
                }

                double heaviestMinusCloudlet = sum[heaviestLoad] -
                        heaviestPosition[cloudletNumber];
                double lightestMinusCloudlet = sum[lightestLoad] -
                        lightestPosition[cloudletNumber];

                if (heaviestMinusCloudlet < lightestMinusCloudlet) {
                    break;
                } else {
                    heaviestPosition[cloudletNumber] = 0;
                    lightestPosition[cloudletNumber] = 1;
                    newPositionsMatrix.set(heaviestLoad, heaviestPosition);
                    newPositionsMatrix.set(lightestLoad, lightestPosition);
                }
            }

            if (counter == 3) {
                done = true;
            }

            counter++;
        }

        return newPositionsMatrix;
    }

    /**
     * Run the particle swarm optimization algorithm.
     *
     * @return vmIds List of Vm ids to match cloudlets with.
     */
    public ArrayList<Integer> run() {
        List<Integer> vmIds = new ArrayList<Integer>();
        /* Initialize swarm with numParticles particles. */
        initSwarm();

        for (int i = 0; i < numIterations; i++) {
            for (int j = 0; j < numParticles; j++) {
                /* Get current particle and attributes from swarm. */
                Particle currentParticle = swarm.get(j);
                // ArrayList<double[]> currentVelocity = currentParticle.velocity;
                ArrayList<int[]> currentPosition = currentParticle.position;

                /* Calculate inertia value. */
                double w = calcInertia(j, averageFitnesses);

                /* Calculate new velocities. */
                calcNewVelocities(w, currentParticle);

                /* Calculate new positions. */
                calcNewPositions(w, currentParticle);

                /* Calculate the fitness.
                 *
                 * @TODO: Modify calcFitness() to return void and take the
                 * particle as a parameter, updating the particle's fitness
                 * within the method itself.
                 */
                double newFitness = calcFitness(runTime, currentPosition);
                currentParticle.fitness = newFitness;

                /* Evaluate the solution and update the current particle and
                 * global bests if appropriate. */
                evaluateSolution(newFitness, currentPosition, currentParticle);

                /* Add the new fitness to the averageFitnesses array. */
                double[] fitnessArray = averageFitnesses.get(j);
                fitnessArray[i] = newFitness;
                averageFitnesses.set(j, fitnessArray);
            }

            iteration++;
        }
        return getCloudletPositions(globalBest);
    }

    /**
     * Sigmoid function.
     *
     * @param x Parameter of type double.
     */
    private double s(double x) {
        return 1 / (1 + Math.exp(-x));
    }

    /**
     * Set fixed inertia weight (FIW).
     *
     * @param fixedInertiaWeight The desired value.
     */
    private void setFixedInertiaWeight(double fixedInertiaWeight) {
        this.fixedInertiaWeight = fixedInertiaWeight;
    }
}
