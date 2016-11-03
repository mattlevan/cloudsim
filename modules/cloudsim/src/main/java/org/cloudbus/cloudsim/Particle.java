/**
 * Particle class implements the particle for building a swarm.
 *
 * @todo: How exactly is the fitness calculated? Does each particle calculate 
 * execution times of the entire power set of cloudlets on only a single Vm or 
 * on every Vm? 
 *
 * @todo: Also, does each particle only calculate execution times of a 
 * partition of the power set of cloudlets?
 *
 * @todo: How does the particle express the personal current and best 
 * solutions (in terms of pairings)?
 */
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

    /** Calculate new velocities. */
    public void calcNewVelocities(double w, double c1, double c2, double g) {
        for (int i = 0; i < v.length; i++) {
            r1 = random.nextInt();
            r2 = random.nextInt();
            v[i] = w*v[i]+c1*r1*(best-current)+c2*r2*(g-current);
        }
    }

    /** Calculate personal current solution. */
    public void calcCurrent() {
        p[current] = current+v;
    }

    /** Sigmoid function. */
    public double s(double current) {
        return 1/(1+Math.exp(-1*current));
    }

    /** Calculate position. */
    public double calcPosition() {
        if (s(p[current]) <= r)
            p[current] = 0;
        else
            p[current]  = 1;
    }

}
