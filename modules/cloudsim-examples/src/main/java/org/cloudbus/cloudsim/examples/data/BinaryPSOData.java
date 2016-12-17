/*
 * Title:        CloudSim Toolkit
 * Description:  CloudSim (Cloud Simulation) Toolkit for Modeling and 
 * Simulation of Clouds
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2009-2012, The University of Melbourne, Australia
 */

package org.cloudbus.cloudsim.examples.data;

import org.cloudbus.cloudsim.BinaryPSO;
import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.Vm;

import java.util.List;

/**
 * BinaryPSO implements the popular population-based metaheuristic algorithm
 * for scheduling cloudlets on VMs.
 *
 * @author Amlan Chatterjee, Crosby Lanham, Matt Levan, Mishael Zerrudo
 */
public class BinaryPSOData extends BinaryPSO {

    public BinaryPSOData(List<Vm> vmList, List<Cloudlet> cloudletList,
                         int numIterations, int numParticles, int inertiaTechnique,
                         double c1, double c2) {
        super(vmList, cloudletList, numIterations, numParticles, inertiaTechnique);
        setC1(c1);
        setC2(c2);

    }

    public BinaryPSOData(List<Vm> vmList, List<Cloudlet> cloudletList,
                         int numIterations, int numParticles, int inertiaTechnique,
                         double c1, double c2, double fixedinertiaweight) {
        super(vmList, cloudletList, numIterations, numParticles, inertiaTechnique);
        setC1(c1);
        setC2(c2);
        fixedInertiaWeight = fixedinertiaweight;
    }

    public void setC1(double c) {
        c1 = c;
    }

    public void setC2(double c) {
        c2 = c;
    }
}
