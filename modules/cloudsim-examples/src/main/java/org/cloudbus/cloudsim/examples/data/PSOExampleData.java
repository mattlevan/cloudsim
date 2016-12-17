package org.cloudbus.cloudsim.examples.data;

import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.examples.PSOExample;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

/**
 * Created by crosbylanham on 12/16/16.
 */

public class PSOExampleData extends PSOExample {
    private static final String FILENAME = "Data.txt";

    public static void main(String[] args) {
        Log.printLine("Starting PSOExample...");
        // First step: Initialize the CloudSim package. It should be called
        // before creating any entities.
/*----------------------------------------------------------------------------------------- */
        //change iterations you can change what it tests
        for (int c1 = 0; c1 <= 300; c1 += 10) {
            for (int c2 = 0; c2 <= 300; c2 += 10) {
                writeToFile(String.format("r1: %7.3f  r2: %7.3f  ", c1 / 100.0, c2 / 100.0));
                for (int tests = 0; tests < 10; tests++) {
                    try {
                        int num_user = 1;   // number of cloud users
                        Calendar calendar = Calendar.getInstance();
                        boolean trace_flag = false;  // mean trace events

                        // Initialize the CloudSim library
                        CloudSim.init(num_user, calendar, trace_flag);

                        // Second step: Create Datacenters
                        // Datacenters are the resource providers in CloudSim.
                        // We need at least one of them to run a CloudSim simulation
                        @SuppressWarnings("unused")
                        Datacenter datacenter0 = createDatacenter("Datacenter_0");


                        //Third step: Create Broker
                        DatacenterBroker broker = createBroker(c1, c2);
                        int brokerId = broker.getId();

                        //Fourth step: Create one virtual machine
                        vmlist = new ArrayList<Vm>();

                        //VM description
                        int vmid = 0;
                        int mips = 128;
                        long size = 10000; //image size (MB)
                        int ram = 512; //vm memory (MB)
                        long bw = 1000;
                        int pesNumber = 1; //number of cpus
                        int numVms = 10; //number of vms
                        String vmm = "Xen"; //VMM name

                        for (int i = 0; i < numVms; i++) {
                            Vm vm = new Vm(i, brokerId, mips, pesNumber, ram, bw,
                                    size, vmm, new CloudletSchedulerTimeShared());

                            vmlist.add(vm);
                        }

                        //submit vm list to the broker
                        broker.submitVmList(vmlist);


                        //Fifth step: Create two Cloudlets
                        cloudletList = new ArrayList<Cloudlet>();

                        //Cloudlet properties
                        long length = 40000;
                        long fileSize = 300;
                        long outputSize = 300;
                        int numJobs = 100;
                        UtilizationModel utilizationModel = new UtilizationModelFull();

                        for (int i = 0; i < numJobs; i++) {
                            Cloudlet cloudlet = new Cloudlet(i, length, pesNumber,
                                    fileSize, outputSize, utilizationModel,
                                    utilizationModel, utilizationModel);

                            cloudlet.setUserId(brokerId);
                            cloudletList.add(cloudlet);
                        }


                        //submit cloudlet list to the broker
                        broker.submitCloudletList(cloudletList);

                        // Sixth step: Starts the simulation
                        CloudSim.startSimulation();


                        // Final step: Print results when simulation is over
                        List<Cloudlet> newList = broker.getCloudletReceivedList();

                        CloudSim.stopSimulation();
                        printCloudletList(newList);

                        Log.printLine("PSOExample finished!");
                    } catch (Exception e) {
                        e.printStackTrace();
                        Log.printLine("The simulation has been terminated due to an unexpected error");
                    }
                }
                writeToFile("\n");
            }
        }
    }

    public static void printCloudletList(List<Cloudlet> list) {
        int size = list.size();
        Cloudlet cloudlet;

        String indent = "    ";
        Log.printLine();
        Log.printLine("========== OUTPUT ==========");
        Log.printLine("Cloudlet ID" + indent + "STATUS" + indent +
                "Data center ID" + indent + "VM ID" + indent + "Time" + indent + "Start Time" + indent + "Finish Time");

        DecimalFormat dft = new DecimalFormat("###.##");
        for (int i = 0; i < size; i++) {
            cloudlet = list.get(i);
            Log.print(indent + cloudlet.getCloudletId() + indent + indent);

            if (cloudlet.getCloudletStatus() == Cloudlet.SUCCESS) {
                Log.print("SUCCESS");

                Log.printLine(indent + indent + cloudlet.getResourceId() + indent + indent + indent + cloudlet.getVmId() +
                        indent + indent + dft.format(cloudlet.getActualCPUTime()) + indent + indent + dft.format(cloudlet.getExecStartTime()) +
                        indent + indent + dft.format(cloudlet.getFinishTime()));


            }
        }
        String pl = String.format(" %10.3f ", list.get(list.size() - 1).getFinishTime());
        writeToFile(pl);
    }

    protected static void writeToFile(String content) {
        BufferedWriter bw = null;
        FileWriter fw = null;
        try {
            fw = new FileWriter(FILENAME, true);
            bw = new BufferedWriter(fw);
            bw.write(content);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (bw != null)
                    bw.close();
                if (fw != null)
                    fw.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    protected static DatacenterBroker createBroker(double r1, double r2) {
        DatacenterBroker broker = null;
        try {
            broker = new DataCenterBrokerData(r1, r2);
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        return broker;
    }

    protected static DatacenterBroker createBroker(double r1, double r2, double inercia) {

        DatacenterBroker broker = null;
        try {
            broker = new DataCenterBrokerData(r1, r2, inercia);
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        return broker;
    }
}
