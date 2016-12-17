package org.cloudbus.cloudsim.examples.data;

import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.core.CloudSimTags;
import org.cloudbus.cloudsim.core.SimEvent;
import org.cloudbus.cloudsim.lists.VmList;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by crosbylanham on 12/16/16.
 */
public class DataCenterBrokerData extends DatacenterBroker {
    protected double c1, c2, inertiaweight;

    public DataCenterBrokerData(double c1, double c2) throws Exception {
        super("Broker");
        this.c1 = c1;
        this.c2 = c2;
    }

    public DataCenterBrokerData(double c1, double c2, double inercia) throws Exception {
        super("Broker");
        this.c1 = c1;
        this.c2 = c2;
        this.inertiaweight = inercia;
    }

    /**
     * Submit cloudlets to the created VMs.
     *
     * @pre $none
     * @post $none
     * @see #submitCloudletList(java.util.List)
     */
    protected void submitCloudlets(double r1, double r2) {
        BinaryPSO pso = new BinaryPSOData(getVmsCreatedList(), getCloudletList(),
                1000, 100, 0, r1, r2);
        List<Integer> vmIds;
        vmIds = pso.run();
        // System.out.println("VM IDs:" + vmIds);
        int vmIndex = vmIds.get(0); // Get the first vmId in the vmIds list.
        List<Cloudlet> successfullySubmitted = new ArrayList<Cloudlet>();
        for (int i = 0; i < getCloudletList().size(); i++) {
            Cloudlet cloudlet = getCloudletList().get(i);
            Vm vm;
            // if user didn't bind this cloudlet and it has not been executed yet
            if (cloudlet.getVmId() == -1) {
                // TODO: Fix this.. it always submits to the first Vm.
                vm = vmsCreatedList.get(vmIds.get(i));
            } else { // submit to the specific vm
                vm = VmList.getById(getVmsCreatedList(), cloudlet.getVmId());
                if (vm == null) { // vm was not created
                    if (!Log.isDisabled()) {
                        Log.printConcatLine(CloudSim.clock(), ": ", getName(),
                                ": Postponing execution of cloudlet ",
                                cloudlet.getCloudletId(), ": bount VM not available");
                    }
                    continue;
                }
            }

            if (!Log.isDisabled()) {
                Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Sending cloudlet ",
                        cloudlet.getCloudletId(), " to VM #", vm.getId());
            }

            cloudlet.setVmId(vm.getId());
            sendNow(getVmsToDatacentersMap().get(vm.getId()), CloudSimTags.CLOUDLET_SUBMIT, cloudlet);
            cloudletsSubmitted++;
            // vmIndex = (vmIndex + 1) % getVmsCreatedList().size();
            getCloudletSubmittedList().add(cloudlet);
            successfullySubmitted.add(cloudlet);
        }

        // remove submitted cloudlets from waiting list
        getCloudletList().removeAll(successfullySubmitted);
    }

    protected void submitCloudlets(double r1, double r2, double inercia) {
        BinaryPSO pso = new BinaryPSOData(getVmsCreatedList(), getCloudletList(),
                1000, 100, 0, r1, r2, inercia);
        List<Integer> vmIds;
        vmIds = pso.run();
        // System.out.println("VM IDs:" + vmIds);
        int vmIndex = vmIds.get(0); // Get the first vmId in the vmIds list.
        List<Cloudlet> successfullySubmitted = new ArrayList<Cloudlet>();
        for (int i = 0; i < getCloudletList().size(); i++) {
            Cloudlet cloudlet = getCloudletList().get(i);
            Vm vm;
            // if user didn't bind this cloudlet and it has not been executed yet
            if (cloudlet.getVmId() == -1) {
                // TODO: Fix this.. it always submits to the first Vm.
                vm = vmsCreatedList.get(vmIds.get(i));
            } else { // submit to the specific vm
                vm = VmList.getById(getVmsCreatedList(), cloudlet.getVmId());
                if (vm == null) { // vm was not created
                    if (!Log.isDisabled()) {
                        Log.printConcatLine(CloudSim.clock(), ": ", getName(),
                                ": Postponing execution of cloudlet ",
                                cloudlet.getCloudletId(), ": bount VM not available");
                    }
                    continue;
                }
            }

            if (!Log.isDisabled()) {
                Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Sending cloudlet ",
                        cloudlet.getCloudletId(), " to VM #", vm.getId());
            }

            cloudlet.setVmId(vm.getId());
            sendNow(getVmsToDatacentersMap().get(vm.getId()), CloudSimTags.CLOUDLET_SUBMIT, cloudlet);
            cloudletsSubmitted++;
            // vmIndex = (vmIndex + 1) % getVmsCreatedList().size();
            getCloudletSubmittedList().add(cloudlet);
            successfullySubmitted.add(cloudlet);
        }

        // remove submitted cloudlets from waiting list
        getCloudletList().removeAll(successfullySubmitted);
    }

    protected void processVmCreate(SimEvent ev) {
        int[] data = (int[]) ev.getData();
        int datacenterId = data[0];
        int vmId = data[1];
        int result = data[2];

        if (result == CloudSimTags.TRUE) {
            getVmsToDatacentersMap().put(vmId, datacenterId);
            getVmsCreatedList().add(VmList.getById(getVmList(), vmId));
            Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": VM #", vmId,
                    " has been created in Datacenter #", datacenterId, ", Host #",
                    VmList.getById(getVmsCreatedList(), vmId).getHost().getId());
        } else {
            Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Creation of VM #", vmId,
                    " failed in Datacenter #", datacenterId);
        }

        incrementVmsAcks();

        // all the requested VMs have been created
        if (getVmsCreatedList().size() == getVmList().size() - getVmsDestroyed()) {
            submitCloudlets(c1, c2);
        } else {
            // all the acks received, but some VMs were not created
            if (getVmsRequested() == getVmsAcks()) {
                // find id of the next datacenter that has not been tried
                for (int nextDatacenterId : getDatacenterIdsList()) {
                    if (!getDatacenterRequestedIdsList().contains(nextDatacenterId)) {
                        createVmsInDatacenter(nextDatacenterId);
                        return;
                    }
                }

                // all datacenters already queried
                if (getVmsCreatedList().size() > 0) { // if some vm were created
                    submitCloudlets(c1, c2);
                } else { // no vms created. abort
                    Log.printLine(CloudSim.clock() + ": " + getName()
                            + ": none of the required VMs could be created. Aborting");
                    finishExecution();
                }
            }
        }
    }
}
