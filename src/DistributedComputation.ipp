//////////////////////////////////////////////////////////////////////
// DistributedComputation.cpp
//
// Class for performing distributed optimization.
//////////////////////////////////////////////////////////////////////

#include "DistributedComputation.hpp"

enum CommandType
{ 
    CommandType_LoadSharedData, 
    CommandType_DoWork, 
    CommandType_SendResultSize,
    CommandType_SendResult, 
    CommandType_ClearResult,
    CommandType_Quit
};

//////////////////////////////////////////////////////////////////////
// DistributedComputationBase::DistributedComputationBase()
//
// Constructor.  Performs MPI initializations if MULTI is
// defined.
//////////////////////////////////////////////////////////////////////

template<class SharedData, class NonSharedData>
DistributedComputationBase<SharedData, NonSharedData>::DistributedComputationBase(bool toggle_verbose) : 
    toggle_verbose(toggle_verbose),
    processing_time(0),
    total_time(0),
    id(0),
    num_procs(1)
{

#ifdef MULTI
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if (id == 0 && toggle_verbose)
    {
        WriteProgressMessage("");
        std::cerr << "Distributed Optimization Library started.  Using " 
                  << num_procs << " processor(s)." << std::endl;
    }
}


//////////////////////////////////////////////////////////////////////
// DistributedComputationBase::RunAsComputeNode()
//
// Turn into a compute node and process work requests from the master
// node until the command to quit is sent.  Should only be called
// ifdef MULTI is defined.
//////////////////////////////////////////////////////////////////////

template<class SharedData, class NonSharedData>
void DistributedComputationBase<SharedData, NonSharedData>::RunAsComputeNode()
{
    Assert(id != 0, "Routine should not be called by master process.");
    
#ifdef MULTI
    MPI_Status status;
    SharedData shared_data;
    NonSharedData nonshared_data;
    std::vector<RealT> result;
    std::vector<RealT> partial_result;
    
    while (true)
    {
        // block until command received
        int command;
        MPI_Recv(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        
        switch (command)
        {
            case CommandType_LoadSharedData:
            {
                // get shared data
                MPI_Bcast(&shared_data, sizeof(SharedData), MPI_BYTE, 0, MPI_COMM_WORLD);
            }
            break;
            
            case CommandType_DoWork:
            {
                // get nonshared data
                MPI_Recv(&nonshared_data, sizeof(NonSharedData), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);

                // perform and time computation
                processing_time = GetSystemTime();
                DoComputation(partial_result, shared_data, nonshared_data);
                processing_time = GetSystemTime() - processing_time;

                // resize results vector as needed
                if (result.size() == 0)
                    result.resize(partial_result.size());
                else if (result.size() != partial_result.size())
                    Error("Encountered return values of different size.");
                
                // accumulate results
                result += partial_result;
                
                // return processing time to main node
                MPI_Send(&processing_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            break;
            
            case CommandType_SendResultSize:
            {
                // send result size to main node
                int size = result.size();
                MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            break;
            
            case CommandType_SendResult:
            {
                // make sure all results are of the same length first
                int size;
                MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

                if (size == 0)
                    Error("Should not call this unless results expected.");
                else if (result.size() == 0)
                    result.resize(size);
                else if (size != int(result.size()))
                    Error("Return values of different size (%d vs %d).", size, int(result.size()));

                // then send result to main node
                MPI_Reduce(&result[0], NULL, result.size(), GetResultMPIDataType(), MPI_SUM, 0, MPI_COMM_WORLD);
            }
            break;
            
            case CommandType_ClearResult:
                result.clear();
                break;
                
            case CommandType_Quit:
                return;
        }      
    }
#endif

}

//////////////////////////////////////////////////////////////////////
// DistributedComputationBase::StopComputeNodes()
//
// Closes down MPI connections for compute nodes.
//////////////////////////////////////////////////////////////////////

template<class SharedData, class NonSharedData>
void DistributedComputationBase<SharedData, NonSharedData>::StopComputeNodes()
{
#ifdef MULTI
    if (id == 0)
    {
        for (int i = 1; i < num_procs; i++)
        {
            int command = CommandType_Quit;
            MPI_Send(&command, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
#endif    
}

//////////////////////////////////////////////////////////////////////
// DistributedComputationBase::DistributeComputation()
//
// Distribute computation tasks among all nodes (other than 0) if
// MULTI is defined; work units are allocated starting from largest
// unit size to smallest unit size.
//////////////////////////////////////////////////////////////////////

const int NOT_ALLOCATED = -1;
const int DO_NOT_ALLOCATE = -2;

template<class SharedData, class NonSharedData>
void DistributedComputationBase<SharedData, NonSharedData>::DistributeComputation(std::vector<RealT> &result,
                                                                                         const SharedData &shared_data,
                                                                                         const std::vector<NonSharedData> &nonshared_data)
{
    Assert(id == 0, "Routine should only be called by master process.");
    Assert(nonshared_data.size() > 0, "Must submit at least one work description for processing.");

    double starting_time = GetSystemTime();
    size_t units_complete = 0;

    result.clear();
    
#ifdef MULTI
    size_t num_procs_in_use = 1;
    size_t curr_unit = 0;
    int command;
    int size;

    MPI_Status status;
    std::string progress;

    // initialize work assignments
    std::vector<int> assignment(num_procs, NOT_ALLOCATED);
    assignment[0] = DO_NOT_ALLOCATE;

    // clear accumulated result on all processors
    if (toggle_verbose) WriteProgressMessage("Clearing accumulated result on all processors.");
    command = CommandType_ClearResult;
    for (int proc = 1; proc < num_procs; proc++)
        MPI_Send(&command, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
    
    // broadcast shared data to all processors
    if (toggle_verbose) WriteProgressMessage("Broadcasting shared data to all processors.");
    command = CommandType_LoadSharedData;
    for (int proc = 1; proc < num_procs; proc++)
        MPI_Send(&command, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
    MPI_Bcast(const_cast<SharedData *>(&shared_data), sizeof(SharedData), MPI_BYTE, 0, MPI_COMM_WORLD);

    // while there is work to be done
    if (toggle_verbose) WriteProgressMessage("Sending work units to all processors.");
    while (num_procs_in_use > 1 || curr_unit < nonshared_data.size())
    {
        // allocate the max number of processors possible
        while (int(num_procs_in_use) < num_procs && curr_unit < nonshared_data.size())
        {
            // find free processor
            size_t proc = 0;
            while (proc < assignment.size() && assignment[proc] != NOT_ALLOCATED) proc++;
            Assert(proc < assignment.size(), "Expected to find free processor.");
            
            // send command            
            command = CommandType_DoWork;
            MPI_Send(&command, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
            
            // send nonshared data
            MPI_Send(const_cast<NonSharedData *>(&nonshared_data[curr_unit]), sizeof(NonSharedData), MPI_BYTE, proc, 0, MPI_COMM_WORLD);

            // update processor allocation table
            num_procs_in_use++;
            assignment[proc] = curr_unit;
            curr_unit++;
        }
        
        // write progress message (at most 1 update per second)
        double current_time = GetSystemTime();
        static double prev_reporting_time = 0;
        if (current_time - prev_reporting_time > 1)
        {
            prev_reporting_time = current_time;
            size_t percent_complete = 100 * units_complete / nonshared_data.size();
            if (toggle_verbose) WriteProgressMessage(SPrintF("%u/%u work units allocated, %d%% complete.", curr_unit, nonshared_data.size(), percent_complete));
        }
        
        // if no processors left, or all work allocated, then wait for results        
        MPI_Recv(&current_time, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        Assert(current_time >= 0, "Expected positive time value for acknowledgment of job completion.");
        processing_time += current_time;

        // update processor allocation table
        num_procs_in_use--;
        assignment[status.MPI_SOURCE] = -1;
        units_complete++;
    }
    
    // get accumulated result size
    if (toggle_verbose) WriteProgressMessage("Computing result size.");
    command = CommandType_SendResultSize;
    MPI_Send(&command, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Recv(&size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
    Assert(size >= 0, "Size should be nonnegative.");
    result.resize(size);
    
    // check if result is expected
    if (size > 0)
    {
        // tell all processors to send results
        if (toggle_verbose) WriteProgressMessage("Requesting results from processors.");
        command = CommandType_SendResult;
        for (int proc = 1; proc < num_procs; proc++)
            MPI_Send(&command, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);

        // check that results are all of the appropriate size
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // retrieving accumulated results
        if (toggle_verbose) WriteProgressMessage("Receiving accumulated results from processors.");
        MPI_Reduce(MPI_IN_PLACE, &result[0], size, GetResultMPIDataType(), MPI_SUM, 0, MPI_COMM_WORLD);
    }
   

#else
    
    // retrieve one result at a time, and accumulate
    std::vector<RealT> partial_result;    
    if (toggle_verbose) WriteProgressMessage("Starting first work unit.");
    for (size_t j = 0; j < nonshared_data.size(); j++)
    {
        DoComputation(partial_result, shared_data, nonshared_data[j]);

        // resize results vector as needed
        if (result.size() == 0)
            result.resize(partial_result.size());
        else if (result.size() != partial_result.size())
            Error("Encountered return values of different size.");
        
        // accumulate results
        result += partial_result;
        units_complete++;
        
        // write progress message (at most 1 update per second)
        double current_time = GetSystemTime();
        static double prev_reporting_time = 0;
        if (current_time - prev_reporting_time > 1)
        {
            prev_reporting_time = current_time;
            size_t percent_complete = 100 * units_complete / nonshared_data.size();
            if (toggle_verbose) WriteProgressMessage(SPrintF("%u/%u work units allocated, %d%% complete.", units_complete, nonshared_data.size(), percent_complete));
        }
    }
    
#endif
    
    if (toggle_verbose) WriteProgressMessage("");
    total_time += (GetSystemTime() - starting_time);
}

//////////////////////////////////////////////////////////////////////
// DistributedComputationBase::GetEfficiency()
//
// Compute the processor usage efficiency.
//////////////////////////////////////////////////////////////////////

template<class SharedData, class NonSharedData>
double DistributedComputationBase<SharedData, NonSharedData>::GetEfficiency() const
{
    Assert(IsMasterNode(), "Should only be called by master node.");
#ifdef MULTI
    return 100.0 * (processing_time / (num_procs - 1)) / (1e-10 + total_time);
#else
    return 100.0;
#endif
}

//////////////////////////////////////////////////////////////////////
// DistributedComputationBase::ResetEfficiency()
//
// Reset efficiency statistics.
//////////////////////////////////////////////////////////////////////

template<class SharedData, class NonSharedData>
void DistributedComputationBase<SharedData, NonSharedData>::ResetEfficiency()
{
    Assert(IsMasterNode(), "Should only be called by master node.");
    processing_time = total_time = 0;
}

/*

//////////////////////////////////////////////////////////////////////
// DistributedComputation::DistributedComputation()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT, class SharedData, class NonSharedData>
DistributedComputation<RealT, SharedData, NonSharedData>::DistributedComputation(bool toggle_verbose) :
    DistributedComputationBase(toggle_verbose)
{}

//////////////////////////////////////////////////////////////////////
// DistributedComputation::~DistributedComputation()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT, class SharedData, class NonSharedData>
DistributedComputation<RealT, SharedData, NonSharedData>::~DistributedComputation()
{}

//////////////////////////////////////////////////////////////////////
// DistributedComputation::GetResultMPIDataType()
//
// Routine for indicating data type to be transmitted.
//////////////////////////////////////////////////////////////////////

template<class RealT, class SharedData, class NonSharedData>
int DistributedComputation<RealT, SharedData, NonSharedData>::GetResultMPIDataType() const
{
    return 0;
}

#ifdef MULTI

template<class SharedData, class NonSharedData>
class DistributedComputation<float, SharedData, NonSharedData> : public DistributedComputationBase<float, SharedData, NonSharedData>
{
    int GetResultMPIDataType() const { return MPI_FLOAT; }
};

template<class SharedData, class NonSharedData>
class DistributedComputation<double> : public DistributedComputationBase<double, SharedData, NonSharedData>
{
    int GetResultMPIDataType() const { return MPI_DOUBLE; }
};

#endif
*/
