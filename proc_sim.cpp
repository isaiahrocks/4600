/* Yosvany Blanco
 * Miguel Martinez
 * CSCE 4600-003
 * project 1
 *
 * COMPILE WITH: g++ -std=c++0x proc_sim.cpp -o project
 * Simulation of different scheduling methods for operating systems
*/


#include <iostream>
#include <vector>	// I love vectors so much
#include <random>

int const QUANTUM = 50; // Quantum cycle
int const CONTEXT = 10; // context switch
int const NUM_OF_PROCS = 50;
// process class
class Proc
{
public:
	int pid; // pid >= 0
	int cycles; // 1000 -> 11,000
	int run_time;
	int mem_size; // in KB 1 -> 100
	int start_time;
	int arrival_time;
	int finish_time;
	int wait_time;
	bool done;
	
	Proc(int mypid, int mycycles, int mymem_size, int myarrival_time)
	{ // simple constructor used to make each process
		this->pid = mypid;
		this->run_time = this->cycles = mycycles;
		this->mem_size = mymem_size;
		this->arrival_time = myarrival_time;
		this->start_time = 0;
		this->finish_time = 0;
		this->wait_time = 0;
		this->done = false;
	}

	void print()
	{
		std::cout << "pid: " << this->pid << " , cycles: " << this->cycles;
		std::cout << " , memory footprint: " << mem_size << ", arrival time: " << arrival_time
		<< std::endl;
	}

};

class Processor
{
public:
	int avg_wait, context_penalty;
	std::vector<Proc> procs; // each processor will have it's collection of processes
	Processor(std::vector<Proc> myprocs, int start, int end)
	{
		this->avg_wait = 0;
		this->context_penalty = 0;
		for (int i = start; i < end; i++ )
		{ //get the sub vector
			procs.push_back(myprocs[i]);
		}
	}
};

class quantum_run
{
public:
	int pid;
	int cycle_run;
	int start_time;
	quantum_run(int mypid, int mycycle_run, int mystart_time)
	{
		this->pid = mypid;
		this->cycle_run = mycycle_run;
		this->start_time = mystart_time;
	}
};

//////////////// random number generator that generates within a provided range
//////////////// I found this method online at the following link 3rd answer down
//////////////// http://stackoverflow.com/questions/28531724/how-do-i-generate-random-numbers-with-a-mean-in-c
//////////////// I don't claim to have come up with this specific function

double randAvg(double minimum, double mean, double maximum) 
{

    // random generator
	static std::random_device rd;
    static std::mt19937 generator(rd());

    //This forces the random number to fit the given range 
    double const average_bound_width = ((mean-minimum) + (maximum-mean)) / 2;
    double const standard_deviation  = average_bound_width / 3;
    std::normal_distribution<double> distribution(mean, standard_deviation);

    // gets rid of bad values
    double value;
    do {
        value = distribution(generator);
    } while (value < minimum || maximum < value);

    return value;
}

bool all_done(std::vector<Proc> &set_of_procs)
{
	for (auto it = set_of_procs.begin(); it != set_of_procs.end(); it++)
	{
		if(it->done == false)
		{
			return false; // it's true that the processes are all done
		}
	}
	return true; // atleast one process is not done
}

void create_procs(std::vector<Proc> &set_of_procs)
{
	// creat number of processes/////////////////////////////////////////
	int num_of_proc;
	// specify number of processes
	num_of_proc = NUM_OF_PROCS;
	// here I'll just generate a random number within each of the 
	// given ranges and distribute evenly
	int cycle_arrival = 0;
	
	for (int pid = 0; pid < num_of_proc; pid++ )
	{
		int rand_cycles = (int)randAvg(1000,6000,11000);
		// 1000,6000,11000 // 50,70,100
		int rand_mem_size = (int)randAvg(1,50,100);
		Proc temp(pid, rand_cycles, rand_mem_size, cycle_arrival);
		set_of_procs.push_back(temp);
		cycle_arrival += 50;
	}


	for (int i = 0; i < num_of_proc; i++)
	{
		set_of_procs[i].print();
	}
	//////////////////////////////////////////////////////////////////////////
	
}

// this function also happens to make a gant chart that we only used for debugging purposes
void RR(std::vector<Proc> set_of_procs, int proc_num)
{// because the instructions have us assume everything comes in every 50 cycles and the quantum also happens to be 50 cycles
	// this oversimplifies the solution for RR
	int i = 0;
	int timing = 0;
	int context_penalty = 0;
	bool check = all_done(set_of_procs);
	std::vector<quantum_run> round_robin; // this is a gant chart made for debugging purposes only
	while ( !check ) // while the processes are not done running
	{
		if (set_of_procs[i].done == true)
		{
			if (i == set_of_procs.size() - 1)
				i = 0;
			else
				i++;
			continue; // skip to next
		}
		int pid;
		int cycle_run;
		int start_time;
		
		if (set_of_procs[i].cycles - QUANTUM < 0)
		{
			cycle_run = set_of_procs[i].cycles;
			pid = set_of_procs[i].pid;
			set_of_procs[i].cycles -= set_of_procs[i].cycles; // adjust time left
			start_time = timing; // start time of this burst is whatever current cycle is
			set_of_procs[i].done = true;
			timing += cycle_run ; //update current time and account for context switch
			set_of_procs[i].finish_time = timing;
			timing += CONTEXT;
			context_penalty += CONTEXT;
		}
		else
		{
			pid = set_of_procs[i].pid; // save pid info
			cycle_run = QUANTUM; // burst time is the quantum time
			start_time = timing; // start time of this burst is whatever current cycle is
			set_of_procs[i].cycles -= QUANTUM; // adjust time left
			timing += (cycle_run + CONTEXT) ; //update current time and account for context switch
			context_penalty += CONTEXT;
		}
		
		quantum_run temp(pid, cycle_run, start_time);
		
		/////////////////////////////////
		if (i == set_of_procs.size() - 1)
			i = 0;
		else
			i++;
		
		check = all_done(set_of_procs);
		round_robin.push_back(temp); // update round robin vector
		//////////////////////////////////
	}

	// GANTT chart
	/*std::cout << "GANTT chart below. Context switch of " << CONTEXT << " cycles should be accounted for below!" << std::endl;
	for (auto it = round_robin.begin(); it!=round_robin.end(); it++)
	{
		std::cout << "pid: " << it->pid << " burst: " << it->cycle_run << " start time: " << it->start_time << std::endl;
	} // debugging
	std::cout << "GANTT chart above" << std::endl;*/
	// GANTT chart ^^^^^^^^
	float average = 0;
	for (auto it = set_of_procs.begin(); it!=set_of_procs.end(); it++)
	{
		average += (it->finish_time - it->arrival_time - it->run_time);
	}
	average = (float)average / proc_num;
	std::cout << "Round Robin average wait time: "<<average << std::endl;
	std::cout << "Round Robin context penalty: "<<context_penalty << std::endl;
}

void SJF(std::vector<Proc> set_of_procs, int proc_num)
{
	std::vector<quantum_run> shortest_job_first; // this is a gant chart made for debugging purposes only
	int i = 0;
	int timing = 0;
	int context_penalty = 0;
	int pid = set_of_procs[0].pid;
	int cycle_run = set_of_procs[0].cycles;
	int start_time = timing;
	quantum_run temp(pid,cycle_run,start_time);
	shortest_job_first.push_back(temp);
	timing += cycle_run;
	set_of_procs[0].finish_time = timing;
	timing += CONTEXT; // account for context switch
	context_penalty += CONTEXT;
	
	// sort the jobs
	for(int i = 1; i < proc_num; i++)
	{
		for(int j = 1; j < proc_num - i - 1; ++j)
		{
			if (set_of_procs[j].cycles > set_of_procs[j+1].cycles)
			{
				std::swap(set_of_procs[j], set_of_procs[j + 1]);
			}
		}
	}
	//////////////////////////////////////////////////////////////
	
	for(int i = 1; i < proc_num; i++)
	{
		pid = set_of_procs[i].pid;
		cycle_run = set_of_procs[i].cycles;
		start_time = timing;
		quantum_run temp1(pid,cycle_run,start_time);
		shortest_job_first.push_back(temp1);
		timing += cycle_run;
		set_of_procs[i].finish_time = timing;
		timing += CONTEXT;
		context_penalty += CONTEXT;
	}
	context_penalty -= 10; // get rid of extra context switch
	// GANTT chart
	/*std::cout << "GANTT chart below. Context switch of " << CONTEXT << " cycles should be accounted for below!" << std::endl;
	for (auto it = shortest_job_first.begin(); it!=shortest_job_first.end(); it++)
	{
		std::cout << "pid: " << it->pid << " burst: " << it->cycle_run << " start time: " << it->start_time << std::endl;
	} // debugging
	std::cout << "GANTT chart above" << std::endl;*/
	// GANTT chart ^^^^^^^^
	float average = 0;
	for (auto it = set_of_procs.begin(); it!=set_of_procs.end(); it++)
	{
		average += (it->finish_time - it->arrival_time - it->run_time);
	}
	average = (float)average / proc_num;
	std::cout << "Shortest job first average wait time: "<<average << std::endl;
	std::cout << "Shortest job first context penalty: "<<context_penalty << std::endl;
}

void FIFO(std::vector<Proc> set_of_procs, int proc_num)
{
	int timing = 0, t2 = 0; //keep track of finish time and start time of each process
	int context_penalty = 0;
	for(int i = 0; i < proc_num; i++)
	{
		if(i == 0)
		{ //process 0 doesn't have a context switch to start
			t2 = timing;
			set_of_procs[i].start_time = timing;
		}
		else
		{ //add context
			t2 = (timing + CONTEXT);
			set_of_procs[i].start_time = (timing+CONTEXT);
			context_penalty += CONTEXT;
			
		}
		timing += set_of_procs[i].cycles; //adds the number of cycles
	//	set_of_procs[i].start_time = t2+CONTEXT;
		set_of_procs[i].finish_time = timing+CONTEXT; //finishing point of process
		// GANTT chart
		/*
		if(i == 0)
			std::cout << "Process " << i << " start: " << t2 << " finish: " << timing << std::endl;
		else
			std::cout << "Process " << i << " start: " << t2 << " finish: " << (timing+CONTEXT) << std::endl;
		*/
		////////////////////////////////////////
	}
  float average = 0;
	for(auto it = set_of_procs.begin(); it != set_of_procs.end(); it++) //calculate the average waiting time
	{
		average += (it->start_time - it->arrival_time);
	}
	average = (float)average/proc_num;
	std::cout << "FIFO average wait time: " << average << std::endl;
	std::cout << "FIFO context penalty: " << context_penalty << std::endl;
}

void run_suite(std::vector<Proc> set_of_procs, int proc_num)
{
	RR(set_of_procs,proc_num);
        SJF(set_of_procs, proc_num);
        FIFO(set_of_procs, proc_num);
}

///////////////
int main()
{
	std::vector<Proc> set_of_procs;
	create_procs(set_of_procs);
	//run the scheduling algorithms on all 50 processes
	RR(set_of_procs,NUM_OF_PROCS );
	SJF(set_of_procs, NUM_OF_PROCS);
	FIFO(set_of_procs, NUM_OF_PROCS);
	///////////////////////////////////////////////////
	// set up the 4 processors to distribute the processes
	std::vector<Processor> processors;// create processor array
	processors.push_back(Processor(set_of_procs, 0,12) ); //13, 13, 13, 11 these are the number of processes that will be used per processor
	processors.push_back(Processor(set_of_procs, 13,25) );
	processors.push_back(Processor(set_of_procs, 26,38) );
	processors.push_back(Processor(set_of_procs, 39,49) );
	//////////////////////////////////////////////////////

	
	// simulate the 4 processors running
	for (int i = 0; i < 4; i++)
	{
		std::cout << "Processor " << i << "#############################" << std::endl;
		run_suite(processors[i].procs,processors[i].procs.size());
		std::cout << "#########################################################" << std::endl;
	}


return 0;
}
