#include "inputoutput.h"

class dolfinVersionContents//dummy class
{
	friend std::ostream &operator<<(std::ostream &os, const dolfinVersionContents &io);
};

	int inputoutput::parseInputLine(int argc, char* argv[])
	{
		std::cout<<"argc= "<<argc<<std::endl;
		if(argc == 1) isLocalComputer = true;
		if(argc > 1)
		{
			std::cout<<"argv[] = ";
			for (int i=0; i<argc; i++)
				std::cout<<argv[i]<<" ";
			std::cout<<std::endl;

		if(std::string(argv[1]) == std::string("cluster.math.uni-augsberg.de"))
		{
		  isAugsbergCluster = true;   
		}
		else if(std::string(argv[1]) == std::string("opuntia.cacds.uh.edu"))
		{
		  isOpuntiaCluster = true;   
		}
			else if(std::string(argv[1]) == std::string("test"))
			{
				std::cout<<"TESTING ACKNOWLEDGED"<<std::endl;
				return 0;		
			}
		}
		if(argc > 2)
		{
			slurmArrayIndex = atoi(argv[2]);
			std::cout<<"SLURM_ARRAY_TASK_ID = "<<slurmArrayIndex<<std::endl;
		}	
		return 1;
	}

    // std::string inputOutput::initOutputFiles(struct initParams &p, std::shared_ptr<jNonlinearVariationalSolver> solver)
    std::string inputoutput::initOutputFiles(struct initParams &p)
    {

	  // sstream <<  dolfinVC;//dummy class to outpu dolfin environment
	//=========================

		std::string dateString = __TIME__;
		std::replace( dateString.begin(), dateString.end(), ':', '_'); // replace all 'x' to 'y'
		dateString += "-"; 
		dateString += __DATE__;
		std::replace( dateString.begin(), dateString.end(), ' ', '_'); // replace all 'x' to 'y'

		std::stringstream sstream;
	// std::cout<<dateString<<std::endl; 

   // sstream << logString;//10/18/17:  added description string.

    sstream << dateString<<std::endl; 
    sstream<<"TimeStepSize: "<<p.TimeStepSize<<std::endl;
    sstream<<"MaximumTimeSteps: "<<p.MaximumTimeSteps<<std::endl; 
    sstream<<"g_slurmArrayIndex: "<<slurmArrayIndex<<std::endl;
    sstream<<"Alpha: "<<p.Alpha<<std::endl;
    sstream<<"Beta: "<<p.Beta<<std::endl; 
    sstream<<"Delta: "<<p.Delta<<std::endl;
    
    std::cout << sstream.str();
    #define IMAGEFILEPATH  "./images"
    #define LOGFILEPATH "./logs"
    std::string fbase;
    std::string flogpath;
		if(isAugsbergCluster){    
		  fbase.assign("/homes/gast/rahul/images");  
		  flogpath.assign("/homes/gast/rahul/logs"); 
		  fbase += "/";  fbase += dateString;  
		  if(slurmArrayIndex > -1)
		    fbase += ("-" + std::to_string(slurmArrayIndex));  
		    flogpath += ("-" + std::to_string(slurmArrayIndex));            
		  fbase += "/";
		  flogpath += "/";
		}   
		else if(isOpuntiaCluster){    
		  fbase.assign("/home/cpbhanda/images");  
		  fbase += "/";  fbase += dateString;  
		  if(slurmArrayIndex > -1)
		    fbase += ("-" + std::to_string(slurmArrayIndex));      
		  //fbase += ("_" + std::to_string(p.timeSinceEpoch));   
		  fbase += "/";
		}   
		else
		{
			fbase.assign(IMAGEFILEPATH);  
			fbase += "/";  fbase += dateString;  
		//	fbase += ("_" + std::to_string(p.timeSinceEpoch));      
			fbase += "/";
			flogpath.assign(LOGFILEPATH);
		    flogpath += "/"; flogpath += dateString; fbase += "/";
		}

	  
		fpath.assign(fbase + "logFile.txt");
	    logFile.open(fpath, std::ios::trunc);
	    	logFile << sstream.str();
	    logFile.flush();
	    logFile.close();

	    return fbase;
	}

	void inputoutput::finalizeLogFile(double timingData[][5])
	{	
	    logFile.open(fpath, std::ios::app);
	    for (int i=0; i< timingCounter;i++)
	    {
	      for (int j=0; j<5; j++)
	        logFile << timingData[i][j] << ",";
	      logFile<<std::endl;
	    }
	    logFile.flush();
	    logFile.close();
	}

  //lambda-function for recording data to array:
  void inputoutput::logTimingData(double t, double steps, double tstep, double dx0, double contraction1,
  			double timingData[][5])
  {
    timingData[timingCounter][0] = t; 
    timingData[timingCounter][1] = steps; 
    timingData[timingCounter][2] = tstep; 
    timingData[timingCounter][3] = dx0; 
    timingData[timingCounter][4] = contraction1; 
    timingCounter++;
  }


