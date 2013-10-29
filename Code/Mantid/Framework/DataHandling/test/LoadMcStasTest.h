
#ifndef LOADMCSTASTEST_H_
#define LOADMCSTASTEST_H_

#include <fstream>
#include <cxxtest/TestSuite.h>

#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidDataHandling/LoadMcStas.h"
// These includes seem to make the difference between initialization of the
// workspace names (workspace2D/1D etc), instrument classes and not for this test case.
#include "MantidDataObjects/WorkspaceSingleValue.h" 
#include "MantidDataHandling/LoadInstrument.h" 
#include <Poco/Path.h>

using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::DataHandling;
using namespace Mantid::DataObjects;

//
// Test checks if number  of workspace equals one
// Test checks if number  getNumberHistograms = 327682x16384. (128x128= 16384 pixels in one detector)
// 
class LoadMcStasTest : public CxxTest::TestSuite
{
public: 
  
  void testInit()
  {
    TS_ASSERT_THROWS_NOTHING(algToBeTested.initialize());
    TS_ASSERT( algToBeTested.isInitialized() );
  }
  

  void testExec()
  {
    if ( !algToBeTested.isInitialized() ) algToBeTested.initialize();
  
    outputSpace="LoadMcStasTest";
    algToBeTested.setPropertyValue("OutputWorkspace", outputSpace);     
    
    // Should fail because mandatory parameter has not been set
    TS_ASSERT_THROWS(algToBeTested.execute(),std::runtime_error);
        
    
    // Now set it... 
    // specify name of file to load workspace from
    inputFile = "mcstas_event_hist.h5";
    algToBeTested.setPropertyValue("Filename", inputFile);

   
    TS_ASSERT_THROWS_NOTHING(algToBeTested.execute());    
    TS_ASSERT( algToBeTested.isExecuted() );
    //
    //  test workspace created by LoadMcStas
    WorkspaceGroup_sptr output = AnalysisDataService::Instance().retrieveWS<WorkspaceGroup>(outputSpace);
    TS_ASSERT_EQUALS( output->getNumberOfEntries(), 5); // 5 NXdata groups
    //
    //
    MatrixWorkspace_sptr outputItem1 = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(outputSpace+"_1");
    TS_ASSERT_EQUALS( outputItem1->getNumberHistograms(), 8192);  
    //
    //
    MatrixWorkspace_sptr outputItem2 = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(outputSpace+"_2");
    TS_ASSERT_EQUALS( outputItem2->getNumberHistograms(), 1);  
    TS_ASSERT_EQUALS( outputItem2->getNPoints(), 1000); 
    //
    //
    MatrixWorkspace_sptr outputItem3 = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(outputSpace+"_3");
    TS_ASSERT_EQUALS( outputItem3->getNumberHistograms(), 128);   
    //
    //
    MatrixWorkspace_sptr outputItem4 = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(outputSpace+"_4");
    TS_ASSERT_EQUALS( outputItem4->getNumberHistograms(), 1);  
    TS_ASSERT_EQUALS( outputItem4->getNPoints(), 100);
    //
    //
    MatrixWorkspace_sptr outputItem5 = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(outputSpace+"_5");
    TS_ASSERT_EQUALS( outputItem5->getNumberHistograms(), 1);  
    TS_ASSERT_EQUALS( outputItem5->getNPoints(), 100);     
  } // testExec()


 
private:
  LoadMcStas algToBeTested;
  std::string inputFile;
  std::string outputSpace;

};

#endif /*LoadMcStasTEST_H_*/