// -----------------------------------------------------------------------------
// -----                                                                   -----
// -----                           R3BCaloUnpack                           -----
// -----                           Version 1.0                             -----
// -----                    Created  11/10/2013 by Y.Gonzalez              -----
// -----                    Modified 03/03/2014 by M. Bendel               -----
// -----                                                                   -----
// -----------------------------------------------------------------------------

//ROOT headers
#include "TClonesArray.h"

//Fair headers
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairLogger.h"

#include <iomanip>

//Califa headers
#include "R3BCaloRawHit.h"
#include "R3BCaloUnpack.h"

//R3BCaloUnpack: Constructor
R3BCaloUnpack::R3BCaloUnpack(char *strCalDir,
                             Short_t type, Short_t subType,
                             Short_t procId,
                             Short_t subCrate, Short_t control)
  : FairUnpack(type, subType, procId, subCrate, control),
    fRawData(new TClonesArray("R3BCaloRawHit")),
    fNHits(0),fCaloUnpackPar(0)
{
}



//Virtual R3BCaloUnpack: Public method
R3BCaloUnpack::~R3BCaloUnpack()
{
  LOG(INFO) << "R3BCaloUnpack: Delete instance" << FairLogger::endl;
  delete fRawData;
}



//Init: Public method
Bool_t R3BCaloUnpack::Init()
{
  Register();
  return kTRUE;
}


void R3BCaloUnpack::SetParContainers()
{
  // Get run and runtime database
  FairRunAna* run = FairRunAna::Instance();
  if (!run) Fatal("R3BCaloUnpack::SetParContainers", "No analysis run");
  
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  if (!rtdb) Fatal("R3BCaloUnpack::SetParContainers", "No runtime database");
  
  fCaloUnpackPar = (R3BCaloUnpackPar*)(rtdb->getContainer("R3BCaloUnpackPar"));
  
  if ( fCaloUnpackPar ) {
    LOG(INFO) << "R3BCaloUnpack::SetParContainers() "<< FairLogger::endl;
    LOG(INFO) << "Container R3BCaloUnpackPar loaded " << FairLogger::endl;
  }  
}


//Register: Protected method
void R3BCaloUnpack::Register()
{
  LOG(DEBUG) << "Registering" << FairLogger::endl;
  FairRootManager *fMan = FairRootManager::Instance();
  if(! fMan) {
    return;
  }
  fMan->Register("CaloRawHit", "Raw data from Califa", fRawData, kTRUE);
}


//DoUnpack: Public method
Bool_t R3BCaloUnpack::DoUnpack(Int_t *data, Int_t size)
{
  //	trace_head_t *califa_trace = new trace_head_t;
  //	UShort_t *trace_data;
  
  UInt_t l_s = 0;
  UInt_t *pl_data = (UInt_t*) data;

  LOG(DEBUG) << "Unpacking" << FairLogger::endl;

  while(l_s < size) {

    // Loops over CALIFA halves have to be implemented!

    // Remove 0xadd... words
    while((data[l_s] & 0xfff00000) == 0xadd00000) {
      l_s++;
    }
    
    LOG(DEBUG) << "At GOSIP memory" << FairLogger::endl;
    
    UShort_t header_size;
    UShort_t trigger;
    UShort_t pc_id = 0;     //Should be implemented for use of two CALIFA halves
    UShort_t sfp_id = 0;
    UShort_t module_id;
    UShort_t submemory_id;
    UInt_t   data_size;

    header_size = pl_data[l_s] & 0xff;
    trigger = (pl_data[l_s] >> 8) & 0xff;
    module_id = (pl_data[l_s] >> 16) & 0xff;
    submemory_id = (pl_data[l_s++] >> 24) & 0xff;
    data_size = pl_data[l_s++];
    
    LOG(DEBUG) << "========== gosip sub " << FairLogger::endl
    << "     header_size " << header_size << FairLogger::endl
    << "         trigger " << trigger << FairLogger::endl
    << "       module_id " << module_id << FairLogger::endl
    << "    submemory_id " << submemory_id << FairLogger::endl
    << "       data_size " << data_size << FairLogger::endl;

    
    if(header_size != 0x34) {
      break;
    }
    
    // Data reduction: size == 0 -> no more events
    if(data_size == 0) {
      continue;
    }
    
    // special channel
    if(submemory_id == 0xff)
    {
      //!! Prepared for use with different SFPs    !!
      //!! not available in current data structure !!
      // l_s++;
      // sfp_id = (pl_data[l_s++] >> 24) & 0xff;
      l_s += data_size / 4;
      continue;
    }
    
    // Real CALIFA data channel
    UShort_t evsize;
    UShort_t magic_affe;
    UInt_t event_id;
    ULong_t timestamp;
    UShort_t cfd_samples[4];
    UShort_t loverflow;
    UShort_t hoverflow;
    UInt_t   overflow;
    UShort_t self_triggered;
    UShort_t num_pileup;
    UShort_t num_discarded;
    UShort_t energy;
    //    UShort_t reserved;
    UShort_t qpid_size;
    UShort_t magic_babe;
    Short_t n_f;
    UShort_t n_s;
    UChar_t error = 0;
    evsize = pl_data[l_s] & 0xffff;
    magic_affe = (pl_data[l_s++] >> 16) & 0xffff;
    event_id = pl_data[l_s++];
    timestamp = pl_data[l_s++];
    timestamp |= (ULong_t)pl_data[l_s++] << 32;
    cfd_samples[0] = pl_data[l_s] & 0xff;
    cfd_samples[1] = pl_data[l_s++] >> 16;
    cfd_samples[2] = pl_data[l_s] & 0xff;
    cfd_samples[3] = pl_data[l_s++] >> 16;
    overflow = pl_data[l_s] & 0xffffff;
    self_triggered = (pl_data[l_s++] >> 24) & 0xff;
    num_pileup = pl_data[l_s] & 0xffff;
    num_discarded = (pl_data[l_s++] >> 16) & 0xffff;
    energy = pl_data[l_s++] & 0xffff;
    qpid_size = pl_data[l_s] & 0xffff;
    magic_babe = (pl_data[l_s++] >> 16) & 0xffff;
    n_f = pl_data[l_s] & 0xffff;
    n_s = (pl_data[l_s++] >> 16) & 0xffff;

    // Set error bits
    // Error flags: [Pileup][PID][Energy][Timing]

    // Timing not valid
    if (overflow & 0x601)  //11000000001
      error |= 1;
    // Energy not valid
    if (overflow & 0x63e)  //11000111110
      error |= 1<<1;
    // PID not valid
    if (overflow & 0x78E)  //11110001110
      error |= 1<<2;
    if (num_pileup)
      error |= 1<<3;
 
    // Generate crystalID
    UShort_t crystal_id = pc_id * (max_submemory_id*max_module_id*max_sfp_id)
      + sfp_id * (max_submemory_id*max_module_id)
      + module_id * max_submemory_id
      + submemory_id;

    LOG(DEBUG) << " --------- event " << FairLogger::endl
    << "        event_id " << event_id << FairLogger::endl
    << "          energy " << energy << FairLogger::endl
    << "       timestamp " << timestamp << FairLogger::endl
    << "=================================" << FairLogger::endl;
    
    new ((*fRawData)[fNHits]) R3BCaloRawHit(crystal_id, 
					    energy, n_f, n_s, timestamp,
					    error);
    fNHits++;
  }
  
  
  LOG(DEBUG) << "End of memory" << FairLogger::endl;
  LOG(DEBUG) << "R3BCaloUnpack: Number of CALIFA raw hits: " << fNHits << FairLogger::endl;
  
  
  return kTRUE;
}



//Reset: Public method
void R3BCaloUnpack::Reset()
{
  LOG(DEBUG) << "Clearing Data Structure" << FairLogger::endl;
  fRawData->Clear();
  fNHits = 0;
}



ClassImp(R3BCaloUnpack)
