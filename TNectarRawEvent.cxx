#include "TNectarRawEvent.h"

#include "TGo4Log.h"

//************************************************************************//
// BASE class of all subboard structures
TNectarBoard::TNectarBoard() :
    TGo4EventElement(), fUniqueId(0)
{

}
TNectarBoard::TNectarBoard(const char* name, UInt_t unid, Short_t index) :
    TGo4EventElement(name, name, index), fUniqueId(unid)
{
  TGo4Log::Info("TNectarBoard: Create instance %s with unique id: %d, index: %d", name, unid, index);

}

TNectarBoard::~TNectarBoard()
{

}


void TNectarBoard::Clear(Option_t *t)
{
 // NOP, implement in subclass

}

////////////////////////////////////////////////////////////////////

TVmmrSlave::TVmmrSlave(): fDeltaTGate(-1)
   {
     fMessages.clear();
  }

TVmmrSlave::~TVmmrSlave()
   {
     Clear();
   }



void TVmmrSlave::Clear(Option_t*)
{

  for (unsigned i = 0; i < fMessages.size(); ++i)
     {
       delete fMessages[i];
     }
  fMessages.clear();
}


  void TVmmrSlave::AddAdcMessage(TVmmrAdcData* msg)
  {
    fMessages.push_back(msg);
  }


    UInt_t TVmmrSlave::NumAdcMessages()
    {
      return fMessages.size();
    }

    TVmmrAdcData* TVmmrSlave::GetAdcMessage(UInt_t i)
    {
      if(i>=fMessages.size()) return 0;
      return fMessages[i];
    }








// subdata for one single VMMR board
TVmmrBoard::TVmmrBoard() : TNectarBoard()
{
  fSlaves.resize(VMMR_CHAINS,0);
}
TVmmrBoard::TVmmrBoard(const char* name, UInt_t unid, Short_t index) :
    TNectarBoard(name, unid, index)
{
  TGo4Log::Info("TVmmrBoard: Create instance %s with unique id: %d, index: %d", name, unid, index);
  fSlaves.resize(VMMR_CHAINS,0); // provide empty slots for each possible slave board
}

TVmmrBoard::~TVmmrBoard()
{

  for (UInt_t sl = 0; sl < fSlaves.size(); ++sl)
  {
    TVmmrSlave* slave = fSlaves[sl];
    if (slave)
      delete slave;
  }
  fSlaves.clear();
}

TVmmrSlave* TVmmrBoard::GetSlave(UInt_t id,  Bool_t addautomatic)
{
  if (id >= VMMR_CHAINS)
    return 0;
  TVmmrSlave* theSlave = 0;
  theSlave = fSlaves[id];
  if ((theSlave == 0) && addautomatic)
    return AddSlave(id);
  return theSlave;
}

UInt_t TVmmrBoard::GetSlaveNumber()
{
  UInt_t fSlaveNum=0;
  for(UInt_t i=0 ; i < fSlaves.size() ;i++){
    if((long) fSlaves[i]>0){
      fSlaveNum++;
    }
  }
  return fSlaveNum;  
}
TVmmrSlave* TVmmrBoard::AddSlave(UInt_t id)
{
  if (id >= VMMR_CHAINS)
    return 0;
  TVmmrSlave* theSlave = new TVmmrSlave();
  fSlaves[id] = theSlave;
  return theSlave;

}

void TVmmrBoard::Clear(Option_t *t)
{
  for (UInt_t sl = 0; sl < fSlaves.size(); ++sl)
  {
    TVmmrSlave* slave = fSlaves[sl];
    if (slave)
      slave->Clear(t);
  }
}


////////////////////////////////////////////////////////


// sub data for one single MDPP board
TMdppBoard::TMdppBoard() : TNectarBoard()
{

}
TMdppBoard::TMdppBoard(const char* name, UInt_t unid, Short_t index) :
    TNectarBoard(name, unid, index)
{
  TGo4Log::Info("TMdppBoard: Create instance %s with unique id: %d, index: %d", name, unid, index);

}

TMdppBoard::~TMdppBoard()
{
  Clear();
}


void TMdppBoard::AddAdcMessage(TMdppAdcData* msg, UChar_t channel)
{
  if (channel < MDPP_CHANNELS)
    {
      fAdcMessages[channel].push_back(msg);
    }
}

 void TMdppBoard::AddTdcMessage(TMdppTdcData* msg, UChar_t channel)
 {

   if (channel < MDPP_CHANNELS)
       {
         fTdcMessages[channel].push_back(msg);
       }
 }

 void TMdppBoard::AddExtDTMessage(TMdppTdcData* msg, UChar_t trigchan)
   {
     if (trigchan < MDPP_EXTDTCHANNELS)
            {
             fExtTrigDT[trigchan].push_back(msg);
            }
   }


 UInt_t TMdppBoard::NumTdcMessages(UChar_t channel)
   {
       if(channel>=MDPP_CHANNELS) return 0;
       return fTdcMessages[channel].size();
   }

   TMdppTdcData* TMdppBoard::GetTdcMessage(UChar_t channel, UInt_t i)
   {
     if(channel>=MDPP_CHANNELS) return 0;
     if(i>=fTdcMessages[channel].size()) return 0;
     return fTdcMessages[channel].at(i);
   }




     UInt_t TMdppBoard::NumExtDTMessages(UChar_t trigchan)
     {
       if(trigchan>=MDPP_EXTDTCHANNELS) return 0;
             return fExtTrigDT[trigchan].size();
     }

     TMdppTdcData* TMdppBoard::GetExtDtMessage(UChar_t trigchan, UInt_t i)
     {
         if(trigchan>=MDPP_EXTDTCHANNELS) return 0;
            if(i>=fExtTrigDT[trigchan].size()) return 0;
            return fTdcMessages[trigchan].at(i);
     }




  UInt_t TMdppBoard::NumAdcMessages(UChar_t channel)
  {
    if(channel>=MDPP_CHANNELS) return 0;
    return fAdcMessages[channel].size();
  }

  TMdppAdcData* TMdppBoard::GetAdcMessage(UChar_t channel, UInt_t i)
  {
    if(channel>=MDPP_CHANNELS) return 0;
    if(i>=fAdcMessages[channel].size()) return 0;
    return fAdcMessages[channel].at(i);
  }

  Short_t TMdppBoard::GetAdcResolution()
  {
    return (fResolution & 0xB);
  }
  Short_t TMdppBoard::GetTdcResolution()
   {
     return ((fResolution >> 3) & 0xB);
   }






void TMdppBoard::Clear(Option_t *t)
{
  for (Int_t ch = 0; ch < MDPP_CHANNELS; ++ch)
  {

    for (unsigned i = 0; i < fAdcMessages[ch].size(); ++i)
    {
      delete fAdcMessages[ch][i];
    }
    fAdcMessages[ch].clear();

    for (unsigned i = 0; i < fTdcMessages[ch].size(); ++i)
    {
      delete fTdcMessages[ch][i];
    }
    fTdcMessages[ch].clear();



  }
  for (Int_t tc = 0; tc < MDPP_EXTDTCHANNELS; ++tc)
   {
    for (unsigned i = 0; i < fExtTrigDT[tc].size(); ++i)
       {
         delete fExtTrigDT[tc][i];
       }
    fExtTrigDT[tc].clear();
   }

}







//***********************************************************

std::vector<UInt_t> TNectarRawEvent::fgConfigMdppBoards;

std::vector<UInt_t> TNectarRawEvent::fgConfigVmmrBoards;

TNectarRawEvent::TNectarRawEvent() :
    TGo4CompositeEvent(), fSequenceNumber(-1)
{
}
//***********************************************************
TNectarRawEvent::TNectarRawEvent(const char* name, Short_t id) :
    TGo4CompositeEvent(name, name, id), fSequenceNumber(-1)
{
  TGo4Log::Info("TNectarRawEvent: Create instance %s with composite id %d", name, id);
  TString modname;
  UInt_t uniqueid;
  UInt_t indexoff=0;
  for (unsigned i = 0; i < TNectarRawEvent::fgConfigMdppBoards.size(); ++i)
  {
    uniqueid = TNectarRawEvent::fgConfigMdppBoards[i];
    modname.Form("MDPP_Board_%02d", uniqueid);
    addEventElement(new TMdppBoard(modname.Data(), uniqueid, i));
    indexoff++;
  }
  // continue composite index:
  for (unsigned i = 0; i < TNectarRawEvent::fgConfigVmmrBoards.size(); ++i)
    {
      uniqueid = TNectarRawEvent::fgConfigVmmrBoards[i] + VMMR_COMPOSITE_ID_OFFSET;
      modname.Form("VMMR_Board_%02d", uniqueid);
      addEventElement(new TVmmrBoard(modname.Data(), uniqueid, i + indexoff));
    }



}
//***********************************************************
TNectarRawEvent::~TNectarRawEvent()
{
}

TNectarBoard* TNectarRawEvent::GetBoard(UInt_t id)
{
  TNectarBoard* theBoard = 0;
  Short_t numBoards = getNElements();
  for (int i = 0; i < numBoards; ++i)
  {
    theBoard = (TNectarBoard*) getEventElement(i);
    if (theBoard->GetBoardId() == id)
      return theBoard;

  }
  return 0;
}


TVmmrBoard* TNectarRawEvent::GetVmmrBoard(UInt_t uniqueid)
    {
      TVmmrBoard* rev=dynamic_cast<TVmmrBoard*>(GetBoard(uniqueid + VMMR_COMPOSITE_ID_OFFSET));
      return rev;
    }

TMdppBoard* TNectarRawEvent::GetMdppBoard(UInt_t uniqueid)
    {
      TMdppBoard* rev=dynamic_cast<TMdppBoard*>(GetBoard(uniqueid));
      return rev;
    }



void TNectarRawEvent::Clear(Option_t *t)
{
  //TGo4Log::Info("TNectarRawEvent: Clear ");
  TGo4CompositeEvent::Clear();
  fSequenceNumber = -1;
}
