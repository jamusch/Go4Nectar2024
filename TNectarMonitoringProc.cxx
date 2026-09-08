#include "Riostream.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TGo4Analysis.h"
#include "TGo4Log.h"
#include "TGo4Picture.h"
#include "TGo4MbsEvent.h"
#include "TGo4UserException.h"
#include "TGo4WinCond.h"
#include "TGo4ShapedCond.h"
#include "TGo4ListCond.h"
#include "TGo4CondArray.h"

#include "TNectarRawEvent.h"
#include "TNectarRawParam.h"
#include "TNectarDisplay.h"
#include "TNectarRawProc.h"
#include "TNectarMonitoringEvent.h"
#include "TNectarMonitoringDisplay.h"
#include "TNectarMonitoringProc.h"

#include "TGo4Fitter.h"
#include "TGo4FitDataHistogram.h"
#include "TGo4FitParameter.h"
#include "TGo4FitModelPolynom.h"

#include "TGo4UserException.h"

#include "TCutG.h"
#include "TGo4WinCond.h"
#include "TGo4ShapedCond.h"
#include "TGo4CondArray.h"


//*********************************************************
TNectarMonitoringProc::TNectarMonitoringProc() :
  TGo4EventProcessor()
{
  TGo4Log::Info("TNectarMonitoringProc: Create");

}


//*********************************************************
TNectarMonitoringProc::TNectarMonitoringProc(const char* name) :
  TGo4EventProcessor(name), fPar(0)
{
    SetMakeWithAutosave(kTRUE);
    TGo4Log::Info("TNectarMonitoringProc: Create1");
    
    //* Initialize user analysis object *//
    fPar = dynamic_cast<TNectarRawParam*>(MakeParameter("NectarRawParam", "TNectarRawParam", "set_NectarRawParam.C"));
    Histo->InitHistograms(false);

    //* Condition Banana selection *//
    //Double_t cutpnts[11][2]={{-10,9000},{10000,9000},{20000,9000},{30000,9000},{40000,9000},{40000,0},{30000,0},{20000,0},{10000,0},{0,0},{-10,9000}};
    Double_t cutpnts[12][2]={{0,8500},{5000,6000},{15000,4000},{25000,3000},{40000,3000},{40000,1000},{30000,1000},{20000,1500},{10000,3000},{1000,6000},{0,7000},{0,8500}};
    fPolyCon1 = MakePolyCond("polycon", 12, cutpnts, "h_DE_E_tot_keV");

    //* Strainless steel window correction *//
    //TF1* SW_Correction = new TF1("SW_Correction","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)",1,10);
    SW_Correction->SetParameter(0,0.0127746);
    SW_Correction->SetParameter(1,0.150414);
    SW_Correction->SetParameter(2,-0.00560874);
    SW_Correction->SetParameter(3,-0.00148719);
    SW_Correction->SetParameter(4,0.000180784);

}

//***********************************************************
//* Event builder : Fills detectors and updates histograms *//
Bool_t TNectarMonitoringProc::BuildEvent(TGo4EventElement *target)
{
  
    TNectarDetectors* Telescope = (TNectarDetectors*) target;
    //TNectarProfileEvent* Composite = (TNectarProfileEvent*) target;
    //TNectarTelescope* Telescope = (TNectarTelescope*) Composite;
    TNectarRawEvent* rawEvent = (TNectarRawEvent*) GetInputEvent(); 
    
    //* Variable initialization *//
   // Bool_t fOutput = kFALSE; //not store the output
    //Int_t counter = 0;
    Telescope->DT_telescope=0;
    
    for(Int_t i=0;i<17;i++)
    {
        Telescope->BB8_V_Energy[i]=0;
        Telescope->BB8_H_Energy[i]=0;
        Telescope->BB8_VStrip[i]=0;
        Telescope->BB8_HStrip[i]=0;
    }
    
    Telescope->ADC_MSX04_1=0;
    Telescope->ADC_MSX04_2=0;
    Telescope->ADC_MSX04_7=0;
    Telescope->ADC_MSX04_8=0;
    Telescope->ADC_MSX04_9=0;
    Telescope->ADC_MSX04_10=0;
  
    for(Int_t i=0;i<129;i++)
    {
        Telescope->BB29_VStrip[i]=0;
    }
    
    Telescope->DT_BB29_V=0;
    
    for(Int_t i=0;i<65;i++)
    {
        Telescope->BB29_HStrip[i]=0;
    }
    
    Telescope->DT_BB29_H=0;
    Telescope->Trigger = rawEvent->ftrigger;
  
    //* Loop to fill detectors depending on the VMMR id *//
    for (UInt_t i = 0; i < rawEvent->fgConfigVmmrBoards.size(); ++i)
    {      
        UInt_t uniqueid = rawEvent->fgConfigVmmrBoards[i];
        TVmmrBoard* theVmmr = rawEvent->GetVmmrBoard(uniqueid);
        
        /* MOFIFICATION MARCH 2024 
        Previous if condition : if(theVmmr->GetSlaveNumber()==5)
        Checks for the number of buses in the VMMR ! if 5 buses then VMMR0 - if 2 buses then VMMR1
        //TGo4Log::Info("number board %d",uniqueid);
        New condition checks for board id instead
        */
        if(uniqueid==0)
        {
        //TGo4Log::Info("TNectarMonitoring: Telescope_Vmmr %d ", counter); 
        //TGo4Log::Info("TNectarMonitoring: Uniqueid %d ", uniqueid); 
        Filling_Ev_Telescope(theVmmr,Telescope);
        Filling_Ev_FFtopbot(theVmmr,Telescope);
        Filling_Ev_FFside(theVmmr,Telescope);
        }

        if(uniqueid==1)
        {
        //TGo4Log::Info("TNectarMonitoring: BB29_Vmmr %d ", counter);  
        Filling_Ev_BB29(theVmmr,Telescope);
        }
      //counter++;
    }

    Telescope->SetValid(kTRUE);
  
    /* MOFIFICATION MARCH 2024 
    Checks for the type of trigger (target or not) 
    The condition has been removed for tests without target

    if((Telescope->Trigger==1))  Update_Histograms_Target_on(Telescope);
    if((Telescope->Trigger==2))  Update_Histograms_Target_off(Telescope);
        
    if((Telescope->Trigger==1))  Online_Analysis_Target_on(Telescope);  
    if((Telescope->Trigger==2))  Online_Analysis_Target_off(Telescope); */
  
    Update_Histograms_Target_on(Telescope);
    //Update_Histograms_Target_off(Telescope);
    Online_Analysis_Target_on(Telescope);
    //Online_Analysis_Target_off(Telescope);

  return kTRUE;	
}	

//***********************************************************
//* Telescope filling function *//
void TNectarMonitoringProc:: Filling_Ev_Telescope(TVmmrBoard* theVmmr, TNectarDetectors* Telescope)
{
    Double_t ADC1=0,ADC2=0,ADC3=0,ADC4=0,ADC5=0,ADC6=0; //,counter=0;
    
    for (int slid=0; slid<VMMR_CHAINS; ++slid)
    {
        TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
        if((long) theslave>0)
        {
            /* MOFIFICATION MARCH 2024 
            Previous if condition :             
            if(maxvmmrmessages==64){ 
            Does not work for E028 since there is another MMR64 for fission on the first VMMR 
            New condition looks for bus id instead */
            
            UInt_t maxvmmrmessages = theslave->NumAdcMessages();
            //* TELESCOPE *//
            if(slid==0)
            {
                //For cycle around the strips values to fill new variables
                //Conversion  Array VMMR channel -> Strips number 
                int Vmmr_strip_C[32] = {15,11,16,14,12,9,13,8,10,5,6,3,7,2,4,1,15,16,10,14,13,12,9,11,4,8,7,6,3,5,1,2};   
                
                for(UInt_t a=0;a<64;a++)
                { 
                    for (UInt_t i = 0; i < maxvmmrmessages; i++)
                    {
                        TVmmrAdcData* dat = theslave->GetAdcMessage(i);
        
                        if(dat->fFrontEndSubaddress==a)
                        {   
                            // Horizontal strips
                            if(a<=15)   // Filling BB8 MonitoringEvent with Convertion (VMMR channel -> Strips number) 
                            {  
                                
                                Int_t id = Vmmr_strip_C[a]; 
                                Telescope->BB8_HStrip[id] = dat->fData;
                                //counter++;
                            }
                            
                            // Vertical strips
                            if ((a>15)&&(a<32))
                            {
                                Int_t id = Vmmr_strip_C[a];
                                Telescope->BB8_VStrip[id] = dat->fData;
                            }                            
                            
                        }
             
                        if (dat->fFrontEndSubaddress==42) ADC1 = dat->fData;    //MSX04->ADC_MSX04_1 = dat->fData;
                        if (dat->fFrontEndSubaddress==43) ADC2 = dat->fData;    //MSX04->ADC_MSX04_2 = dat->fData; 
                        if (dat->fFrontEndSubaddress==44) ADC3 = dat->fData;    //MSX04->ADC_MSX04_7 = dat->fData;
                        if (dat->fFrontEndSubaddress==45) ADC4 = dat->fData;    //MSX04->ADC_MSX04_8 = dat->fData; 
                        if (dat->fFrontEndSubaddress==46) ADC5 = dat->fData;    //MSX04->ADC_MSX04_9 = dat->fData;
                        if (dat->fFrontEndSubaddress==47) ADC6 = dat->fData;    //MSX04->ADC_MSX04_10 = dat->fData;  
                    }
                }
                
                Telescope->DT_telescope = theslave->fDeltaTGate;
                //TGo4Log::Info("Trigger RAW: %d ", Telescope->DT_telescope);     
                Telescope->ADC_MSX04_1 = ADC1;
                Telescope->ADC_MSX04_2 = ADC2;
                Telescope->ADC_MSX04_7 = ADC3;
                Telescope->ADC_MSX04_8 = ADC4;
                Telescope->ADC_MSX04_9 = ADC5;
                Telescope->ADC_MSX04_10 = ADC6;
            }
	
        }
    }
}

//***********************************************************
//* FF top and bottom filling function *//
void TNectarMonitoringProc:: Filling_Ev_FFtopbot(TVmmrBoard* theVmmr, TNectarDetectors* Telescope)
{
    for (int slid=0; slid<VMMR_CHAINS; ++slid)
    {
        TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
        if((long) theslave>0)
        {
            /* MOFIFICATION MARCH 2024 
            Previous if condition :             
            if(maxvmmrmessages==64){ 
            Does not work for E028 since there is another MMR64 for fission on the first VMMR 
            New condition looks for bus id instead */
            
            UInt_t maxvmmrmessages = theslave->NumAdcMessages();
            
            //* TELESCOPE *//
            if(slid==7 && maxvmmrmessages==32)
            {
                //TGo4Log::Info("max vmmr messages fftop: %d ", maxvmmrmessages); 
                //TGo4Log::Info("slave id: %d ", slid); 
                
                //For cycle around the strips values to fill new variables
                //Conversion  Array VMMR channel -> Strips number 
                int Vmmr_strip_C[32] = {15,11,16,14,12,9,13,8,10,5,6,3,7,2,4,1,15,16,10,14,13,12,9,11,4,8,7,6,3,5,1,2};   
                
                for(UInt_t a=0;a<32;a++)
                { 
                    for (UInt_t i = 0; i < maxvmmrmessages; i++)
                    {
                        TVmmrAdcData* dat = theslave->GetAdcMessage(i);
        
                        if(dat->fFrontEndSubaddress==a)
                        {   
                            // Horizontal strips
                            if(a<=15)   // Filling BB8 MonitoringEvent with Convertion (VMMR channel -> Strips number) 
                            {  
                                
                                Int_t id = Vmmr_strip_C[a]; 
                                Telescope->BB36_top_HStrip[id] = dat->fData;
                                //counter++;
                            }
                            
                            // Vertical strips
                            if ((a>15)&&(a<32))
                            {
                                Int_t id = Vmmr_strip_C[a];
                                Telescope->BB36_top_VStrip[id] = dat->fData;
                            }                            
                            
                        }

                    }
                }
                
                Telescope->DT_BB36_top = theslave->fDeltaTGate;
                //TGo4Log::Info("Trigger RAW: %d ", Telescope->DT_telescope);     
            }
            
            if(slid==4 && maxvmmrmessages==32)
            {
                //TGo4Log::Info("max vmmr messages fftop: %d ", maxvmmrmessages); 
                //TGo4Log::Info("slave id: %d ", slid); 
                
                //For cycle around the strips values to fill new variables
                //Conversion  Array VMMR channel -> Strips number 
                int Vmmr_strip_C[32] = {15,11,16,14,12,9,13,8,10,5,6,3,7,2,4,1,15,16,10,14,13,12,9,11,4,8,7,6,3,5,1,2};   
                
                for(UInt_t a=0;a<32;a++)
                { 
                    for (UInt_t i = 0; i < maxvmmrmessages; i++)
                    {
                        TVmmrAdcData* dat = theslave->GetAdcMessage(i);
        
                        if(dat->fFrontEndSubaddress==a)
                        {   
                            // Horizontal strips
                            if(a<=15)   // Filling BB8 MonitoringEvent with Convertion (VMMR channel -> Strips number) 
                            {  
                                
                                Int_t id = Vmmr_strip_C[a]; 
                                Telescope->BB36_bot_HStrip[id] = dat->fData;
                                //counter++;
                            }
                            
                            // Vertical strips
                            if ((a>15)&&(a<32))
                            {
                                Int_t id = Vmmr_strip_C[a];
                                Telescope->BB36_bot_VStrip[id] = dat->fData;
                            }                            
                            
                        }

                    }
                }
                
                Telescope->DT_BB36_bot = theslave->fDeltaTGate;
                //TGo4Log::Info("Trigger RAW: %d ", Telescope->DT_telescope);     
            }
        }
    }
}

//***********************************************************
//* FF side filling function *//
void TNectarMonitoringProc:: Filling_Ev_FFside(TVmmrBoard* theVmmr, TNectarDetectors* Telescope)
{
    // MMR- 128 channels  
    for (int slid=0; slid<VMMR_CHAINS; ++slid)
    {
        TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
        
        if((long) theslave>0)
        {
            UInt_t maxvmmrmessages = theslave->NumAdcMessages();
            
            if(slid==5)
            {
                int Vmmr_strip_C[128] = {0,0,40,14,36,16,30,10,34,12,38,6,44,8,42,2,48,4,46,20,52,18,50,24,56,22,54,28,60,26,58,32,57,31,59,25,53,27,55,21,49,23,51,17,45,19,47,3,41,1,43,7,37,5,33,11,29,9,35,15,39,13,0,0,62,0,66,64,112,110,116,114,120,118,68,122,70,72,74,76,78,80,82,84,88,86,92,90,96,94,100,98,104,102,108,106,107,105,103,101,99,97,95,93,91,89,87,85,81,83,77,79,73,75,69,71,67,121,119,117,115,113,111,109,65,63,61};    //Convension Array VMMR channel -> Strips number 
                
                for(UInt_t a=0;a<128;a++)
                {    
                    for (UInt_t i = 0; i < maxvmmrmessages; i++)
                    {
                        TVmmrAdcData* dat = theslave->GetAdcMessage(i);
                        
                        if(dat->fFrontEndSubaddress==a)
                        {
                            Int_t id = Vmmr_strip_C[a];
		                    Telescope->BB29_side_VStrip[id] = dat->fData;
                            // TGo4Log::Info("TNectarMonitoring V: %d ", dat->fData);
                            //Telescope->BB29_VStrip[a] = dat->fData;
                        }
                    }
                }
                
                Telescope->DT_BB29_side_V = theslave->fDeltaTGate;
            }
     
      
            // MMR - 64 channels
            //TGo4Log::Info("maxvmmrmessages: %d ", maxvmmrmessages);
            if(slid==6)
            {
                TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
                UInt_t maxvmmrmessages = theslave->NumAdcMessages();
                int Vmmr_strip_C[64] = {0,0,0,0,0,0,0,0,0,0,0,0,39,40,37,38,35,36,33,34,31,32,29,30,27,28,25,26,23,24,21,22,19,20,17,18,15,16,13,14,11,12,9,10,7,8,5,6,3,4,1,2,0,0,0,0,0,0,0,0,0,0,0,0};    //Convension Array VMMR channel -> Strips number 
	  
                for(UInt_t a=0;a<64;a++)
                {    
                    for (UInt_t i = 0; i < maxvmmrmessages; i++)
                    {
                        TVmmrAdcData* dat = theslave->GetAdcMessage(i);
                        if(dat->fFrontEndSubaddress==a)
                        {
                            Int_t id = Vmmr_strip_C[a]; 
                            Telescope->BB29_side_HStrip[id] = dat->fData;
                            //Telescope->BB29_HStrip[a] = dat->fData;
                        }
                    }
                }
                
                Telescope->DT_BB29_side_H = theslave->fDeltaTGate;
            }
        }
    }
}

//***********************************************************
//* HR filling function *//
void TNectarMonitoringProc:: Filling_Ev_BB29(TVmmrBoard* theVmmr, TNectarDetectors* Telescope)
{
  // MMR- 128 channels  
  for (int slid=0; slid<VMMR_CHAINS; ++slid)
    {
      TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
      if((long) theslave>0){
	UInt_t maxvmmrmessages = theslave->NumAdcMessages();
	if(slid==5){
	  int Vmmr_strip_C[128] = {0,0,40,14,36,16,30,10,34,12,38,6,44,8,42,2,48,4,46,20,52,18,50,24,56,22,54,28,60,26,58,32,57,31,59,25,53,27,55,21,49,23,51,17,45,19,47,3,41,1,43,7,37,5,33,11,29,9,35,15,39,13,0,0,62,0,66,64,112,110,116,114,120,118,68,122,70,72,74,76,78,80,82,84,88,86,92,90,96,94,100,98,104,102,108,106,107,105,103,101,99,97,95,93,91,89,87,85,81,83,77,79,73,75,69,71,67,121,119,117,115,113,111,109,65,63,61};    //Convension Array VMMR channel -> Strips number 
	  for(UInt_t a=0;a<128;a++){    
	    for (UInt_t i = 0; i < maxvmmrmessages; i++)
	      {
		TVmmrAdcData* dat = theslave->GetAdcMessage(i);
		if(dat->fFrontEndSubaddress==a){
		  Int_t id = Vmmr_strip_C[a];
		  Telescope->BB29_VStrip[id] = dat->fData;
		  // TGo4Log::Info("TNectarMonitoring V: %d ", dat->fData);
		  //Telescope->BB29_VStrip[a] = dat->fData;
		}
	      }
	  }
	  Telescope->DT_BB29_V = theslave->fDeltaTGate;
	}
     
      
	// MMR - 64 channels
	//TGo4Log::Info("maxvmmrmessages: %d ", maxvmmrmessages);
	if(slid==6){
	  TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
	  UInt_t maxvmmrmessages = theslave->NumAdcMessages();
	  int Vmmr_strip_C[64] = {0,0,0,0,0,0,0,0,0,0,0,0,39,40,37,38,35,36,33,34,31,32,29,30,27,28,25,26,23,24,21,22,19,20,17,18,15,16,13,14,11,12,9,10,7,8,5,6,3,4,1,2,0,0,0,0,0,0,0,0,0,0,0,0};    //Convension Array VMMR channel -> Strips number 
	  for(UInt_t a=0;a<64;a++){    
	    for (UInt_t i = 0; i < maxvmmrmessages; i++)
	      {
		TVmmrAdcData* dat = theslave->GetAdcMessage(i);
		if(dat->fFrontEndSubaddress==a){
		  Int_t id = Vmmr_strip_C[a]; 
		  Telescope->BB29_HStrip[id] = dat->fData;
		  //Telescope->BB29_HStrip[a] = dat->fData;
		}
	      }
	  }
	  Telescope->DT_BB29_H = theslave->fDeltaTGate;
	}
      }
    }
}

//***********************************************************
//* Online analysis target on function *//
    void TNectarMonitoringProc:: Online_Analysis_Target_on(TNectarDetectors* Telescope){
    Int_t S_E_id=0, S_E_id_V=0  /*(modified to test mapping BB29)*/, S_E_id_H=0, nV=0, nV1=0, nV2=0, nH=0, Total_ch=0, strip_BB29_V=0, strip_BB29_H=0,  S_E_id_V_HR=0, S_E_id_H_HR=0;
    Double_t Total_E=0, E_1=0, E_2=0, E_7=0, E_8=0, E_9=0, E_10=0; //, E_BB8_Mult=0;
    Int_t mult_BB8_V=0, mult_BB8_H=0;
    Int_t DT_m = Telescope->DT_telescope + fPar->DT_min;
    Int_t DT_M = Telescope->DT_telescope + fPar->DT_max;  
    Int_t DT = Telescope->DT_BB29_V;
  
    fPolyCon1->Enable();
    
    Histo->BB8_V8_H8->Fill(Telescope->BB8_VStrip[8],Telescope->BB8_HStrip[8]);
    
  
    for(Int_t i=1; i<17; i++){
    //Total_ch=0;
    //Total_E=0;
        if(Telescope->BB8_VStrip[i]>fPar->BB8_V_threshold)
        {
            nV=i; S_E_id++; S_E_id_V++;
            Telescope->BB8_V_Energy[nV] = Calib_BB8_V(Telescope->BB8_VStrip[nV],nV) ;
            Histo->BB8_DT_E_V[i-1]->Fill(Telescope->DT_telescope,Telescope->BB8_VStrip[nV]);
            mult_BB8_V++;
            for(Int_t j=1; j<17; j++)
            {
                if(Telescope->BB8_HStrip[j]>fPar->BB8_H_threshold)
                {
                    mult_BB8_H++;
                    nH=j; S_E_id++; S_E_id_H++;
                    Telescope->BB8_H_Energy[nH] = Calib_BB8_H(Telescope->BB8_VStrip[nH],nH);
                    Histo->h_BB8_Map->Fill(nV,-nH);
                    Histo->BB8_1V_allH[nV-1]->Fill(Telescope->BB8_VStrip[nV],Telescope->BB8_HStrip[nH]);
                    
                }	  	  
            }
        } 
    }
  
  // Loop to check bidim energy of adjacent vertical strips
  for(Int_t i=1; i<17; i++){
    if(Telescope->BB8_VStrip[i]>fPar->BB8_V_threshold)
    {
        nV1=i; S_E_id++; S_E_id_V++;
        
        
        for(Int_t j=1; j<17; j++){
            if(Telescope->BB8_VStrip[i]>fPar->BB8_V_threshold)
            {
                nV2=j; S_E_id++; S_E_id_H++;
                
                if (j==i+1){
                    Histo->BB8_adj_V[nV1-1]->Fill(Telescope->BB8_VStrip[nV1],Telescope->BB8_VStrip[nV2]);
                }
                
            }
        }
        
        
    } 
  }
  
    if((mult_BB8_V==1)&&(mult_BB8_H==1)){
      Histo->BB8_V8_H8_mult1->Fill(Telescope->BB8_VStrip[8],Telescope->BB8_HStrip[8]);
    }
  
  // Loop to check bidim energy of adjacent horizontal strips
  for(Int_t i=1; i<17; i++){
    if(Telescope->BB8_HStrip[i]>fPar->BB8_H_threshold)
    {
        nV1=i; S_E_id++; S_E_id_V++;        
        for(Int_t j=1; j<17; j++){
            if(Telescope->BB8_HStrip[i]>fPar->BB8_H_threshold)
            {
                nV2=j; S_E_id++; S_E_id_H++;
                if (j==i+1){
                    Histo->BB8_adj_H[nV1-1]->Fill(Telescope->BB8_HStrip[nV1],Telescope->BB8_HStrip[nV2]);
                }
                
            }
        }
    } 
  }
  
    for(Int_t i=1; i<17; i++){
    //Total_ch=0;
    //Total_E=0;
    if(Telescope->BB36_top_VStrip[i]>fPar->BB8_V_threshold)
    {
        nV=i; S_E_id++; S_E_id_V++;
        Histo->BB36_FFtop_DT_E_V[i-1]->Fill(Telescope->DT_BB36_top,Telescope->BB36_top_VStrip[nV]);
        
        for(Int_t j=1; j<17; j++){
            if(Telescope->BB36_top_HStrip[j]>fPar->BB8_H_threshold)
            {
                nH=j; S_E_id++; S_E_id_H++;
                Histo->h_FFtop_Map->Fill(nV,-nH);                
            }	  	  
        }
    } 
    
    }
    for(Int_t i=1; i<17; i++){
    //Total_ch=0;
    //Total_E=0;
    if(Telescope->BB36_bot_VStrip[i]>fPar->BB8_V_threshold)
    {
        nV=i; S_E_id++; S_E_id_V++;
        Histo->BB36_FFbot_DT_E_V[i-1]->Fill(Telescope->DT_BB36_bot,Telescope->BB36_bot_VStrip[nV]);
        
        for(Int_t j=1; j<17; j++){
            if(Telescope->BB36_bot_HStrip[j]>fPar->BB8_H_threshold)
            {
                nH=j; S_E_id++; S_E_id_H++;
                Histo->h_FFbot_Map->Fill(nV,-nH);                
            }	  	  
        }
    } 
  }

    for(Int_t i=1; i<129; i++){
    //Total_ch=0;
    //Total_E=0;
    if(Telescope->BB29_VStrip[i]>fPar->BB8_V_threshold)
    {
        nV=i; S_E_id++; S_E_id_V++;        
        for(Int_t j=1; j<65; j++){
            if(Telescope->BB29_HStrip[j]>fPar->BB8_H_threshold)
            {
                nH=j; S_E_id++; S_E_id_H++;
                Histo->h_HR_Map->Fill(i,-j);
                
                if (i==1)   Histo->BB29_HR_V1_allH->Fill(Telescope->BB29_VStrip[i],Telescope->BB29_HStrip[j]);
                if (i==2)   Histo->BB29_HR_V2_allH->Fill(Telescope->BB29_VStrip[i],Telescope->BB29_HStrip[j]);
                if (i==60)   Histo->BB29_HR_V60_allH->Fill(Telescope->BB29_VStrip[i],Telescope->BB29_HStrip[j]);
                if (i==61)   Histo->BB29_HR_V61_allH->Fill(Telescope->BB29_VStrip[i],Telescope->BB29_HStrip[j]);
                if (i==119)   Histo->BB29_HR_V119_allH->Fill(Telescope->BB29_VStrip[i],Telescope->BB29_HStrip[j]);
                if (i==120)   Histo->BB29_HR_V120_allH->Fill(Telescope->BB29_VStrip[i],Telescope->BB29_HStrip[j]);
              
            }	  	  
        }
    } 
  }
  
  Histo->BB29_HR_DT_V1->Fill(Telescope->DT_BB29_V,Telescope->BB29_VStrip[1]);
  Histo->BB29_HR_DT_V2->Fill(Telescope->DT_BB29_V,Telescope->BB29_VStrip[2]);
  Histo->BB29_HR_DT_V60->Fill(Telescope->DT_BB29_V,Telescope->BB29_VStrip[60]);
  Histo->BB29_HR_DT_V61->Fill(Telescope->DT_BB29_V,Telescope->BB29_VStrip[61]);
  Histo->BB29_HR_DT_V119->Fill(Telescope->DT_BB29_V,Telescope->BB29_VStrip[119]);
  Histo->BB29_HR_DT_V120->Fill(Telescope->DT_BB29_V,Telescope->BB29_VStrip[120]);
  
  Histo->BB29_HR_DT_H1->Fill(Telescope->DT_BB29_H,Telescope->BB29_HStrip[1]);
  Histo->BB29_HR_DT_H2->Fill(Telescope->DT_BB29_H,Telescope->BB29_HStrip[2]);
  Histo->BB29_HR_DT_H20->Fill(Telescope->DT_BB29_H,Telescope->BB29_HStrip[20]);
  Histo->BB29_HR_DT_H21->Fill(Telescope->DT_BB29_H,Telescope->BB29_HStrip[21]);
  Histo->BB29_HR_DT_H39->Fill(Telescope->DT_BB29_H,Telescope->BB29_HStrip[39]);
  Histo->BB29_HR_DT_H40->Fill(Telescope->DT_BB29_H,Telescope->BB29_HStrip[40]);
  
  Histo->BB29_HR_V1_V2->Fill(Telescope->BB29_VStrip[1],Telescope->BB29_VStrip[2]);
  Histo->BB29_HR_V2_V3->Fill(Telescope->BB29_VStrip[2],Telescope->BB29_VStrip[3]);
  Histo->BB29_HR_V60_V61->Fill(Telescope->BB29_VStrip[60],Telescope->BB29_VStrip[61]);
  Histo->BB29_HR_V61_V62->Fill(Telescope->BB29_VStrip[61],Telescope->BB29_VStrip[62]);
  Histo->BB29_HR_V118_V119->Fill(Telescope->BB29_VStrip[118],Telescope->BB29_VStrip[119]);
  Histo->BB29_HR_V119_V120->Fill(Telescope->BB29_VStrip[119],Telescope->BB29_VStrip[120]);
  
    for(Int_t i=1; i<129; i++){
    //Total_ch=0;
    //Total_E=0;
    if(Telescope->BB29_side_VStrip[i]>fPar->BB8_V_threshold)
    {
        nV=i; S_E_id++; S_E_id_V++;        
        for(Int_t j=1; j<65; j++){
            if(Telescope->BB29_side_HStrip[j]>fPar->BB8_H_threshold)
            {
                nH=j; S_E_id++; S_E_id_H++;
                Histo->h_FFside_Map->Fill(nV,-nH);
                if (i==1)   Histo->BB29_FFside_V1_allH->Fill(Telescope->BB29_side_VStrip[i],Telescope->BB29_side_HStrip[j]);
                if (i==2)   Histo->BB29_FFside_V2_allH->Fill(Telescope->BB29_side_VStrip[i],Telescope->BB29_side_HStrip[j]);
                if (i==60)   Histo->BB29_FFside_V60_allH->Fill(Telescope->BB29_side_VStrip[i],Telescope->BB29_side_HStrip[j]);
                if (i==61)   Histo->BB29_FFside_V61_allH->Fill(Telescope->BB29_side_VStrip[i],Telescope->BB29_side_HStrip[j]);
                if (i==119)   Histo->BB29_FFside_V119_allH->Fill(Telescope->BB29_side_VStrip[i],Telescope->BB29_side_HStrip[j]);
                if (i==120)   Histo->BB29_FFside_V120_allH->Fill(Telescope->BB29_side_VStrip[i],Telescope->BB29_side_HStrip[j]);
                
            }	  	  
        }
    } 
  }
  
  Histo->BB29_FFside_V1_V2->Fill(Telescope->BB29_side_VStrip[1],Telescope->BB29_side_VStrip[2]);
  Histo->BB29_FFside_V2_V3->Fill(Telescope->BB29_side_VStrip[2],Telescope->BB29_side_VStrip[3]);
  Histo->BB29_FFside_V60_V61->Fill(Telescope->BB29_side_VStrip[60],Telescope->BB29_side_VStrip[61]);
  Histo->BB29_FFside_V61_V62->Fill(Telescope->BB29_side_VStrip[61],Telescope->BB29_side_VStrip[62]);
  Histo->BB29_FFside_V118_V119->Fill(Telescope->BB29_side_VStrip[118],Telescope->BB29_side_VStrip[119]);
  Histo->BB29_FFside_V119_V120->Fill(Telescope->BB29_side_VStrip[119],Telescope->BB29_side_VStrip[120]);
  
  Histo->BB29_FFside_DT_V1->Fill(Telescope->DT_BB29_side_V,Telescope->BB29_side_VStrip[1]);
  Histo->BB29_FFside_DT_V2->Fill(Telescope->DT_BB29_side_V,Telescope->BB29_side_VStrip[2]);
  Histo->BB29_FFside_DT_V60->Fill(Telescope->DT_BB29_side_V,Telescope->BB29_side_VStrip[60]);
  Histo->BB29_FFside_DT_V61->Fill(Telescope->DT_BB29_side_V,Telescope->BB29_side_VStrip[61]);
  Histo->BB29_FFside_DT_V119->Fill(Telescope->DT_BB29_side_V,Telescope->BB29_side_VStrip[119]);
  Histo->BB29_FFside_DT_V120->Fill(Telescope->DT_BB29_side_V,Telescope->BB29_side_VStrip[120]);
  
  Histo->BB29_FFside_DT_H1->Fill(Telescope->DT_BB29_side_H,Telescope->BB29_side_HStrip[1]);
  Histo->BB29_FFside_DT_H2->Fill(Telescope->DT_BB29_side_H,Telescope->BB29_side_HStrip[2]);
  Histo->BB29_FFside_DT_H20->Fill(Telescope->DT_BB29_side_H,Telescope->BB29_side_HStrip[20]);
  Histo->BB29_FFside_DT_H21->Fill(Telescope->DT_BB29_side_H,Telescope->BB29_side_HStrip[21]);
  Histo->BB29_FFside_DT_H39->Fill(Telescope->DT_BB29_side_H,Telescope->BB29_side_HStrip[39]);
  Histo->BB29_FFside_DT_H40->Fill(Telescope->DT_BB29_side_H,Telescope->BB29_side_HStrip[40]);
  
  Histo->DTtel_ADCFFtop_V1->Fill(Telescope->DT_telescope,Telescope->BB36_top_VStrip[1]);
  Histo->DTtel_ADCFFtop_V8->Fill(Telescope->DT_telescope,Telescope->BB36_top_VStrip[8]);
  Histo->DTtel_ADCFFbot_V1->Fill(Telescope->DT_telescope,Telescope->BB36_bot_VStrip[1]);
  Histo->DTtel_ADCFFbot_V8->Fill(Telescope->DT_telescope,Telescope->BB36_bot_VStrip[8]);
  Histo->DTtel_ADCFFside_V1->Fill(Telescope->DT_telescope,Telescope->BB29_side_VStrip[1]);
  Histo->DTtel_ADCFFside_V60->Fill(Telescope->DT_telescope,Telescope->BB29_side_VStrip[60]);
  Histo->DTtel_ADCFFside_H1->Fill(Telescope->DT_telescope,Telescope->BB29_side_HStrip[1]);
  Histo->DTtel_ADCFFside_H20->Fill(Telescope->DT_telescope,Telescope->BB29_side_HStrip[20]);  
  

  
    //Condition on Events with Multiplicity larger than 1
  if(((S_E_id_V>1)&&(S_E_id_H>=1))||((S_E_id_V>=1)&&(S_E_id_H>1))){
      
      if(Telescope->ADC_MSX04_1>fPar->MSX04_1_threshold) 
	{
	  Total_ch=Total_ch+(Telescope->ADC_MSX04_1);
	  E_1=Calib_MSX04(Telescope->ADC_MSX04_1, 1);
	  Total_E=Total_E+E_1;
	  if(Telescope->ADC_MSX04_2>fPar->MSX04_2_threshold)
	    {
	      Total_ch=Total_ch+(Telescope->ADC_MSX04_2);
	      E_2=Calib_MSX04(Telescope->ADC_MSX04_2, 2);
	      Total_E=Total_E+E_2;
	      if(Telescope->ADC_MSX04_7>fPar->MSX04_7_threshold)
		{
		  Total_ch=Total_ch+(Telescope->ADC_MSX04_7);
		  E_7=Calib_MSX04(Telescope->ADC_MSX04_7, 7);
		  Total_E=Total_E+E_7;
		  if(Telescope->ADC_MSX04_8>fPar->MSX04_8_threshold)
		    {
		      Total_ch=Total_ch+(Telescope->ADC_MSX04_8);
		      E_8=Calib_MSX04(Telescope->ADC_MSX04_8, 8);
		      Total_E=Total_E+E_8;
		      if(Telescope->ADC_MSX04_9>fPar->MSX04_9_threshold)
			{
			  Total_ch=Total_ch+(Telescope->ADC_MSX04_9);
			  E_9=Calib_MSX04(Telescope->ADC_MSX04_9, 9);
			  Total_E=Total_E+E_9;
			  if(Telescope->ADC_MSX04_10>fPar->MSX04_10_threshold)
			    {
			      Total_ch=Total_ch+(Telescope->ADC_MSX04_10);
			      E_10=Calib_MSX04(Telescope->ADC_MSX04_10, 10);
			      Total_E=Total_E+E_10;
			    }
			}
		    }
		}
	    }
	}
      
      Histo->h_DE_E_tot_Mult->Fill(Total_E,Telescope->BB8_V_Energy[nV]);
      
    }
  
  
  
  
  
  //Condition on the events with 1 V Strip and 1 H Strip on the BB8
  if((S_E_id_V==1)&&(S_E_id_H==1)){ // it must be setted equal to 1 in order to consider only one strip events
    Total_ch=0;
    Total_E=0;
    Histo->h_DE_E_1[nV-1]->Fill(Telescope->ADC_MSX04_1,Telescope->BB8_VStrip[nV]);
    Histo->h_DE_E_2[nV-1]->Fill(Telescope->ADC_MSX04_2,Telescope->BB8_VStrip[nV]);
    Histo->h_DE_E_7[nV-1]->Fill(Telescope->ADC_MSX04_7,Telescope->BB8_VStrip[nV]);
    Histo->h_DE_E_8[nV-1]->Fill(Telescope->ADC_MSX04_8,Telescope->BB8_VStrip[nV]);
    Histo->h_DE_E_9[nV-1]->Fill(Telescope->ADC_MSX04_9,Telescope->BB8_VStrip[nV]);
    Histo->h_DE_E_10[nV-1]->Fill(Telescope->ADC_MSX04_10,Telescope->BB8_VStrip[nV]);
    
    //Energy Conversion??
    
    //DSSD_E=C_DSSSD(Telescope->BB8_VStrip[nV]);
    
    if(Telescope->ADC_MSX04_1>fPar->MSX04_1_threshold) 
      {
	Total_ch=Total_ch+(Telescope->ADC_MSX04_1);
	E_1=Calib_MSX04(Telescope->ADC_MSX04_1, 1);
	Total_E=Total_E+E_1;
	if(Telescope->ADC_MSX04_2>fPar->MSX04_2_threshold)
	  {
	    Total_ch=Total_ch+(Telescope->ADC_MSX04_2);
	    E_2=Calib_MSX04(Telescope->ADC_MSX04_2, 2);
	    Total_E=Total_E+E_2;
	    if(Telescope->ADC_MSX04_7>fPar->MSX04_7_threshold)
	      {
		Total_ch=Total_ch+(Telescope->ADC_MSX04_7);
		E_7=Calib_MSX04(Telescope->ADC_MSX04_7, 7);
		Total_E=Total_E+E_7;
		if(Telescope->ADC_MSX04_8>fPar->MSX04_8_threshold)
		  {
		    Total_ch=Total_ch+(Telescope->ADC_MSX04_8);
		    E_8=Calib_MSX04(Telescope->ADC_MSX04_8, 8);
		    Total_E=Total_E+E_8;
		    if(Telescope->ADC_MSX04_9>fPar->MSX04_9_threshold)
		      {
			Total_ch=Total_ch+(Telescope->ADC_MSX04_9);
			E_9=Calib_MSX04(Telescope->ADC_MSX04_9, 9);
			Total_E=Total_E+E_9;
			if(Telescope->ADC_MSX04_10>fPar->MSX04_10_threshold)
			  {
			    Total_ch=Total_ch+(Telescope->ADC_MSX04_10);
			    E_10=Calib_MSX04(Telescope->ADC_MSX04_10, 10);
			    Total_E=Total_E+E_10;
			  }
		      }
		  }
	      }
	  }
      }
    //DE_E Energy Plots in ch for strip
    Histo->h_DE_E_tot_ch[nV-1]->Fill(Total_ch, Telescope->BB8_VStrip[nV]);
    //DE_E Energy Plots in KeV for strip
    Histo->h_DE_E_tot_E[nV-1]->Fill(Total_E,Telescope->BB8_V_Energy[nV]);
    //BB8_Mapping  
    //Histo->h_BB8_Map->Fill(nV,-nH);
    //Calibration E detectors
      Histo->h_DE_E_1_Calib->Fill(E_1,Telescope->BB8_V_Energy[nV]);
      Histo->h_DE_E_2_Calib->Fill(E_2,Telescope->BB8_V_Energy[nV]);
      Histo->h_DE_E_7_Calib->Fill(E_7,Telescope->BB8_V_Energy[nV]);
      Histo->h_DE_E_8_Calib->Fill(E_8,Telescope->BB8_V_Energy[nV]);
      Histo->h_DE_E_9_Calib->Fill(E_9,Telescope->BB8_V_Energy[nV]);
      Histo->h_DE_E_10_Calib->Fill(E_10,Telescope->BB8_V_Energy[nV]);
    //Strip vs Total Energy plot
    Histo->h_V_E_tot->Fill(nV,Total_E+Telescope->BB8_V_Energy[nV]);
    Histo->h_V_ch_tot->Fill(nV,Total_ch+Telescope->BB8_VStrip[nV]);
    Histo->h_DE_E_tot_keV->Fill(Total_E,Telescope->BB8_V_Energy[nV]);
    Histo->h_DE_E_tot_keV_cond->Fill(Total_E,Telescope->BB8_V_Energy[nV]);
    if(fPolyCon1->Test(Total_E,Telescope->BB8_V_Energy[nV])) Ex_Calculation(nV,nH,Telescope,Total_E);   //Excitation Energy calculation for Banana Selection
  }

  //Multiplicity BB8
 if(S_E_id_V>0) Histo->h_Mult_BB8_V->Fill(S_E_id_V);
 if(S_E_id_H>0) Histo->h_Mult_BB8_H->Fill(S_E_id_H);
  
  
  
  for(Int_t a=1;a<=128;a++){
    if(Telescope->BB29_VStrip[a]>fPar->BB29_V_threshold){
      strip_BB29_V=a; S_E_id_V_HR++;
      for(Int_t b=1;b<=64;b++){
	if(Telescope->BB29_HStrip[b]>fPar->BB29_H_threshold)
	  { 
	    strip_BB29_H=-b; S_E_id_H_HR++;
	    //Histo->h_Map_HR->Fill(strip_BB29_V,strip_BB29_H);
	    if(fPolyCon1->Test(Total_E,Telescope->BB8_V_Energy[nV]))    Histo->h_cond_Map_BB29_DE->Fill(strip_BB29_V,strip_BB29_H);
	    if((DT>DT_m)&&(DT<DT_M))
	      {
		Histo->h_cond_Map_BB29_DT->Fill(strip_BB29_V,strip_BB29_H);  
		if(fPolyCon1->Test(Total_E,Telescope->BB8_V_Energy[nV]))    Histo->h_cond_Map_BB29->Fill(strip_BB29_V,strip_BB29_H);
	      }
	  }
      }
      //BB29 Conditioned Mapping
    }
  }
  //Histo->h_Map_BB29->Fill(strip_BB29_V,strip_BB29_H);
  
  //BB29 Conditioned Mapping
  // if(fPolyCon1->Test(DSSSD_E_V,Total_E))        Histo->h_cond_Map_BB29->Fill(strip_BB29_V,strip_BB29_H); 
  
  
  //Multiplicity BB29
  if(S_E_id_V_HR>0) Histo->h_Mult_BB29_V->Fill(S_E_id_V_HR);
  if(S_E_id_H_HR>0) Histo->h_Mult_BB29_H->Fill(S_E_id_H_HR);
  
  /* Fission */
    
  
  
  
}



void TNectarMonitoringProc:: Online_Analysis_Target_off(TNectarDetectors* Telescope){

  Int_t nV=0, nH=0;
  Int_t DT_m = Telescope->DT_telescope + fPar->DT_min;
  Int_t DT_M = Telescope->DT_telescope + fPar->DT_max;  
  Int_t DT = Telescope->DT_BB29_V;
 
  
  for(Int_t i=1; i<17; i++){
    if(Telescope->BB8_VStrip[i]>fPar->BB8_V_threshold)
      {
	nV=i;
	for(Int_t a=1; a<17; a++){
	  if(Telescope->BB8_HStrip[i]>fPar->BB8_H_threshold)
	    {
	      nH=a;
	      Histo->h_BB8_Map_TO->Fill(nV,-nH);
	    }
	}
      }
  }

  for(Int_t a=0;a<=128;a++){
    if(Telescope->BB29_VStrip[a]>fPar->BB29_V_threshold){
      for(Int_t b=0;b<=64;b++){
	if(Telescope->BB29_HStrip[b]>fPar->BB29_H_threshold)
	  { 
	    Histo->h_Map_BB29_TO->Fill(a,-b);
	    if((DT>DT_m)&&(DT<DT_M)){
	      Histo->h_cond_Map_BB29_DT_TO->Fill(a,-b);  
	    }
	  }
      }
      //BB29 Conditioned Mapping
    }
  }
}



Double_t TNectarMonitoringProc:: Calib_BB8_V(Double_t ADC, Int_t nV)
{
  Double_t Energy_V=0;
  Double_t m_calib_V[16]={4.12604,4.06538,4.15413,4.12212,4.06069,3.92988,3.9986,4.02184,4.10397,4.06809,3.98872,3.97404,4.02592,4.09213,4.06796,4.16795};//Michele Calib
  Double_t q_calib_V[16]={454.266,318.196,474.151,475.37,406.817,300.797,454.835,424.663,458.116,463.144,315.981,418.269,466.007,319.59,478.288,427.156};//Michele Calib
  // Double_t m_calib_V[16]={4.072,4.017,4.114,4.095,4.034,3.913,3.95,3.989,4.06,4.021,3.967,3.949,4.007,4.07,4.024,4.127}; //Camille Calib
  //Double_t q_calib_V[16]={365.995,549.775,399.506,418.903,347.597,255.607,367.192,354.707,375.9,372.696,261.65,363.78,481.929,269.406,400.966,356.675}; //Camille Calib
  
  Energy_V = m_calib_V[nV-1] * ADC - q_calib_V[nV-1];
  return Energy_V;
  }

Double_t TNectarMonitoringProc:: Calib_BB8_H(Double_t ADC, Int_t nH)
{
  Double_t Energy_H=0;
 Double_t m_calib_H[16]={3.92803,3.737,3.73528,3.94152,3.76983,3.94507,3.82148,3.80815,3.79017,3.91643,3.80309,3.75779,3.97353,3.91514,3.79044,3.90423};//Michele Calib
 Double_t q_calib_H[16]={294.683,381.02,292.063,437.912,379.667,436.76,438.463,292.005,383.706,456.686,373.67,430.678,430.87,322.961,396.87,410.609};//Michele Calib
 //Double_t m_calib_H[16]={3.89,3.7,3.702,3.892,3.75,3.902,3.793,3.784,3.757,3.881,3.788,3.732,3.93,3.865,3.79,3.862}; //Camille Calib
 //Double_t q_calib_H[16]={235.689,313.648,230.026,351.356,327.141,351.703,375.287,238.374,317.911,385.829,330.806,369.062,348.021,227.849,377.552,333.667}; //Camille Calib
  
 
Energy_H = m_calib_H[nH-1] * ADC - q_calib_H[nH-1];
return Energy_H;
}

Double_t TNectarMonitoringProc:: Calib_MSX04(Double_t ADC, Int_t N_MSX04)
{

  Double_t Energy=0;
  Double_t m_calib[11]={0,4.23,4.21,0,0,0,0,4.324,4.457,4.334,4.347}; //Camille Calib
  Double_t q_calib[11]={0,326.089,402.757,0,0,0,0,324.731,483.684,332.537,417.8}; //Camille Calib
  
  Energy = m_calib[N_MSX04] * ADC - q_calib[N_MSX04];
  return Energy;


}


void TNectarMonitoringProc:: Ex_Calculation(Int_t n_V, Int_t n_H, TNectarDetectors* Telescope, Double_t Total_E){
  
  Double_t a=0,b=0,c=0, Ex=0, z=96, x1=0,y1=0,z1=0, r1=0, Total_Energy=0, H_Theta=0;
  Double_t tfix = (fPar->tfix*TMath::Pi())/180; 
  Total_Energy = Total_E/1000 + SW_Correction->Eval((Telescope->BB8_V_Energy[n_V])/1000)+Telescope->BB8_V_Energy[n_V]/1000;
  
  //Theta angle calculation
  
  x1 = z*TMath::Sin(tfix)+(fPar->A_Pos[n_V]+fPar->Beam_shift)*TMath::Cos(tfix);
  y1 = fPar->A_Pos[17-n_H];
  z1 = z*TMath::Cos(tfix)-(fPar->A_Pos[n_V]+fPar->Beam_shift)*TMath::Sin(tfix);
  r1 = TMath::Sqrt(TMath::Power(x1,2)+TMath::Power(y1,2));
  H_Theta = TMath::ATan(r1/z1); 
 
  //Excitation Energy Equation 
  
  a = TMath::Power((fPar->E_Pb + fPar->M_Pb)-(Total_Energy+fPar->M_H)+fPar->M_H,2);
  b = TMath::Power(fPar->E_Pb, 2)+2*fPar->E_Pb*fPar->M_Pb;
  c = TMath::Power(Total_Energy,2)+2*Total_Energy*fPar->M_H; 
  Ex = TMath::Sqrt(a-b-c+2*TMath::Sqrt(b)*TMath::Sqrt(c)*TMath::Cos(H_Theta))-fPar->M_Pb;
  
  Histo->h_Ex->Fill(Ex);
  Histo->h_Kinematics_plot->Fill(H_Theta*180/TMath::Pi(),Total_Energy*1000);
  
}




void TNectarMonitoringProc:: Update_Histograms_Target_on(TNectarDetectors* Telescope)
{
  
  Int_t Adc=0;
  //Fill the vertical Strips Histograms::  BB8
  for(Int_t i=0; i<16; i++){
    Adc=Telescope->BB8_VStrip[i+1];
    if(Adc>0)  Histo->BB8_V_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB8_ADC_V_strip->Fill(i+1,Adc);
  }
  //Fill the horizontal Strips Histograms:: BB8
  for(Int_t i=0; i<16; i++){
    Adc=Telescope->BB8_HStrip[i+1]; 
    if(Adc>0) Histo->BB8_H_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB8_ADC_H_strip->Fill(i+1,Adc); 
  }    
  //Fill the ADC_histogram:  MSX04 Detectors
  // MSX04_1 & 2
  if(Telescope->ADC_MSX04_1>0) Histo->h_ADC_MSX04_1->Fill(Telescope->ADC_MSX04_1);
  if(Telescope->ADC_MSX04_1>0) Histo->MSXO4->Fill(1,Telescope->ADC_MSX04_1);

  if(Telescope->ADC_MSX04_2>0) Histo->h_ADC_MSX04_2->Fill(Telescope->ADC_MSX04_2);
  if(Telescope->ADC_MSX04_2>0) Histo->MSXO4->Fill(2,Telescope->ADC_MSX04_2);
  // MSX04_7 & 8
  if(Telescope->ADC_MSX04_7>0) Histo->h_ADC_MSX04_7->Fill(Telescope->ADC_MSX04_7);
  if(Telescope->ADC_MSX04_7>0) Histo->MSXO4->Fill(7,Telescope->ADC_MSX04_7);

  if(Telescope->ADC_MSX04_8>0) Histo->h_ADC_MSX04_8->Fill(Telescope->ADC_MSX04_8);
  if(Telescope->ADC_MSX04_8>0) Histo->MSXO4->Fill(8,Telescope->ADC_MSX04_8);
  // MSX04_9 & 10
  if(Telescope->ADC_MSX04_9>0) Histo->h_ADC_MSX04_9->Fill(Telescope->ADC_MSX04_9);
  if(Telescope->ADC_MSX04_9>0) Histo->MSXO4->Fill(9,Telescope->ADC_MSX04_9);
  
  if(Telescope->ADC_MSX04_10>0) Histo->h_ADC_MSX04_10->Fill(Telescope->ADC_MSX04_10);
  if(Telescope->ADC_MSX04_10>0) Histo->MSXO4->Fill(10,Telescope->ADC_MSX04_10);
  
  //Delta time histograms
  
  if(Telescope->DT_telescope>0) Histo->h_DT_telescope->Fill(Telescope->DT_telescope);
  
  //*******************
  //Histograms HR 
  //*******************
  
  Adc=0;
  //Fill the vertical Strips Histograms::  B29
  for(Int_t i=0; i<128; i++){
    Adc=Telescope->BB29_VStrip[i+1];
    if(Adc>0) Histo->BB29_V_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB29_HR_ADC_V_STRIP->Fill(i+1,Adc);
  }
  //Fill the horizontal Strips Histograms:: B29
  for(Int_t i=0; i<64; i++){
    Adc=Telescope->BB29_HStrip[i+1]; 
    if(Adc>0) Histo->BB29_H_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB29_HR_ADC_H_STRIP->Fill(i+1,Adc);
  }    
  //Delta time histograms
  
  if(Telescope->DT_BB29_V>0) Histo->h_DT_BB29_V->Fill(Telescope->DT_BB29_V);
  if(Telescope->DT_BB29_H>0) Histo->h_DT_BB29_H->Fill(Telescope->DT_BB29_H);
  
  
    //*******************
    //Histograms fission
    //*******************
    
    /* TOP */
    //Fill the vertical Strips Histograms
  for(Int_t i=0; i<16; i++){
    Adc=Telescope->BB36_top_VStrip[i+1];
    if(Adc>0) Histo->BB36_FFtop_V_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB36_FFtop_ADC_V_STRIP->Fill(i+1,Adc);

  }
  //Fill the horizontal Strips Histograms
  for(Int_t i=0; i<16; i++){
    Adc=Telescope->BB36_top_HStrip[i+1]; 
    if(Adc>0) Histo->BB36_FFtop_H_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB36_FFtop_ADC_H_STRIP->Fill(i+1,Adc);
  } 
  
    /* BOT */
    //Fill the vertical Strips Histograms
  for(Int_t i=0; i<16; i++){
    Adc=Telescope->BB36_bot_VStrip[i+1];
    if(Adc>0) Histo->BB36_FFbot_V_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB36_FFbot_ADC_V_STRIP->Fill(i+1,Adc);

  }
  //Fill the horizontal Strips Histograms
  for(Int_t i=0; i<16; i++){
    Adc=Telescope->BB36_bot_HStrip[i+1]; 
    if(Adc>0) Histo->BB36_FFbot_H_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB36_FFbot_ADC_H_STRIP->Fill(i+1,Adc);
  }
  
    /* SIDE */
    
  //Fill the vertical Strips Histograms::  B29
  for(Int_t i=0; i<128; i++){
    Adc=Telescope->BB29_side_VStrip[i+1];
    if(Adc>0) Histo->BB29_FFside_V_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB29_FFside_ADC_V_STRIP->Fill(i+1,Adc);
  }
  //Fill the horizontal Strips Histograms:: B29
  for(Int_t i=0; i<64; i++){
    Adc=Telescope->BB29_side_HStrip[i+1]; 
    if(Adc>0) Histo->BB29_FFside_H_S[i]->Fill(Adc);
    if(Adc>0) Histo->BB29_FFside_ADC_H_STRIP->Fill(i+1,Adc);
  }  
  
  
  
}


void TNectarMonitoringProc:: Update_Histograms_Target_off(TNectarDetectors* Telescope)
{
  
  //Int_t Adc=0;
  //Fill the vertical Strips Histograms::  BB8
  for(Int_t i=0; i<16; i++){
    if(Telescope->BB8_HStrip[i+1]>0)  Histo->BB8_ADC_V_strip_TO->Fill(i+1,Telescope->BB8_VStrip[i+1]);
  }
  //Fill the horizontal Strips Histograms:: BB8
  for(Int_t i=0; i<16; i++){
    if(Telescope->BB8_HStrip[i+1]>0)  Histo->BB8_ADC_H_strip->Fill(i+1,Telescope->BB8_HStrip[i+1]); 
  }    
  //Fill the ADC_histogram:  MSX04 Detectors
  // MSX04_1 & 2
  if(Telescope->ADC_MSX04_1>0) Histo->MSXO4_TO->Fill(1,Telescope->ADC_MSX04_1);
  if(Telescope->ADC_MSX04_2>0) Histo->MSXO4_TO->Fill(2,Telescope->ADC_MSX04_2);
  // MSX04_7 & 8
  if(Telescope->ADC_MSX04_7>0) Histo->MSXO4_TO->Fill(7,Telescope->ADC_MSX04_7);
  if(Telescope->ADC_MSX04_8>0) Histo->MSXO4_TO->Fill(8,Telescope->ADC_MSX04_8);
  // MSX04_9 & 10
  if(Telescope->ADC_MSX04_9>0) Histo->MSXO4_TO->Fill(9,Telescope->ADC_MSX04_9);
  if(Telescope->ADC_MSX04_10>0) Histo->MSXO4_TO->Fill(10,Telescope->ADC_MSX04_10);

  
   //Delta time histograms

  if(Telescope->DT_telescope>0) Histo->h_DT_telescope_TO->Fill(Telescope->DT_telescope);

  //*******************
  //Histograms BB29 
  //*******************
  
  //Adc=0;
  //Fill the vertical Strips Histograms::  B29
  for(Int_t i=0; i<128; i++){
    if(Telescope->BB29_VStrip[i+1]>0)  Histo->BB29_ADC_V_strip_TO->Fill(i+1,Telescope->BB29_VStrip[i+1]);
  }
  //Fill the horizontal Strips Histograms:: B29
  for(Int_t i=0; i<64; i++){
    if(Telescope->BB29_HStrip[i+1]>0)  Histo->BB29_ADC_H_strip_TO->Fill(i+1,Telescope->BB29_HStrip[i+1]);
  }    
   //Delta time histograms

  if(Telescope->DT_BB29_V>0) Histo->h_DT_BB29_V_TO->Fill(Telescope->DT_BB29_V);
  if(Telescope->DT_BB29_H>0) Histo->h_DT_BB29_H_TO->Fill(Telescope->DT_BB29_H);

}



