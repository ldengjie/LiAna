#include "LiAna.h"
#include "LafKernel/AlgFactory.h"
#include  <algorithm>
#include  <fstream>

DECLARE_ALGORITHM(LiAna);

LiAna::LiAna(const std::string& name)
: AlgBase(name)
{
	OptionParser::setOption(name,"EnergyLow",ELow);
	OptionParser::setOption(name,"EnergyHigh",EHigh);
	OptionParser::setOption(name,"PromptEnergyLow4Li",promptELow4Li);
	OptionParser::setOption(name,"PromptEnergyHigh4Li",promptEHigh4Li);
	OptionParser::setOption(name,"DelayedEnergyLow4Li",delayedELow4Li);
	OptionParser::setOption(name,"DelayedEnergyHigh4Li",delayedEHigh4Li);
	OptionParser::setOption(name,"AdMuonELow",AdMuonELow);
	OptionParser::setOption(name,"AdMuonEHigh",AdMuonEHigh);
	OptionParser::setOption(name,"LiIntervalMin",LiIntervalMin);
	OptionParser::setOption(name,"LiIntervalMax",LiIntervalMax);
	OptionParser::setOption(name,"Time2LastBufEvent",Time2LastBufEvent);
	OptionParser::setOption(name,"DelayedTrigTime2AdMuon4Li",DelayedTrigTime2AdMuon4Li=1.e-3);
	OptionParser::setOption(name,"DelayedTrigTime2IWpMuon4Li",DelayedTrigTime2IWpMuon4Li=6.e-4);
	OptionParser::setOption(name,"DelayedTrigTime2OWpMuon4Li",DelayedTrigTime2OWpMuon4Li=6.e-4);
	OptionParser::setOption(name,"DelayedTrigTime2AdShowerMuon4Li",DelayedTrigTime2AdShowerMuon4Li=1.e-3);

    lastOwpMuonTrigtime.SetSec(0);
    lastOwpMuonTrigtime.SetNanoSec(0);
    lastIwpMuonTrigtime.SetSec(0);
    lastIwpMuonTrigtime.SetNanoSec(0);

	for( int i=0 ; i<GlobalVar::NumADs ; i++ )
	{
		AdEvtBuf[i].clear();
        time2MuonBuf[i].clear();
        for( int j=0 ; j<6 ; j++ )
        {
            lastshowermuonTrigtime[j][i].SetSec(0);
            lastshowermuonTrigtime[j][i].SetNanoSec(0);
            lastshowermuonTrigtimeVec[j][i].clear();
        }
        
        lastAdMuonTrigtime[i].SetSec(0);
        lastAdMuonTrigtime[i].SetNanoSec(0);
        lastShowerMuonTrigtime[i].SetSec(0);
        lastShowerMuonTrigtime[i].SetNanoSec(0);
	}

	saveSingle=false;
}

bool LiAna::initialize()
{
	EvtBuffer = dynamic_cast<PhyEventBuf*>(service("Cycler"));
	liveTimeSvc = dynamic_cast<LiveTimeSvc*>(service("LiveTimeSvc"));
	muonVeto_l = MuonVeto::instance();

    for( int i=0 ; i<7 ; i++ )
    {
        ADMuonNum[i]=0;
    }
    
    for( int i=0 ; i<6 ; i++ )
    {
        histName="time2lastshowermuon";
        histName+=i+1;
        time2lastshowermuon[i] = new TH1F(histName,"time2lastshowermuon",100000,0.,100.); 
	    ntupleSvc()->attach("FILE1/lidj",time2lastshowermuon[i]);
        histName="showermuonNum";
        histName+=i+1;
	    showermuonNum[i] = new TH1F(histName,"number of AD showermuon",liveTimeSvc->nBins(),liveTimeSvc->startTime().AsDouble(),liveTimeSvc->endTime().AsDouble());
	    ntupleSvc()->attach("FILE1/lidj",showermuonNum[i]);
        histName="time2lastshowermuon";
        histName+=i+1;
        histName+="4Li";
        time2lastshowermuon4Li[i] = new TH1F(histName,"time2lastshowermuon Li",100000,0.,100.); 
	    ntupleSvc()->attach("FILE1/lidj",time2lastshowermuon4Li[i]);

    }
    time2Allmuon=new TH1F("time2Allmuon","time2Allmuon",100000,0.,100.);
	ntupleSvc()->attach("FILE1/lidj",time2Allmuon);
    time2Allmuon4Li=new TH1F("time2Allmuon4Li","time2Allmuon4Li",100000,0.,100.);
	ntupleSvc()->attach("FILE1/lidj",time2Allmuon4Li);
	if( saveSingle )
	{
		SingleTree = ntupleSvc()->bookTree("FILE2/Tree/SingleTree","SingleTree");
		SingleTree->Branch("energy",&energy,"energy/F");
		SingleTree->Branch("x",&x,"x/F");
		SingleTree->Branch("y",&y,"y/F");
		SingleTree->Branch("z",&z,"z/F");
		SingleTree->Branch("det",&det,"det/I");
		//SingleTree->Branch("t2muon",t2muon,"t2muon[10]/F");
		SingleTree->Branch("t2lastshowermuon",&t2lastshowermuon,"t2lastshowermuon/D");
	}

	LiTree = ntupleSvc()->bookTree("FILE2/Tree/LiTree","LiTree");
	LiTree->Branch("det",&det_l,"det_l/I");
	LiTree->Branch("energy",energy_l,"energy_l[2]/F");
	LiTree->Branch("x",x_l,"x_l[2]/F");
	LiTree->Branch("y",y_l,"y_l[2]/F");
	LiTree->Branch("z",z_l,"z_l[2]/F");
	LiTree->Branch("timeInterval",&timeInterval,"timeInterval/D");
	//LiTree->Branch("promptT2Muon",promptT2Muon,"promptT2Muon[10]/F");
	LiTree->Branch("t2lastshowermuon",&t2lastshowermuon,"t2lastshowermuon/D");
	LiTree->Branch("t2lastshowermuonWithMR",&t2lastshowermuonWithMR,"t2lastshowermuonWithMR/D");

	//MutilTree = ntupleSvc()->bookTree("","");
	//MutilTree->Branch("",&,"");


	std::cout<<"LiAna initializing..."<<endl;
	return true;
}

bool LiAna::execute()
{
	
	PhyEvent* CurEvent= EvtBuffer->curEvt();
	
	if( !CurEvent->m_good )
	{
		return true;
	}
    if( CurEvent->isAdMuon() )
    {
        ADMuonNum[0]++;
    }
	if( CurEvent->isMuon())
	{
        //need to delete this function "updateVetoWindow",in "NtupleAna/Algorithms/src/MuonTagAlg.cc" ,it has been executed once.If you don't use the project"MuonTagAlg",you should not delete this line.
		//muonVeto_l->updateVetoWindow(CurEvent);

        if( CurEvent->isWpMuon() && CurEvent->m_det==5 )
        {
            lastIwpMuonTrigtime.SetSec(CurEvent->m_trigTime.GetSec());
            lastIwpMuonTrigtime.SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
        } else if(CurEvent->isWpMuon() && CurEvent->m_det==6)
        {
            lastOwpMuonTrigtime.SetSec(CurEvent->m_trigTime.GetSec());
            lastOwpMuonTrigtime.SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
        } else if( CurEvent->isAdMuon() )
        {
            lastAdMuonTrigtime[CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
            lastAdMuonTrigtime[CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
        }else if(CurEvent->isShowerMuon())
        {
            lastShowerMuonTrigtime[CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
            lastShowerMuonTrigtime[CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
        }else
        {
            // 
        }
		if( !(CurEvent->isAD()) )
		{
			return true;
		}

        double AdMuonEnergy=0.;//CurEvent->energy();
        AdMuonEnergy=CurEvent->energy();

		if( AdMuonEnergy>=20. && AdMuonEnergy<500. )
		{
            lastshowermuonTrigtimeVec[0][CurEvent->m_det-1].push_back(CurEvent->m_trigTime);
            ADMuonNum[1]++;

		} else if( AdMuonEnergy>=500. && AdMuonEnergy<1500. )
		{
            lastshowermuonTrigtimeVec[1][CurEvent->m_det-1].push_back(CurEvent->m_trigTime);
            ADMuonNum[2]++;
		} else if( AdMuonEnergy>=1500. && AdMuonEnergy<2500. )
		{
            lastshowermuonTrigtimeVec[2][CurEvent->m_det-1].push_back(CurEvent->m_trigTime);
            ADMuonNum[3]++;
		} else if( AdMuonEnergy>=2500. && AdMuonEnergy<3500. )
		{
            lastshowermuonTrigtimeVec[3][CurEvent->m_det-1].push_back(CurEvent->m_trigTime);
            ADMuonNum[4]++;
		} else if( AdMuonEnergy>=3500. && AdMuonEnergy<4500. )
		{
            lastshowermuonTrigtimeVec[4][CurEvent->m_det-1].push_back(CurEvent->m_trigTime);
            ADMuonNum[5]++;
			
		}else if( AdMuonEnergy>=4500. )
		{
            lastshowermuonTrigtimeVec[5][CurEvent->m_det-1].push_back(CurEvent->m_trigTime);
            ADMuonNum[6]++;
			
		}else
		{
			//continue;
		}
	}

	if( !CurEvent->isAD())
	{
		return true;
	}

    bool isNeulike=false;
    if( CurEvent->energy()>1.8 && CurEvent->energy()<12.0 )
    {
        isNeulike=true;
    } 
    /*
    else if( CurEvent->energy()>1.8 && CurEvent->energy()<2.8 )
    {
            isNeulike=true;
    }else
    {
        //
    }
    */
    if(isNeulike) 
    {

        //for( int i=0 ; i<6 ; i++ )
        for( int i=0 ; i<6 ; i++ )
        {
            if( !lastshowermuonTrigtimeVec[i][CurEvent->m_det-1].empty() )
            {
                vector<TTimeStamp>::iterator it=lastshowermuonTrigtimeVec[i][CurEvent->m_det-1].end();
                it--;
                for(; it>=lastshowermuonTrigtimeVec[i][CurEvent->m_det-1].begin();it--) 
                {
                    if( (CurEvent->m_trigTime-*it)>1.e-5 &&(CurEvent->m_trigTime-*it)<2.e-4)
                    {
                        lastshowermuonTrigtime[i][CurEvent->m_det-1].SetSec((*it).GetSec());
                        lastshowermuonTrigtime[i][CurEvent->m_det-1].SetNanoSec((*it).GetNanoSec());
			            showermuonNum[i]->Fill(lastshowermuonTrigtime[i][CurEvent->m_det-1]);
                        it=lastshowermuonTrigtimeVec[i][CurEvent->m_det-1].erase(it);
                    } else if( (CurEvent->m_trigTime-*it)>=2.e-4 )
                    {
                        it++;
                        lastshowermuonTrigtimeVec[i][CurEvent->m_det-1].erase(lastshowermuonTrigtimeVec[i][CurEvent->m_det-1].begin(),it); 
                        break;
                    } else
                    {
                        //
                    }
                }
            }
        }
        
    }


	int vetotag=muonVeto_l->veto(CurEvent);

	if( vetotag==1 )
	{
		return true;
	} 

	if( !(CurEvent->energy()>0.7) )
	{
		return true;
	}

    dump(CurEvent->m_det,0,CurEvent);

	AdEvtBuf[CurEvent->m_det-1].push_back(CurEvent->GrabInstance());
    CalTime2Muon(CurEvent);
    time2MuonBuf[CurEvent->m_det-1].push_back(time2MuonVec);

	return true;

}

bool LiAna::finalize()
{
    PhyEvent* CurEvent=0;
    for( int i=1 ; i<5 ; i++ )
    {
        dump(i,1,CurEvent);
    }
    
	std::cout<<"LiAna finalize..."<<endl;
    for( int i=0 ; i<7 ; i++ )
    {
        std::cout<<"ADmuon number "<<i <<" : "<<ADMuonNum[i]<<endl;
    }
    
	return true;
}

bool LiAna::FillLi(vector<PhyEvent*> EvtGroup)
{
    if( time2MuonBuf[EvtGroup[0]->m_det-1].size()!=2 )
    {
        std::cout<<"time2MuonBuf size is wrong in FillLi(),size is "<<time2MuonBuf[EvtGroup[0]->m_det-1].size()<<endl;
    }

	if(!( (EvtGroup[0]->energy()<promptEHigh4Li)&&(EvtGroup[1]->energy()<delayedEHigh4Li)&&(EvtGroup[1]->energy()>delayedELow4Li)&&((EvtGroup[1]->m_trigTime-EvtGroup[0]->m_trigTime)<LiIntervalMax)&&((EvtGroup[1]->m_trigTime-EvtGroup[0]->m_trigTime)>=LiIntervalMin) ))
	{
		return true;
	}
	for( int i=0 ; i<10 ; i++ )
	{
		delayedT2Muon[i]=time2MuonBuf[EvtGroup[1]->m_det-1][1][i];	
	}
	if( !((delayedT2Muon[0]>=DelayedTrigTime2AdMuon4Li)&&(delayedT2Muon[2]>=DelayedTrigTime2IWpMuon4Li)&&(delayedT2Muon[1]>=DelayedTrigTime2AdShowerMuon4Li)&&(delayedT2Muon[3]>=DelayedTrigTime2OWpMuon4Li)) )
	{
	    return true;
	}
    /*
    float Dis=(EvtGroup[1]->m_x-EvtGroup[0]->m_x)*(EvtGroup[1]->m_x-EvtGroup[0]->m_x)+(EvtGroup[1]->m_y-EvtGroup[0]->m_y)*(EvtGroup[1]->m_y-EvtGroup[0]->m_y)+(EvtGroup[1]->m_z-EvtGroup[0]->m_z)*(EvtGroup[1]->m_z-EvtGroup[0]->m_z);
    if(!(Dis<1.e6))
    {
        return true;
    }
    */
    /*
    if( time2MuonBuf[EvtGroup[1]->m_det-1][0][6]>2. )
    {
        return true;
    }

    */
	//fill IBD tree AND histograms
	det_l = EvtGroup[1]->m_det;
	timeInterval=(EvtGroup[1]->m_trigTime-EvtGroup[0]->m_trigTime)*1.e6;
	for( int i=0 ; i<(int)EvtGroup.size() ; i++ )
	{
		energy_l[i]=EvtGroup[i]->energy();
		x_l[i]=EvtGroup[i]->m_x;
		y_l[i]=EvtGroup[i]->m_y;
		z_l[i]=EvtGroup[i]->m_z;
	}
    /*
	for( int i=0 ; i<10 ; i++ )
	{
		promptT2Muon[i]=time2MuonBuf[EvtGroup[0]->m_det-1][0][i];	
	}
	*/
    t2lastshowermuon=time2MuonBuf[EvtGroup[0]->m_det-1][0][1];
    t2lastshowermuonWithMR=0.;
    double t1=time2MuonBuf[EvtGroup[0]->m_det-1][0][7];
    double t2=time2MuonBuf[EvtGroup[0]->m_det-1][0][8];
    double t3=time2MuonBuf[EvtGroup[0]->m_det-1][0][9];
    double t4=t1>t2?t2:t1;
    double t5=t4>t3?t3:t4;
    t2lastshowermuonWithMR=t5;
	LiTree->Fill();



    if( !((EvtGroup[0]->energy()>=3.5)&&(EvtGroup[1]->m_trigTime-EvtGroup[0]->m_trigTime)<1.e-4)) 
    {
        return true;
    }
    std::cout<<"det  : "<<EvtGroup[0]->m_det<<" trigtime  : "<<EvtGroup[0]->m_trigTime<<endl;
    for( int i=0 ; i<6 ; i++ )
    {
        time2lastshowermuon4Li[i]->Fill(time2MuonBuf[EvtGroup[0]->m_det-1][0][i+4]);
    }
    time2Allmuon4Li->Fill(time2MuonBuf[EvtGroup[0]->m_det-1][0][1]);

	return true;
}


bool LiAna::FillSingle(PhyEvent* Evt)
{
    if( time2MuonBuf[Evt->m_det-1].size()!=1 )
    {
        std::cout<<"wrong !! time2MuonBuf["<<Evt->m_det-1<<"].size()="<<time2MuonBuf[Evt->m_det-1].size() <<endl;
    }
	energy=Evt->energy();
	//if( energy>EHigh || energy<ELow)
	if( energy>12. || energy<6.)
	{
		return true;
	}

    /*
	if( ((time2MuonBuf[Evt->m_det-1][0][0]>=DelayedTrigTime2AdMuon)&&(time2MuonBuf[Evt->m_det-1][0][2]>=DelayedTrigTime2IWpMuon)&&(time2MuonBuf[Evt->m_det-1][0][1]>=DelayedTrigTime2AdShowerMuon)&&(time2MuonBuf[Evt->m_det-1][0][3]>=DelayedTrigTime2OWpMuon)) )
	{
	}
    */
    for( int i=0 ; i<6 ; i++ )
    {
        time2lastshowermuon[i]->Fill(time2MuonBuf[Evt->m_det-1][0][i+4]);
    }
    time2Allmuon->Fill(time2MuonBuf[Evt->m_det-1][0][1]);
    if( ((Evt->m_quadrant>0.5)&&(Evt->m_maxQ>0.03)) )
    {
        return true;
    }
	if( saveSingle )
	{
		det=Evt->m_det;
		x=Evt->m_x;
		y=Evt->m_y;
		z=Evt->m_z;
		//CalTime2Muon(Evt);
		//for(int i=0;i<10;i++) t2muon[i]=time2Muon[i];
        t2lastshowermuon=time2MuonBuf[Evt->m_det-1][0][1];
		SingleTree->Fill();
	}
	return true;
}

bool LiAna::CalTime2Muon(PhyEvent* event)
{
    //can not be used in dump().must only be used for current event,because these last*muonTrigtimes are latest muons!!!
    time2MuonVec.clear();
	for(int i=0;i<10;i++) time2Muon[i]=10.e6;
	//pre ADMuon
    if( lastAdMuonTrigtime[event->m_det-1].GetSec()!=0. )
    {
        time2Muon[0]=event->m_trigTime-lastAdMuonTrigtime[event->m_det-1];
    }
	//pre ADshower muon
    if( lastShowerMuonTrigtime[event->m_det-1].GetSec()!=0. )
    {
        time2Muon[1]=event->m_trigTime - lastShowerMuonTrigtime[event->m_det-1];
    }
	//pre Iwp muon
    if( lastIwpMuonTrigtime.GetSec()!=0. )
    {
        time2Muon[2]=event->m_trigTime - lastIwpMuonTrigtime;
    }
	//pre Owp muon
    if( lastOwpMuonTrigtime.GetSec()!=0. )
    {
        time2Muon[3]=event->m_trigTime - lastOwpMuonTrigtime;
    }
    //for He8/Li9
    for( int i=0 ; i<6 ; i++ )
    {
        if( lastshowermuonTrigtime[i][event->m_det-1].GetSec()!=0. )
        {
            time2Muon[i+4]=event->m_trigTime - lastshowermuonTrigtime[i][event->m_det-1];
        }
    }
    

    for( int i=0 ; i<10 ; i++ )
    {
        time2MuonVec.push_back(time2Muon[i]);
    }
	return true;
}


bool LiAna::printEvt(PhyEvent* CurEvent)
{
	std::cout<<"===> info <==="<<endl;
	std::cout<<"entry = "<<CurEvent->m_entry<<endl;
	std::cout<<"fileNum = "<<CurEvent->m_fileNum<<endl;
	std::cout<<"localentry = "<<CurEvent->m_localEntry<<endl;
	std::cout<<"trigtime = "<<CurEvent->m_trigTime<<endl;
	std::cout<<"trigType = "<<CurEvent->m_trigType<<endl;
	std::cout<<"det = "<<CurEvent->m_det<<endl;
	std::cout<<"energy = "<<CurEvent->m_energy<<endl;
	std::cout<<"x = "<<CurEvent->m_x<<endl;
	std::cout<<"y = "<<CurEvent->m_y<<endl;
	std::cout<<"z = "<<CurEvent->m_z<<endl;
	std::cout<<"rawEvis = "<<CurEvent->m_rawEvis<<endl;
	std::cout<<"energy() = "<<CurEvent->energy()<<endl;
	std::cout<<"nPmt = "<<CurEvent->m_nPmt<<endl;
	std::cout<<"flasherTag = "<<CurEvent->m_flasherTag<<endl;
	std::cout<<"forceTrigTag = "<<CurEvent->m_forceTrigTag<<endl;
	std::cout<<"crossTrigTag = "<<CurEvent->m_crossTrigTag<<endl;
	std::cout<<"rpcNoiseTag = "<<CurEvent->m_rpcNoiseTag<<endl;
	std::cout<<"adLowEnergyTag = "<<CurEvent->m_adLowEnergyTag<<endl;
	std::cout<<"goodevent = "<<CurEvent->m_good<<endl;
	std::cout<<"muontag = "<<CurEvent->m_muonTag<<endl;
	return true;
}

void LiAna::dump(int i_det,bool IsFinal,PhyEvent* CurEvent)
{
	if( !AdEvtBuf[i_det-1].empty() )
	{
		PhyEvent* LastEvent=AdEvtBuf[i_det-1][AdEvtBuf[i_det-1].size()-1];
		if( IsFinal || (CurEvent->m_trigTime-LastEvent->m_trigTime)>Time2LastBufEvent )
		{
			if( AdEvtBuf[i_det-1].size()==1 )
			{
				//FillSingle(AdEvtBuf[i_det-1][0]);
			} else if(AdEvtBuf[i_det-1].size()==2)
			{
				FillLi(AdEvtBuf[i_det-1]);
			}else
            {
                //
            }
			for( int i=0 ; i<(int)AdEvtBuf[i_det-1].size();i++)
			{
				AdEvtBuf[i_det-1][i]->ReleaseInstance();
			}
			
			AdEvtBuf[i_det-1].clear(); 
            time2MuonBuf[i_det-1].clear();
		}
	} 
    
}
