#include "LiAna.h"
#include "LafKernel/AlgFactory.h"
#include  <algorithm>
#include  <fstream>

DECLARE_ALGORITHM(LiAna);

    LiAna::LiAna(const std::string& name)
: AlgBase(name)
{
    OptionParser::setOption(name,"PromptEnergyLow4Li",promptELow4Li);
    OptionParser::setOption(name,"PromptEnergyHigh4Li",promptEHigh4Li);
    OptionParser::setOption(name,"DelayedEnergyLow4Li",delayedELow4Li);
    OptionParser::setOption(name,"DelayedEnergyHigh4Li",delayedEHigh4Li);
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
        for( int j=0 ; j<5 ; j++ )
        {
            lastshowermuonTrigtime[j][i].SetSec(0);
            lastshowermuonTrigtime[j][i].SetNanoSec(0);
        }

        lastAdMuonTrigtime[i].SetSec(0);
        lastAdMuonTrigtime[i].SetNanoSec(0);
        lastShowerMuonTrigtime[i].SetSec(0);
        lastShowerMuonTrigtime[i].SetNanoSec(0);
        finalTestADMuon[i]=NULL;
        adMuonTriggerTimeBuf[i].clear();
        adMuonEnergyBuf[i].clear();
        nextImuonTriggerTime[i]=0.;
    }

}

bool LiAna::initialize()
{
    EvtBuffer = dynamic_cast<PhyEventBuf*>(service("Cycler"));
    liveTimeSvc = dynamic_cast<LiveTimeSvc*>(service("LiveTimeSvc"));
    muonVeto_l = MuonVeto::instance();


    for( int i=0 ; i<5 ; i++ )
    {
        histName="time2lastshowermuon";
        histName+=i+1;
        time2lastshowermuon[i] = new TH1F(histName,"time2lastshowermuon",100000,0.,100.); 
        ntupleSvc()->attach("FILE1/lidj",time2lastshowermuon[i]);
        histName="showermuonNum";
        histName+=i+1;
        showermuonNum[i] = new TH1F(histName,"number of AD showermuon",liveTimeSvc->nBins(),liveTimeSvc->startTime().AsDouble(),liveTimeSvc->endTime().AsDouble());
        ntupleSvc()->attach("FILE1/lidj",showermuonNum[i]);
        histName="muonEnergy";
        histName+=i+1;
        muonEnergy[i] = new TH1F(histName,"energy of AD muon",10000,0,10000);
        ntupleSvc()->attach("FILE1/lidj",muonEnergy[i]);

    }

    LiTree = ntupleSvc()->bookTree("FILE2/Tree/LiTree","LiTree");
    LiTree->Branch("det",&det_l,"det_l/I");
    LiTree->Branch("energy",energy_l,"energy_l[2]/F");
    LiTree->Branch("x",x_l,"x_l[2]/F");
    LiTree->Branch("y",y_l,"y_l[2]/F");
    LiTree->Branch("z",z_l,"z_l[2]/F");
    LiTree->Branch("timeInterval",&timeInterval,"timeInterval/D");
    LiTree->Branch("promptT2Muon",promptT2Muon,"promptT2Muon[9]/D");


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


        if( AdMuonEnergy>=20. && AdMuonEnergy<1500. )
        {
            lastshowermuonTrigtime[0][CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
            lastshowermuonTrigtime[0][CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
            showermuonNum[0]->Fill(lastshowermuonTrigtime[0][CurEvent->m_det-1]);
            muonEnergy[0]->Fill(AdMuonEnergy);
            //pre AD muon
            PhyEvent* preADMuon=muonVeto_l->preMuon(CurEvent,CurEvent->m_det);
            double time2preADMuon=1000000.;
            if( preADMuon )
            {
                time2preADMuon=CurEvent->m_trigTime - preADMuon->m_trigTime;
            } 
            //next AD muon
            PhyEvent* nextADMuon=muonVeto_l->nextMuon(CurEvent,CurEvent->m_det);
            double time2nextADMuon=1000000.;
            if( nextADMuon )
            {
                time2nextADMuon=nextADMuon->m_trigTime - CurEvent->m_trigTime;
            } 
            if( time2preADMuon>=0.1 && time2nextADMuon>=0.1  )
            {
                lastshowermuonTrigtime[3][CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
                lastshowermuonTrigtime[3][CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
                showermuonNum[3]->Fill(lastshowermuonTrigtime[3][CurEvent->m_det-1]);
                muonEnergy[3]->Fill(AdMuonEnergy);
            }
            if(  time2preADMuon>=0.2 )
            {
                lastshowermuonTrigtime[4][CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
                lastshowermuonTrigtime[4][CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
                showermuonNum[4]->Fill(lastshowermuonTrigtime[4][CurEvent->m_det-1]);
                muonEnergy[4]->Fill(AdMuonEnergy);
            }

        } else if( AdMuonEnergy>=1500. && AdMuonEnergy<2500. )
        {
            lastshowermuonTrigtime[1][CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
            lastshowermuonTrigtime[1][CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
            showermuonNum[1]->Fill(lastshowermuonTrigtime[1][CurEvent->m_det-1]);
            muonEnergy[1]->Fill(AdMuonEnergy);

        }else if( AdMuonEnergy>=2500. )
        {
            lastshowermuonTrigtime[2][CurEvent->m_det-1].SetSec(CurEvent->m_trigTime.GetSec());
            lastshowermuonTrigtime[2][CurEvent->m_det-1].SetNanoSec(CurEvent->m_trigTime.GetNanoSec());
            showermuonNum[2]->Fill(lastshowermuonTrigtime[2][CurEvent->m_det-1]);
            muonEnergy[2]->Fill(AdMuonEnergy);

        }else
        {
            //continue;
        }
    }

    if( !CurEvent->isAD())
    {
        return true;
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
    for( int i=0 ; i<9 ; i++ )
    {
        delayedT2Muon[i]=time2MuonBuf[EvtGroup[1]->m_det-1][1][i];
    }
    if( !((delayedT2Muon[0]>=DelayedTrigTime2AdMuon4Li)&&(delayedT2Muon[2]>=DelayedTrigTime2IWpMuon4Li)&&(delayedT2Muon[1]>=DelayedTrigTime2AdShowerMuon4Li)&&(delayedT2Muon[3]>=DelayedTrigTime2OWpMuon4Li)) )
    {
        return true;
    }

    //fill IBD tree AND histograms
    //************************
    for(int i=0;i<9;i++) promptT2Muon[i]=10.e6;
    for(int i=0;i<8;i++) promptT2Muon[i]=time2MuonBuf[EvtGroup[0]->m_det-1][0][i];

    int detTmp=EvtGroup[0]->m_det-1;
    if( nextImuonTriggerTime[detTmp]<=EvtGroup[0]->m_trigTime.AsDouble() )
    {
        PhyEvent* nextADMuon=NULL;
        double lastAdmuonTrigTime=0.;
        if( finalTestADMuon[detTmp] && finalTestADMuon[detTmp]->m_trigTime>EvtGroup[0]->m_trigTime &&  finalTestADMuon[detTmp]->m_trigTime-EvtGroup[0]->m_trigTime<=0.2 )
        {
            nextADMuon=muonVeto_l->nextMuon(finalTestADMuon[detTmp],EvtGroup[0]->m_det);
            lastAdmuonTrigTime=finalTestADMuon[detTmp]->m_trigTime.AsDouble();
        }else
        {
            adMuonTriggerTimeBuf[detTmp].clear();
            adMuonEnergyBuf[detTmp].clear();
            double lastTrigtime=EvtGroup[0]->m_trigTime.AsDouble()-(promptT2Muon[0]>promptT2Muon[1]?promptT2Muon[1]:promptT2Muon[0]);
            adMuonTriggerTimeBuf[detTmp].push_back(lastTrigtime);
            lastAdmuonTrigTime=lastTrigtime;
            adMuonEnergyBuf[detTmp].push_back(15);
            nextADMuon=muonVeto_l->nextMuon(EvtGroup[0],EvtGroup[0]->m_det);
        }
        bool doLoop=1;
        bool findOut=0;
        while(doLoop)
        {
            if( !adMuonTriggerTimeBuf[detTmp].empty() )
            {
                if(adMuonTriggerTimeBuf[detTmp].size()!=adMuonEnergyBuf[detTmp].size()) cout<<"Error :  adMuonTriggerTimeBuf[detTmp].size()!=adMuonEnergyBuf[detTmp].size()"<<endl;
                double LastEventTriggerTime=adMuonTriggerTimeBuf[detTmp][adMuonTriggerTimeBuf[detTmp].size()-1];
                if( adMuonTriggerTimeBuf[detTmp].size()==1 && LastEventTriggerTime>EvtGroup[0]->m_trigTime.AsDouble())
                {
                    //find one
                    if( adMuonEnergyBuf[detTmp][0]>=20 && adMuonEnergyBuf[detTmp][0]<1500)
                    {
                        nextImuonTriggerTime[detTmp]=LastEventTriggerTime;
                        findOut=1;
                    }
                }
                adMuonTriggerTimeBuf[detTmp].clear();
                adMuonEnergyBuf[detTmp].clear();
            }
            doLoop=0;
            if( nextADMuon)
            {
                if( nextADMuon->m_trigTime.AsDouble()-lastAdmuonTrigTime>=0.2 )
                {
                    adMuonTriggerTimeBuf[detTmp].push_back(nextADMuon->m_trigTime.AsDouble());
                    adMuonEnergyBuf[detTmp].push_back(nextADMuon->energy());
                }
                finalTestADMuon[detTmp]=nextADMuon->GrabInstance();
                if(  nextADMuon->m_trigTime-EvtGroup[0]->m_trigTime<=0.2 )
                {
                    lastAdmuonTrigTime=nextADMuon->m_trigTime.AsDouble();
                    nextADMuon=muonVeto_l->nextMuon(nextADMuon,EvtGroup[0]->m_det);
                    doLoop=!findOut;
                }
            }
        } 
    }
    if(nextImuonTriggerTime[detTmp]!=0) promptT2Muon[8]=EvtGroup[0]->m_trigTime.AsDouble()-nextImuonTriggerTime[detTmp];

    //************************

    det_l = EvtGroup[1]->m_det;
    timeInterval=EvtGroup[1]->m_trigTime-EvtGroup[0]->m_trigTime;
    for( int i=0 ; i<(int)EvtGroup.size() ; i++ )
    {
        energy_l[i]=EvtGroup[i]->energy();
        x_l[i]=EvtGroup[i]->m_x;
        y_l[i]=EvtGroup[i]->m_y;
        z_l[i]=EvtGroup[i]->m_z;
    }
    LiTree->Fill();

    if( !((EvtGroup[0]->energy()>=3.5)&&(EvtGroup[1]->m_trigTime-EvtGroup[0]->m_trigTime)<1.e-4)) 
    {
        return true;
    }
    std::cout<<"det  : "<<EvtGroup[0]->m_det<<" trigtime  : "<<EvtGroup[0]->m_trigTime<<endl;
    for( int i=0 ; i<5 ; i++ )
    {
        time2lastshowermuon[i]->Fill(time2MuonBuf[EvtGroup[0]->m_det-1][0][i+4]);
    }

    return true;
}


bool LiAna::CalTime2Muon(PhyEvent* event)
{
    //can not be used in dump().must only be used for current event,because these last*muonTrigtimes are latest muons!!!
    time2MuonVec.clear();
    for(int i=0;i<9;i++) time2Muon[i]=10.e6;
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
    for( int i=0 ; i<5 ; i++ )
    {
        if( lastshowermuonTrigtime[i][event->m_det-1].GetSec()!=0. )
        {
            time2Muon[i+4]=event->m_trigTime - lastshowermuonTrigtime[i][event->m_det-1];
        }
    }


    for( int i=0 ; i<9 ; i++ )
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
