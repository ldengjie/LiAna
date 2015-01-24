#include  <iostream>
#include  "TH1F.h"
#include  "TTree.h"
#include  "TMath.h"
#include  <string>
#include  "MuonVeto/MuonVeto.h"
#include  "LafKernel/PhyEvent.h"
#include  <vector>
#include  "LafKernel/OptionParser.h"
#include  "LafKernel/AlgBase.h"
#include  "LafKernel/DataBuffer.h"
#include  "LafKernel/GlobalVar.h"
#include  "LafKernel/LafLog.h"
#include  "LiveTimeSvc/LiveTimeSvc.h"
#include    "LafKernel/PhyEvent/CalibReadout.h"
using namespace std;  

class LiAna : public AlgBase
{
    public:
        LiAna(const std::string& name);
        virtual ~LiAna(){}

        virtual bool initialize();
        virtual bool execute();
        virtual bool finalize();
    private:
        bool FillLi(vector<PhyEvent*> EvtGroup);
        bool FillSingle(PhyEvent* Evt);
        bool CalTime2Muon(PhyEvent* event);
        bool printEvt(PhyEvent* CurEvent);
        void dump(int i_det,bool IsFinal,PhyEvent* CurEvent);

        TString histName;
        double time2Muon[9];
        double t2muon[9];

        TTree* LiTree;
        int det;
        int    det_l;
        float  energy_l[2];
        float  x_l[2];
        float  y_l[2];
        float  z_l[2];
        double timeInterval; //us
        double promptT2Muon[9];
        double delayedT2Muon[9];

        TTimeStamp lastOwpMuonTrigtime;
        TTimeStamp lastIwpMuonTrigtime;
        TTimeStamp lastRpcMuonTrigtime;
        TTimeStamp lastAdMuonTrigtime[4];
        TTimeStamp lastShowerMuonTrigtime[4];

        TTimeStamp lastshowermuonTrigtime[5][4];
        vector<TTimeStamp> lastshowermuonTrigtimeVec[5][4];
        TH1F* time2lastshowermuon[5];
        TH1F* showermuonNum[5];
        TH1F* muonEnergy[5];

        LiveTimeSvc* liveTimeSvc;
        PhyEventBuf* EvtBuffer;
        PhyEvent* CurEvent;
        vector<PhyEvent*> AdEvtBuf[4];
        vector<double> time2MuonVec;
        vector< vector<double> > time2MuonBuf[4];
        MuonVeto* muonVeto_l;

        vector<double> adMuonTriggerTimeBuf[4];
        vector<double> adMuonEnergyBuf[4];
        PhyEvent* finalTestADMuon[4];
        double nextImuonTriggerTime[4];

        double promptELow4Li;//PromptEnergyLow
        double promptEHigh4Li;//PromptEnergyHigh
        double delayedELow4Li;//DelayedEnergyLow
        double delayedEHigh4Li;//DelayedEnergyHigh

        double LiIntervalMin;//LiIntervalMin
        double LiIntervalMax;//LiIntervalMax
        double Time2LastBufEvent;//Time2LastBufEvent
        double DelayedTrigTime2AdMuon4Li;
        double DelayedTrigTime2IWpMuon4Li;
        double DelayedTrigTime2OWpMuon4Li;
        double DelayedTrigTime2AdShowerMuon4Li;
};
