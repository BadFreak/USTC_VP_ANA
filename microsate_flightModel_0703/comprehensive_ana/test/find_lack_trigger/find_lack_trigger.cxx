void find_lack_trigger(){
    //TString fileDir = "/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/microsate_flightModel_260525/comprehensive_ana/result_microsat/CsI_Calo_260609_1638/cosmic";
    //TString fileName = "result_CsI_Calo_260609_1638.root";
    //TString filePath = fileDir + "/" + fileName;
	TString filePath = "/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/microsate_flightModel_260525/comprehensive_ana/result_microsat/CsI_Calo_260609_1638/cosmic/result_CsI_Calo_260609_1638.root";
    TFile *f = new TFile(filePath,"READ");
    if (!f->IsOpen()) {
        std::cout << "Failed to open file!" << std::endl;
        return;
    };
    TTree *t = (TTree*)f->Get("csiTree");
    if (!t) {
        std::cout << "Failed to get tree!" << std::endl;
        return;
    }

    int triggerID;
    t->SetBranchAddress("TriggerID", &triggerID);

    int preTriggerID = 0;
    for (int i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        // std::cout << "i:" << i << " triggerID: " << triggerID << std::endl;
        if (i==0) {
            std::cout << "First triggerID: " << triggerID << std::endl;
            if (triggerID == 0) {
                std::cout << "OKAY: NEXT" << std::endl;
            } else {
                std::cout << "ERROR: NEXT" << std::endl;
            } 
        } else {
            if (triggerID != preTriggerID + 1) {
                std::cout << "i:" << i << " pretriggerID: " << preTriggerID << " triggerID: " << triggerID << std::endl;

                // std::cout << "i:" << i << " triggerID: " << triggerID << std::endl;
            }
        }
        // Process each entry
        preTriggerID = triggerID;
    }

}
