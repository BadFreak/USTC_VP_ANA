#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <set>
#include <vector>

// DetectorID: 板上 FEE 5–8；否则从 CellID 取 fee 1–4
// CellID = cry*100000 + feeid*10000 + ...
static int feeKeyFromPacket(int detectorId, const std::vector<int>* cellIds)
{
    if (detectorId >= 5 && detectorId <= 8)
        return detectorId;
    if (cellIds && !cellIds->empty()) {
        const int cid = (*cellIds)[0];
        if (cid > 0)
            return (cid % 100000) / 10000; // 1–4
    }
    return -1;
}

void pack_calo_by_trigger(const char* infile  = "result_Calo_cosmic_20260403_1422.root",
                                 const char* outfile = "packed_calo.root")
{
    TFile* fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << infile << std::endl;
        return;
    }

    TTree* caloTree = dynamic_cast<TTree*>(fin->Get("caloTree"));
    if (!caloTree) {
        std::cerr << "Error: cannot find caloTree in " << infile << std::endl;
        fin->Close();
        return;
    }

    int TriggerID = 0;
    int DetectorID = 0;
    std::vector<int>* inCellID   = nullptr;
    std::vector<int>* inCellADC  = nullptr;
    std::vector<int>* inCellPLAT = nullptr;

    caloTree->SetBranchAddress("TriggerID", &TriggerID);
    caloTree->SetBranchAddress("CellID",    &inCellID);
    caloTree->SetBranchAddress("CellADC",   &inCellADC);
    caloTree->SetBranchAddress("CellPLAT",  &inCellPLAT);
    if (caloTree->GetBranch("DetectorID"))
        caloTree->SetBranchAddress("DetectorID", &DetectorID);

    TFile* fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << outfile << std::endl;
        fin->Close();
        return;
    }

    TTree* outTree = new TTree("caloTree",
        "packed events: 4 packets same TriggerID, screened (4 distinct FEEs)");

    int EventID = 0;
    int outTriggerID = 0;
    std::vector<int> outCellID;
    std::vector<int> outCellADC;
    std::vector<int> outCellPLAT;

    outTree->Branch("EventID",   &EventID);
    outTree->Branch("TriggerID", &outTriggerID);
    outTree->Branch("CellID",    &outCellID);
    outTree->Branch("CellADC",   &outCellADC);
    outTree->Branch("CellPLAT",  &outCellPLAT);

    Long64_t nentries = caloTree->GetEntries();
    std::cout << "Input entries = " << nentries << std::endl;

    int currentTrigger = -1;
    int packetCount = 0;
    std::vector<int> bufFee;
    std::vector<int> bufCellID;
    std::vector<int> bufCellADC;
    std::vector<int> bufCellPLAT;

    int nPacked = 0;
    int nDroppedIncomplete = 0;
    int nDroppedScreen = 0;

    auto resetBuf = [&]() {
        packetCount = 0;
        currentTrigger = -1;
        bufFee.clear();
        bufCellID.clear();
        bufCellADC.clear();
        bufCellPLAT.clear();
    };

    auto startGroup = [&]() {
        currentTrigger = TriggerID;
        packetCount = 1;
        bufFee = {feeKeyFromPacket(DetectorID, inCellID)};
        bufCellID   = *inCellID;
        bufCellADC  = *inCellADC;
        bufCellPLAT = *inCellPLAT;
    };

    auto appendPacket = [&]() {
        ++packetCount;
        bufFee.push_back(feeKeyFromPacket(DetectorID, inCellID));
        bufCellID.insert(bufCellID.end(), inCellID->begin(), inCellID->end());
        bufCellADC.insert(bufCellADC.end(), inCellADC->begin(), inCellADC->end());
        bufCellPLAT.insert(bufCellPLAT.end(), inCellPLAT->begin(), inCellPLAT->end());
    };

    // pack 后筛查：4 个不同 FEE（5–8 或 1–4）
    auto screenGroup = [&]() -> bool {
        if (packetCount != 4)
            return false;
        std::set<int> fees;
        for (int f : bufFee) {
            if (f < 0)
                return false;
            fees.insert(f);
        }
        return fees.size() == 4;
    };

    auto tryFillOrDrop = [&]() {
        if (screenGroup()) {
            EventID = nPacked;
            outTriggerID = currentTrigger;
            outCellID = bufCellID;
            outCellADC = bufCellADC;
            outCellPLAT = bufCellPLAT;
            outTree->Fill();
            ++nPacked;
        } else {
            ++nDroppedScreen;
        }
        resetBuf();
    };

    auto dropIncomplete = [&]() {
        if (packetCount > 0 && packetCount < 4)
            ++nDroppedIncomplete;
        resetBuf();
    };

    for (Long64_t i = 0; i < nentries; ++i) {
        caloTree->GetEntry(i);

        if (!inCellID || !inCellADC || !inCellPLAT) {
            std::cerr << "Warning: null vector pointer at entry " << i << std::endl;
            continue;
        }

        if (inCellID->size() != inCellADC->size() ||
            inCellID->size() != inCellPLAT->size()) {
            std::cerr << "Warning: size mismatch at entry " << i
                      << "  CellID="   << inCellID->size()
                      << "  CellADC="  << inCellADC->size()
                      << "  CellPLAT=" << inCellPLAT->size()
                      << std::endl;
        }

        if (packetCount == 0) {
            startGroup();
        } else if (TriggerID == currentTrigger) {
            appendPacket();
        } else {
            dropIncomplete();
            startGroup();
        }

        if (packetCount == 4)
            tryFillOrDrop();
    }

    dropIncomplete();

    outTree->Write();
    fout->Close();
    fin->Close();

    std::cout << "Packed events          = " << nPacked << std::endl;
    std::cout << "Dropped incomplete     = " << nDroppedIncomplete
              << "  (<4 packets)" << std::endl;
    std::cout << "Dropped by screen      = " << nDroppedScreen
              << "  (not 4 distinct FEEs)" << std::endl;
    std::cout << "Output written to " << outfile << std::endl;
}
