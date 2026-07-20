#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <iostream>

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

    // ===== 输入分支 =====
    int TriggerID = 0;
    std::vector<int>* inCellID   = nullptr;
    std::vector<int>* inCellADC  = nullptr;
    std::vector<int>* inCellPLAT = nullptr;

    caloTree->SetBranchAddress("TriggerID", &TriggerID);
    caloTree->SetBranchAddress("CellID",    &inCellID);
    caloTree->SetBranchAddress("CellADC",   &inCellADC);
    caloTree->SetBranchAddress("CellPLAT",  &inCellPLAT);

    // ===== 输出文件 =====
    TFile* fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << outfile << std::endl;
        fin->Close();
        return;
    }

    TTree* outTree = new TTree("caloTree", "packed events from 4 consecutive packets with same TriggerID");

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

    // ===== 当前正在收集的一组 =====
    int currentTrigger = -1;
    int packetCount = 0;

    std::vector<int> bufCellID;
    std::vector<int> bufCellADC;
    std::vector<int> bufCellPLAT;

    int nPacked = 0;
    int nDroppedGroups = 0;

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

        // 第一个包，直接开始
        if (packetCount == 0) {
            currentTrigger = TriggerID;
            packetCount = 1;

            bufCellID   = *inCellID;
            bufCellADC  = *inCellADC;
            bufCellPLAT = *inCellPLAT;
        }
        else {
            if (TriggerID == currentTrigger) {
                // 同一个 TriggerID，继续累积
                ++packetCount;

                bufCellID.insert(bufCellID.end(), inCellID->begin(), inCellID->end());
                bufCellADC.insert(bufCellADC.end(), inCellADC->begin(), inCellADC->end());
                bufCellPLAT.insert(bufCellPLAT.end(), inCellPLAT->begin(), inCellPLAT->end());
            } else {
                // TriggerID 变化：前面没凑满4个就丢掉，从当前这个重新开始
                if (packetCount < 4) ++nDroppedGroups;

                currentTrigger = TriggerID;
                packetCount = 1;

                bufCellID   = *inCellID;
                bufCellADC  = *inCellADC;
                bufCellPLAT = *inCellPLAT;
            }
        }

        // 凑满4个包，输出一个event
        if (packetCount == 4) {
            EventID = nPacked;
            outTriggerID = currentTrigger;
            outCellID = bufCellID;
            outCellADC = bufCellADC;
            outCellPLAT = bufCellPLAT;

            outTree->Fill();
            ++nPacked;

            // 清空，准备收下一个event
            packetCount = 0;
            currentTrigger = -1;
            bufCellID.clear();
            bufCellADC.clear();
            bufCellPLAT.clear();
        }
    }

    // 文件结束时如果还剩 1/2/3 个未满4个，也丢掉
    if (packetCount > 0 && packetCount < 4) {
        ++nDroppedGroups;
    }

    outTree->Write();
    fout->Close();
    fin->Close();

    std::cout << "Packed events   = " << nPacked << std::endl;
    std::cout << "Dropped groups  = " << nDroppedGroups << std::endl;
    std::cout << "Output written to " << outfile << std::endl;
}
