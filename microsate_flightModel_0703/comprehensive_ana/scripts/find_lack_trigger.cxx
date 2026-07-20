#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <string>

static std::string out_dir_of(const char* rootfile)
{
    std::string path = rootfile;
    const auto slash = path.rfind('/');
    if (slash != std::string::npos)
        return path.substr(0, slash + 1);
    return "";
}

static std::string out_path_for(const char* rootfile, const char* tree_tag)
{
    return out_dir_of(rootfile) + std::string(tree_tag) + "_lack_trigger.txt";
}

static bool check_tree_gaps(TTree* tree, const char* tree_name, const char* rootfile,
                            const std::string& out_path)
{
    std::ofstream ofs(out_path);
    if (!ofs) {
        std::cerr << "find_lack_trigger: cannot write " << out_path << std::endl;
        return false;
    }

    ofs << "# find_lack_trigger\n";
    ofs << "# rootfile " << rootfile << "\n";
    ofs << "# tree " << tree_name << "\n";

    if (!tree) {
        ofs << "# " << tree_name << ": (not found)\n";
        ofs.close();
        std::cout << "Saved: " << out_path << std::endl;
        return true;
    }

    int trigger_id = 0;
    tree->SetBranchAddress("TriggerID", &trigger_id);
    const Long64_t nent = tree->GetEntries();
    ofs << "# entries=" << nent << "\n";
    if (nent < 1) {
        ofs << "# empty\n";
        ofs.close();
        std::cout << "Saved: " << out_path << std::endl;
        return true;
    }

    int pre_trigger_id = 0;
    int ngap = 0;
    for (Long64_t i = 0; i < nent; ++i) {
        tree->GetEntry(i);
        if (i == 0) {
            ofs << "first TriggerID: " << trigger_id << "\n";
            if (trigger_id == 0)
                ofs << "OKAY: first TriggerID is 0\n";
            else
                ofs << "ERROR: first TriggerID is not 0\n";
        } else if (trigger_id != pre_trigger_id + 1
                   && !(pre_trigger_id == 4095 && trigger_id == 0)) {
            ofs << "gap Entry=" << i << " pre=" << pre_trigger_id
                << " cur=" << trigger_id << " expected=" << (pre_trigger_id + 1) << "\n";
            ++ngap;
        }
        pre_trigger_id = trigger_id;
    }
    ofs << "# gap_count=" << ngap << "\n";
    ofs.close();
    std::cout << "Saved: " << out_path << std::endl;
    return true;
}

// resultfile: csiTree；packfile: pack 后 caloTree（如 pack.root）
void find_lack_trigger(const char* resultfile, const char* packfile = nullptr)
{
    if (!resultfile || !resultfile[0]) {
        std::cerr << "find_lack_trigger: empty resultfile" << std::endl;
        return;
    }

    TFile* fin_result = TFile::Open(resultfile, "READ");
    if (!fin_result || fin_result->IsZombie()) {
        std::cerr << "find_lack_trigger: cannot open " << resultfile << std::endl;
        return;
    }

    const std::string out_csi = out_path_for(resultfile, "csi");
    check_tree_gaps(dynamic_cast<TTree*>(fin_result->Get("csiTree")), "csiTree", resultfile,
                    out_csi);
    fin_result->Close();

    if (!packfile || !packfile[0]) {
        std::cout << "find_lack_trigger: no packfile, skip caloTree" << std::endl;
        return;
    }

    TFile* fin_pack = TFile::Open(packfile, "READ");
    if (!fin_pack || fin_pack->IsZombie()) {
        std::cerr << "find_lack_trigger: cannot open " << packfile << std::endl;
        return;
    }

    const std::string out_calo = out_path_for(packfile, "calo");
    check_tree_gaps(dynamic_cast<TTree*>(fin_pack->Get("caloTree")), "caloTree", packfile,
                    out_calo);
    fin_pack->Close();
}
