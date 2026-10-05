#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// 与 vdecode ComReader 当前时基一致：2024-10-22 13:02:27
static const int kTimeBaseY = 2024, kTimeBaseMo = 10, kTimeBaseD = 22;
static const int kTimeBaseH = 13, kTimeBaseMi = 2, kTimeBaseS = 27;

static int daysInYear(int y)
{
    return ((y % 4 == 0 && y % 100 != 0) || (y % 400 == 0)) ? 366 : 365;
}

static double secFromYearStart(int y, int mo, int d, int H, int M, int S)
{
    static const int cum[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int doy = cum[mo - 1] + (d - 1);
    if (mo > 2 && daysInYear(y) == 366)
        doy++;
    return doy * 86400. + H * 3600. + M * 60. + S;
}

// 时间码（秒）-> 北京时间，含毫秒
static std::string timeCodeToBeijing(double sec)
{
    double t = secFromYearStart(kTimeBaseY, kTimeBaseMo, kTimeBaseD, kTimeBaseH, kTimeBaseMi, kTimeBaseS) + sec;
    int y = kTimeBaseY;
    while (t >= daysInYear(y) * 86400.) {
        t -= daysInYear(y) * 86400.;
        y++;
    }
    while (t < 0.) {
        y--;
        t += daysInYear(y) * 86400.;
    }
    int doy = static_cast<int>(t / 86400.);
    double sod = t - doy * 86400.;
    int H = static_cast<int>(sod / 3600.);
    sod -= H * 3600.;
    int M = static_cast<int>(sod / 60.);
    double sec_left = sod - M * 60.;
    int S = static_cast<int>(sec_left);
    int ms = static_cast<int>(std::llround((sec_left - S) * 1000.0));
    if (ms >= 1000) {
        ms -= 1000;
        S++;
        if (S >= 60) {
            S -= 60;
            M++;
            if (M >= 60) {
                M -= 60;
                H++;
            }
        }
    }
    static const int cum[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int mo = 1;
    for (int m = 12; m >= 1; --m) {
        int start = cum[m - 1];
        if (m > 2 && daysInYear(y) == 366)
            start++;
        if (doy >= start) {
            mo = m;
            break;
        }
    }
    int d = doy - cum[mo - 1] + 1;
    if (mo > 2 && daysInYear(y) == 366)
        d--;
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d.%03d", y, mo, d, H, M, S, ms);
    return std::string(buf);
}

static std::string fmtTimeCode(bool has_tc, double tc)
{
    if (!has_tc)
        return "(无TimeCode)";
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << tc;
    return os.str();
}

static std::string fmtBeijing(bool has_tc, double tc)
{
    if (!has_tc)
        return "(无TimeCode)";
    return timeCodeToBeijing(tc);
}

struct EvPt {
    Long64_t entry = 0;
    int trigger = 0;
    double tc = 0.;
    bool has_tc = false;
};

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

static int feeKeyFromPacket(int detectorId, const std::vector<int>* cellIds)
{
    if (detectorId >= 5 && detectorId <= 8)
        return detectorId;
    if (cellIds && !cellIds->empty()) {
        const int cid = (*cellIds)[0];
        if (cid > 0)
            return (cid % 100000) / 10000;
    }
    return -1;
}

// pack.root 不存 TimeCode：按与 pack_calo_by_trigger 相同的规则，
// 从解包 caloTree 取出每个保留事例第一包的时间码
static std::vector<double> packedTimeCodesFromRaw(TTree* raw)
{
    std::vector<double> out;
    if (!raw || !raw->GetBranch("TimeCode") || !raw->GetBranch("TriggerID"))
        return out;

    int trigger_id = 0;
    int detector_id = 0;
    double time_code = 0.;
    std::vector<int>* cell_id = nullptr;
    raw->SetBranchAddress("TriggerID", &trigger_id);
    raw->SetBranchAddress("TimeCode", &time_code);
    raw->SetBranchAddress("CellID", &cell_id);
    if (raw->GetBranch("DetectorID"))
        raw->SetBranchAddress("DetectorID", &detector_id);

    int current_trigger = -1;
    int packet_count = 0;
    std::vector<int> fees;
    double group_tc = 0.;

    auto reset = [&]() {
        packet_count = 0;
        current_trigger = -1;
        fees.clear();
    };
    auto start = [&]() {
        current_trigger = trigger_id;
        packet_count = 1;
        fees = {feeKeyFromPacket(detector_id, cell_id)};
        group_tc = time_code;
    };
    auto keep_if_four_fees = [&]() {
        if (packet_count != 4)
            return;
        std::set<int> uniq;
        for (int f : fees) {
            if (f < 0) {
                reset();
                return;
            }
            uniq.insert(f);
        }
        if (uniq.size() == 4)
            out.push_back(group_tc);
        reset();
    };

    const Long64_t n = raw->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        raw->GetEntry(i);
        if (!cell_id)
            continue;
        if (packet_count == 0)
            start();
        else if (trigger_id == current_trigger) {
            ++packet_count;
            fees.push_back(feeKeyFromPacket(detector_id, cell_id));
        } else {
            reset();
            start();
        }
        if (packet_count == 4)
            keep_if_four_fees();
    }
    return out;
}

static void write_point(std::ofstream& ofs, const char* tag, const EvPt& p)
{
    ofs << "  " << tag
        << " 时间码=" << fmtTimeCode(p.has_tc, p.tc)
        << " 北京=" << fmtBeijing(p.has_tc, p.tc) << "\n";
}

static bool check_tree_gaps(TTree* tree, const char* tree_name, const char* rootfile,
                            const std::string& out_path, const std::vector<double>* ext_tc = nullptr)
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
    double time_code = 0.;
    const bool has_tc = tree->GetBranch("TimeCode") != nullptr;
    tree->SetBranchAddress("TriggerID", &trigger_id);
    if (has_tc)
        tree->SetBranchAddress("TimeCode", &time_code);

    const Long64_t nent = tree->GetEntries();
    ofs << "# entries=" << nent << "\n";
    if (nent < 1) {
        ofs << "# empty\n";
        ofs << "# segment_count=0\n";
        ofs << "# gap_count=0\n";
        ofs.close();
        std::cout << "Saved: " << out_path << std::endl;
        return true;
    }

    auto load_pt = [&](Long64_t i) {
        tree->GetEntry(i);
        EvPt p;
        p.entry = i;
        p.trigger = trigger_id;
        if (has_tc) {
            p.tc = time_code;
            p.has_tc = true;
        } else if (ext_tc && i >= 0 && static_cast<size_t>(i) < ext_tc->size()) {
            p.tc = (*ext_tc)[static_cast<size_t>(i)];
            p.has_tc = true;
        }
        return p;
    };

    EvPt seg_begin = load_pt(0);
    EvPt prev = seg_begin;
    int ngap = 0;
    int nseg = 0;

    ofs << "first TriggerID: " << seg_begin.trigger << "\n";
    if (seg_begin.trigger == 0)
        ofs << "OKAY: first TriggerID is 0\n";
    else
        ofs << "ERROR: first TriggerID is not 0\n";

    ofs << "# segments\n";

    auto close_segment = [&](const EvPt& seg_end) {
        ++nseg;
        ofs << "段" << nseg
            << " Entry=" << seg_begin.entry << ".." << seg_end.entry
            << " n=" << (seg_end.entry - seg_begin.entry + 1)
            << " TriggerID=" << seg_begin.trigger << ".." << seg_end.trigger << "\n";
        write_point(ofs, "起", seg_begin);
        write_point(ofs, "止", seg_end);
    };

    for (Long64_t i = 1; i < nent; ++i) {
        EvPt cur = load_pt(i);
        const bool wrap = (prev.trigger == 4095 && cur.trigger == 0);
        if (cur.trigger != prev.trigger + 1 && !wrap) {
            close_segment(prev);
            ++ngap;
            seg_begin = cur;
        }
        prev = cur;
    }
    close_segment(prev);

    ofs << "# segment_count=" << nseg << "\n";
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
    std::vector<double> packed_tc;
    const std::vector<double>* ext_tc = nullptr;
    TTree* calo_packed = dynamic_cast<TTree*>(fin_pack->Get("caloTree"));
    if (calo_packed && !calo_packed->GetBranch("TimeCode")) {
        TFile* fin_raw = TFile::Open(resultfile, "READ");
        if (fin_raw && !fin_raw->IsZombie()) {
            packed_tc = packedTimeCodesFromRaw(dynamic_cast<TTree*>(fin_raw->Get("caloTree")));
            fin_raw->Close();
            if (static_cast<Long64_t>(packed_tc.size()) == calo_packed->GetEntries())
                ext_tc = &packed_tc;
            else
                std::cerr << "find_lack_trigger: raw TimeCode count " << packed_tc.size()
                          << " != packed entries " << calo_packed->GetEntries()
                          << ", calo times omitted\n";
        }
    }
    check_tree_gaps(calo_packed, "caloTree", packfile, out_calo, ext_tc);
    fin_pack->Close();
}
