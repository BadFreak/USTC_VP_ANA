#ifndef COMREADER_HH
#define COMREADER_HH
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdint>
#include <set>
#include <bitset>
#include <unordered_map>
#include <algorithm>
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "yaml-cpp/yaml.h"

enum class WorkMode{
    WAVE,
    AMP,
    TEMP,
    MIX
};

enum class DetType{
    CALO,
    ITK
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        
        // 简单的哈希组合方法
        return h1 ^ h2;
        
        // 更好的哈希组合方法（来自Boost）
        // return h1 ^ (h2 << 1);
    }
};
struct WaveDrawEvent {
    unsigned int trigger_id = 0;
    long long time_code_key8 = 0;
};
struct CsIFirstCrossLate {
    unsigned int trigger_id = 0;
    long long time_code_key8 = 0;
};
class ComReader{
public:
static ComReader& getInstance() {
        // C++11 起，局部静态变量在首次调用时初始化，并自动保证线程安全
        static ComReader instance;
        return instance;
    }
    ComReader(const ComReader&) = delete;
    ComReader& operator=(const ComReader&) = delete;
    //Decode
    void decode(const std::string& filename, const std::string& yamlname="../config/triggerIDToBeDrawn.yaml");
    void decodeMix(const std::string& filename);
private:
    // 私有构造函数，防止外部实例化
    ComReader()=default;
    ~ComReader();
    // wave plot state
    TCanvas *caloWaveCanvas_[4] = {nullptr, nullptr, nullptr, nullptr};
    bool caloWaveCanvasInit_ = false;
    int caloWaveTriggerID_ = -1;
    long long caloWaveTimeCodeKey_ = -1;

    //Logger
    std::shared_ptr<spdlog::logger> logger = spdlog::stdout_color_mt("ComReader");
    //IO
    bool openFile(const std::string& filename);
    std::ifstream *file;
    TFile *fout;
    std::string oname; //output file name
    YAML::Node config;
    unsigned int yaml_package_mode_id;
    std::vector<WaveDrawEvent> yaml_draw_wave_events_;
    std::vector<unsigned int> yaml_calo_trigger_id; 
    std::vector<unsigned int> yaml_csi_trigger_id;
    std::vector<long long> yaml_calo_draw_time_code_;
    std::vector<long long> yaml_csi_draw_time_code_;
    bool yaml_time_filter_enable = false;
    bool yaml_csi_first_cross_hist_enable_ = false;
    double yaml_time_code_min = 0.;
    double yaml_time_code_max = 0.;
    TH1F* csi_first_cross_hist_[8] = {nullptr, nullptr, nullptr, nullptr,
                                       nullptr, nullptr, nullptr, nullptr};
    std::vector<CsIFirstCrossLate> csi_first_cross_late_;
    int nSkippedTime = 0;
    bool time_code_pass_valid_ = false;
    double time_code_pass_min_ = 0.;
    double time_code_pass_max_ = 0.;
    bool time_code_calo_pass_valid_ = false;
    double time_code_calo_pass_min_ = 0.;
    double time_code_calo_pass_max_ = 0.;
    bool time_code_csi_pass_valid_ = false;
    double time_code_csi_pass_min_ = 0.;
    double time_code_csi_pass_max_ = 0.;
    bool timeCodeInRange() const;
    void noteTimeCodePass(int feeid);
    void logTimeCodePassRange() const;
    void getOutputName(const std::string& prefix,const std::string& filename);

    //Initialize output tree for Calo, CsITK and HouseKeeping
    void initCaloTree();
    TTree *caloTree;
    void initCsITree();
    TTree *csiTree;
    void initHKTree(); //HouseKeeping Tree
    TTree *hkTree;

    //Find head
    bool findHead();
    bool readFEE();
    bool readCalo();
    bool readBuffer(char *b,const int& NCHN);
    bool readCaloWave(char *b);
    bool readCsI();
    bool readCsIWave(char *b);
    bool drawCsIWave(char *b, int triggerID);
    bool drawCaloWave(char *b, int triggerID);
    void flushCaloWavePlots();
    void fillCsIFirstCrossHist(const char* buf);
    void writeCsIFirstCrossHist();
    bool readHK();

    //Get set methods
    int getCellID(const int& chn_i); // Get Cell ID based on channel id and FEEID(private member)
    // 累加校验：buf 按 16bit 大端相加（Asum）
    void dosum(int &sum, char *buf, const int& size){
		for (int i = 0; i < size/2; i++) {
			auto c = static_cast<unsigned char>(buf[2*i]) << 8 | static_cast<unsigned char>(buf[2*i+1]);
			sum = (sum + static_cast<unsigned int>(c));
		}
    }
    // CRC-16/CCITT: poly x^16+x^12+x^5+1 (0x1021), init 0xFFFF
    // 校验范围：时间码之后(科学标志位起)到 CRC 字段之前
    static uint16_t crc16Ccitt(const char* data, int size, uint16_t crc = 0xFFFF);
    uint16_t calcPacketCrc(const char* body, int body_len) const;
    void checkPacketCrc(uint16_t recv_crc, const char* body, int body_len, const char* tag);
    char feeHeader_[12] = {0};

    //Package attributes
    unsigned int PackageLength; //Package length
    // int PackageID=-1;
    int TriggerID=-1;
    int TriggerID_csi=-1;
    int FEEID=-1;
    int nCaloHead=0;
    int nCsIHead=0;
    int nHKCaloHead=0;
    int nHKCsIHead=0;
    int nCaliHead=0;
    int nhead=0;
    int m_statusFeID=-1; // raw FE ID (high nibble of status byte): 5-8=Calo, 2=CsITK

    int Calo_EventID=0;
    int CsI_EventID=0;
    int Calo_EventCount=0;
    int asum=0; //Accumulation summation check
    
    int Package_Count=0, Package_Length=0;
    double Time_Code = 0.;
    int Science_Symbol=0;
    int DetectorID=0;
    int PackageID=-1;
    //Scientific Data
    //Calo
    // CryID-FEEID-MBID-GID-CHNID
    std::vector<int> Calo_CellID;
    std::vector<int> Calo_CellADC;
    std::vector<int> Calo_CellPLAT;
    int Calo_Trigger_Status=0;
    int Calo_TriggerID=0;
    int Calo_CRC=0;
    int Calo_Asum=0;

    //CsI
    // CryID-FEEID-MBID-GID-CHNID
    std::vector<int> CsI_CellID;
    std::vector<int> CsI_CellADC;
    std::vector<int> CsI_CellPLAT;
    int CsI_Trigger_Status=0;
    int CsI_TriggerID=0;
    int CsI_CRC=0;
    int CsI_Asum=0;

    //HouseKeeping Data（hkTree 每行 = 一个 HK 包：标量，不按 FEE 拼 vector）
    int csi_byte0=0, csi_byte1=0;
    //byte2~7 current, byte8~15 temperature
    int csi_correct_count=0, csi_incorrect_count=0;
    float hkC0=0.f, hkC1=0.f, hkC2=0.f;
    float hkT0=0.f, hkT1=0.f, hkT2=0.f, hkT3=0.f;

    //Calo data buffer
    std::vector<std::vector<uint8_t>> calo_buffer;

    //Clear function
    void clear();
    void clearHK();

    bool findTempHead();
    int findMixHead(); // return 0 for scientific, 1 for temperature
    bool readWave();
    bool readAmp();
    bool readTemp();
    bool readMixTemp();
    bool readAmpCsI();
    bool readWaveCsI();

    static constexpr unsigned char uc_0F=0x0F;
    static constexpr unsigned char uc_F0=0xF0;

    //Get output name
    
    //Input file and output TFile and TTree


    //Data
    // static constexpr int NChn=26;
    // float chn[NChn];
    // float chnplat[NChn];

    // float npackages=0.;
    // int nhead=0;
    // int nevents_skipped=0;

    // float MixTemperature; // Temperature for MixReader
    // float MixCurrent; // Current for MixReader

    // int MixTempEventCount=0;
    // int MixTime;
    // std::set<int> set_FEEID;

    // //Temperature Current Data
    
    // static constexpr int NFEE = 4;

    // int MixFEEID=-1;
    

    

    //Constants
    static constexpr double e=1.602176634e-19;
    static constexpr double lowgain=0.549; //ADC/fC
    static constexpr double highgain=5.614; //ADC/fC
    static constexpr double MIP_Edep=36.; //MeV
    static constexpr double gain=50.;

    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHash> channelMap_Calo = { //<FEEID,ChnID> <GID,CrystalID>
    {{3, 0}, {1, 1}}, {{3, 13}, {0, 1}}, {{4, 0}, {1, 1}}, {{4, 13}, {0, 1}},
    {{3, 1}, {1, 2}}, {{3, 14}, {0, 2}}, {{4, 1}, {1, 2}}, {{4, 14}, {0, 2}},
    {{1, 0}, {1, 3}},  {{1, 13}, {0, 3}}, {{2, 0}, {1, 3}},  {{2, 13}, {0, 3}},
    {{1, 1}, {1, 4}},  {{1, 14}, {0, 4}},  {{2, 1}, {1, 4}},  {{2, 14}, {0, 4}},
    {{1, 2}, {1, 5}},  {{1, 15}, {0, 5}},  {{2, 2}, {1, 5}},  {{2, 15}, {0, 5}},
    {{3, 2}, {1, 6}}, {{3, 15}, {0, 6}}, {{4, 2}, {1, 6}}, {{4, 15}, {0, 6}},
    {{3, 3}, {1, 7}}, {{3, 16}, {0, 7}}, {{4, 3}, {1, 7}}, {{4, 16}, {0, 7}},
    {{3, 4}, {1, 8}},  {{3, 17}, {0, 8}}, {{4, 4}, {1, 8}},  {{4, 17}, {0, 8}},
    {{1, 3}, {1, 9}},  {{1, 16}, {0, 9}}, {{2, 3}, {1, 9}},  {{2, 16}, {0, 9}},
    {{1, 4}, {1, 10}}, {{1, 17}, {0, 10}}, {{2, 4}, {1, 10}}, {{2, 17}, {0, 10}},
    {{3, 5}, {1, 11}},{{3, 18}, {0, 11}}, {{4, 5}, {1, 11}},{{4, 18}, {0, 11}},
    {{3, 6}, {1, 12}},{{3, 19}, {0, 12}},{{4, 6}, {1, 12}},{{4, 19}, {0, 12}},
    {{1, 5}, {1, 13}},{{1, 18}, {0, 13}},{{2, 5}, {1, 13}},{{2, 18}, {0, 13}},
    {{1, 6}, {1, 14}},{{1, 19}, {0, 14}},{{2, 6}, {1, 14}},{{2, 19}, {0, 14}},
    {{1, 7}, {1, 15}},{{1, 20}, {0, 15}}, {{2, 7}, {1, 15}},{{2, 20}, {0, 15}},
    {{3, 7}, {1, 16}},{{3, 20}, {0, 16}}, {{4, 7}, {1, 16}},{{4, 20}, {0, 16}},
    {{3, 8}, {1, 17}}, {{3, 21}, {0, 17}},{{4, 8}, {1, 17}}, {{4, 21}, {0, 17}},
    {{3, 9}, {1, 18}}, {{3, 22}, {0, 18}},{{4, 9}, {1, 18}}, {{4, 22}, {0, 18}},
    {{1, 8}, {1, 19}},{{1, 21}, {0, 19}},{{2, 8}, {1, 19}},{{2, 21}, {0, 19}},
    {{1, 9}, {1, 20}},{{1, 22}, {0, 20}}, {{2, 9}, {1, 20}},{{2, 22}, {0, 20}},
    {{3, 10}, {1, 21}}, {{3, 23}, {0, 21}}, {{4, 10}, {1, 21}}, {{4, 23}, {0, 21}},
    {{3, 11}, {1, 22}}, {{3, 24}, {0, 22}}, {{4, 11}, {1, 22}}, {{4, 24}, {0, 22}},
    {{1, 10}, {1, 23}},{{1, 23}, {0, 23}},{{2, 10}, {1, 23}},{{2, 23}, {0, 23}},
    {{1, 11}, {1, 24}},{{1, 24}, {0, 24}},{{2, 11}, {1, 24}},{{2, 24}, {0, 24}},
    {{1, 12}, {1, 25}},{{1, 25}, {0, 25}},{{2, 12}, {1, 25}},{{2, 25}, {0, 25}}
};

    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHash> channelMap_CsI = { //<FEEID,ChnID> <GID,CrystalID>
    {{3, 22}, {1, 1}}, {{3, 25}, {0, 1}}, {{4, 22}, {1, 1}}, {{4, 25}, {0, 1}},
    {{3, 14}, {1, 2}}, {{3, 21}, {0, 2}}, {{4, 14}, {1, 2}}, {{4, 21}, {0, 2}},
    {{1, 4}, {1, 3}},  {{1, 17}, {0, 3}}, {{2, 4}, {1, 3}},  {{2, 17}, {0, 3}},
    {{1, 2}, {1, 4}},  {{1, 9}, {0, 4}},  {{2, 2}, {1, 4}},  {{2, 9}, {0, 4}},
    {{1, 0}, {1, 5}},  {{1, 1}, {0, 5}},  {{2, 0}, {1, 5}},  {{2, 1}, {0, 5}},
    {{3, 20}, {1, 6}}, {{3, 23}, {0, 6}}, {{4, 20}, {1, 6}}, {{4, 23}, {0, 6}},
    {{3, 12}, {1, 7}}, {{3, 19}, {0, 7}}, {{4, 12}, {1, 7}}, {{4, 19}, {0, 7}},
    {{3, 6}, {1, 8}},  {{3, 17}, {0, 8}}, {{4, 6}, {1, 8}},  {{4, 17}, {0, 8}},
    {{1, 8}, {1, 9}},  {{1, 11}, {0, 9}}, {{2, 8}, {1, 9}},  {{2, 11}, {0, 9}},
    {{1, 6}, {1, 10}}, {{1, 3}, {0, 10}}, {{2, 6}, {1, 10}}, {{2, 3}, {0, 10}},
    {{3, 18}, {1, 11}},{{3, 7}, {0, 11}}, {{4, 18}, {1, 11}},{{4, 7}, {0, 11}},
    {{3, 10}, {1, 12}},{{3, 13}, {0, 12}},{{4, 10}, {1, 12}},{{4, 13}, {0, 12}},
    {{1, 10}, {1, 13}},{{1, 19}, {0, 13}},{{2, 10}, {1, 13}},{{2, 19}, {0, 13}},
    {{1, 14}, {1, 14}},{{1, 13}, {0, 14}},{{2, 14}, {1, 14}},{{2, 13}, {0, 14}},
    {{1, 20}, {1, 15}},{{1, 5}, {0, 15}}, {{2, 20}, {1, 15}},{{2, 5}, {0, 15}},
    {{3, 16}, {1, 16}},{{3, 5}, {0, 16}}, {{4, 16}, {1, 16}},{{4, 5}, {0, 16}},
    {{3, 8}, {1, 17}}, {{3, 11}, {0, 17}},{{4, 8}, {1, 17}}, {{4, 11}, {0, 17}},
    {{3, 4}, {1, 18}}, {{3, 15}, {0, 18}},{{4, 4}, {1, 18}}, {{4, 15}, {0, 18}},
    {{1, 16}, {1, 19}},{{1, 15}, {0, 19}},{{2, 16}, {1, 19}},{{2, 15}, {0, 19}},
    {{1, 22}, {1, 20}},{{1, 7}, {0, 20}}, {{2, 22}, {1, 20}},{{2, 7}, {0, 20}},
    {{3, 2}, {1, 21}}, {{3, 3}, {0, 21}}, {{4, 2}, {1, 21}}, {{4, 3}, {0, 21}},
    {{3, 0}, {1, 22}}, {{3, 9}, {0, 22}}, {{4, 0}, {1, 22}}, {{4, 9}, {0, 22}},
    {{1, 12}, {1, 23}},{{1, 25}, {0, 23}},{{2, 12}, {1, 23}},{{2, 25}, {0, 23}},
    {{1, 18}, {1, 24}},{{1, 23}, {0, 24}},{{2, 18}, {1, 24}},{{2, 23}, {0, 24}},
    {{1, 24}, {1, 25}},{{1, 21}, {0, 25}},{{2, 24}, {1, 25}},{{2, 21}, {0, 25}}
};

};

#endif
