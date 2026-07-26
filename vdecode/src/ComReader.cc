#include "ComReader.hh"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TApplication.h>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <TDirectory.h>
#include <cstdint>
#include <cstring>

static int daysInYear(int y) {
	return ((y % 4 == 0 && y % 100 != 0) || (y % 400 == 0)) ? 366 : 365;
}

// ITU-IEEE CRC-16/CCITT: poly 0x1021, init 0xFFFF
uint16_t ComReader::crc16Ccitt(const char* data, int size, uint16_t crc)
{
	static constexpr uint16_t kPoly = 0x1021;
	for (int i = 0; i < size; ++i) {
		crc ^= static_cast<uint16_t>(static_cast<uint8_t>(data[i])) << 8;
		for (int bit = 0; bit < 8; ++bit) {
			if (crc & 0x8000)
				crc = static_cast<uint16_t>((crc << 1) ^ kPoly);
			else
				crc = static_cast<uint16_t>(crc << 1);
		}
	}
	return crc;
}

uint16_t ComReader::calcPacketCrc(const char* body, int body_len) const
{
	// 协议：CRC 范围 = 时间码之后 ~ CRC 之前
	// feeHeader_[10]=科学数据包标志位, [11]=从机|包类型；body 含科学数据+触发号，不含 CRC/Asum
	uint16_t crc = crc16Ccitt(feeHeader_ + 10, 2);
	return crc16Ccitt(body, body_len, crc);
}

void ComReader::checkPacketCrc(uint16_t recv_crc, const char* body, int body_len, const char* tag)
{
	const uint16_t calc_crc = calcPacketCrc(body, body_len);
	if (recv_crc != calc_crc)
		logger->error("CRC check failed [{}]: recv=0x{:04X} calc=0x{:04X}", tag, recv_crc, calc_crc);
}

// 基准时刻（北京时间）：2025-01-01 00:00:00；TimeCode 为相对该时刻的秒数
static const int kTimeBaseY = 2025, kTimeBaseMo = 1, kTimeBaseD = 1;

static long long timeCodeKey8(double timeCode)
{
	// TimeCode 毫秒整数的前八位 = 秒（与 caloTimeCode.png 横轴整数部分一致）
	const auto tc_ms = static_cast<long long>(std::llround(timeCode * 1000.0));
	return tc_ms / 1000;
}

static TString waveTag(int triggerID, double timeCode)
{
	return Form("TriggerID_%d_TimeCode_%lld", triggerID, timeCodeKey8(timeCode));
}

static bool waveDrawMatch(const std::vector<WaveDrawEvent>& events,
                          const std::vector<unsigned int>& trigger_ids,
                          const std::vector<long long>& time_code_keys, unsigned int trigger_id,
                          double time_code)
{
	const long long key = timeCodeKey8(time_code);
	if (!events.empty()) {
		for (const auto& e : events) {
			if (e.trigger_id == trigger_id && e.time_code_key8 == key)
				return true;
		}
		return false;
	}
	if (!time_code_keys.empty()) {
		return std::find(time_code_keys.begin(), time_code_keys.end(), key) !=
		       time_code_keys.end();
	}
	if (!trigger_ids.empty()) {
		return std::find(trigger_ids.begin(), trigger_ids.end(), trigger_id) !=
		       trigger_ids.end();
	}
	return false;
}

static constexpr double kWaveYHalfRange = 200.;

static void styleWaveAxes(TGraph* gr)
{
	if (!gr)
		return;
	gr->GetXaxis()->SetLabelSize(0.05);
	gr->GetXaxis()->SetNdivisions(520);
	gr->GetYaxis()->SetLabelSize(0.055);
}

// static const int kMhWaveThreshold[25] = {
//     919, 1187, 1191, 1321, 1062, 1088, 1145, 1070, 1270, 1335, 1003, 1208, 1296, 1286, 1269,
//     919, 1164, 985, 1311, 1103, 1217, 1106, 1209, 1128, 1138};
static const int kMhWaveThreshold[25] = {
    941, 1209, 1198, 1340, 1094, 1107, 1156, 1091, 1287, 1350, 1031, 1222, 1297, 1297, 1287,
    939, 1186, 1012, 1329, 1118, 1247, 1128, 1221, 1143, 1154};
// static const int kBhWaveThreshold[25] = {
//     1094, 1156, 1462, 1476, 1595, 1139, 1078, 1151, 1450, 1317, 1301, 926, 1501, 1317, 1449,
//     1145, 1165, 1369, 1369, 1470, 951, 1040, 1310, 1177, 1698};
static const int kBhWaveThreshold[25] = {
    1111, 1170, 1447, 1453, 1582, 1144, 1080, 1159, 1426, 1297, 1320, 938, 1478, 1297, 1434,
    1158, 1180, 1386, 1346, 1451, 957, 1064, 1305, 1158, 1687};
// static const int kCsIWaveThreshold[8] = {1205, 1264, 1077, 1244, 1098, 1254, 1157, 1151};
static const int kCsIWaveThreshold[8] = {1168, 1207, 990, 1215, 1034, 1209, 1120, 1121};
static constexpr int kCsIFirstCrossRunLength = 10; // 首过阈点 + 其后 9 点共 10 点须全部过阈
static constexpr int kCsIFirstCrossLateIndex = 65; // 8 通道首点位置取 max，> 65 则记录该事例

static unsigned int csIWaveSample(const char* b, size_t chn, size_t sp)
{
	const auto d1 = static_cast<unsigned char>(b[chn * 128 * 2 + sp * 2]);
	const auto d2 = static_cast<unsigned char>(b[chn * 128 * 2 + sp * 2 + 1]);
	return static_cast<unsigned int>(d1 << 8 | d2);
}

// 每通道独立：找第一个 ADC > threshold 的点；该点及其后 9 点须全部过阈，否则不记录（-1）
static int csIFirstConsecutiveCrossIndex(const char* b, size_t chn, int threshold)
{
	const auto thr = static_cast<unsigned int>(threshold);
	for (size_t sp = 0; sp < 128; ++sp) {
		if (csIWaveSample(b, chn, sp) <= thr)
			continue;
		if (sp + kCsIFirstCrossRunLength > 128)
			return -1;
		for (int k = 1; k < kCsIFirstCrossRunLength; ++k) {
			if (csIWaveSample(b, chn, sp + k) <= thr)
				return -1;
		}
		return static_cast<int>(sp);
	}
	return -1;
}

static void drawThresholdLineAndLegend(int threshold)
{
	auto* line = new TLine(0., threshold, 127., threshold);
	line->SetLineColor(kRed);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->Draw("same");

	auto* leg = new TLegend(0.55, 0.78, 0.95, 0.94);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.06);
	leg->AddEntry(static_cast<TObject*>(nullptr), Form("threshold: %d", threshold), "");
	leg->Draw();
}

static void drawWaveThresholdLegend(int canvasIdx, int cry)
{
	if (canvasIdx != 0 && canvasIdx != 2)
		return;
	const int threshold =
	    (canvasIdx == 0) ? kMhWaveThreshold[cry] : kBhWaveThreshold[cry];
	drawThresholdLineAndLegend(threshold);
}

static void drawCsIWaveThresholdLegend(size_t chn)
{
	if (chn >= 8)
		return;
	drawThresholdLineAndLegend(kCsIWaveThreshold[chn]);
}

static bool cdWavePlotDir(TFile* fout, const char* parent)
{
	if (!fout->GetDirectory(parent) && !fout->mkdir(parent))
		return false;
	fout->cd(parent);
	return true;
}
static const int kTimeBaseH = 0, kTimeBaseMi = 0, kTimeBaseS = 0;
static double secFromYearStart(int y, int mo, int d, int H, int M, int S) {
	static const int cum[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
	int doy = cum[mo - 1] + (d - 1);
	if (mo > 2 && daysInYear(y) == 366)
		doy++;
	return doy * 86400. + H * 3600. + M * 60. + S;
}

static double secSinceTimeBase(int y, int mo, int d, int H, int M, int S) {
	double t = secFromYearStart(y, mo, d, H, M, S);
	if (y > kTimeBaseY) {
		for (int yr = kTimeBaseY; yr < y; ++yr)
			t += daysInYear(yr) * 86400.;
	} else if (y < kTimeBaseY) {
		for (int yr = y; yr < kTimeBaseY; ++yr)
			t -= daysInYear(yr) * 86400.;
	}
	static const double kBase =
	    secFromYearStart(kTimeBaseY, kTimeBaseMo, kTimeBaseD, kTimeBaseH, kTimeBaseMi, kTimeBaseS);
	return t - kBase;
}

static std::string secSinceBaseToBeijingString(double sec) {
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
	int S = static_cast<int>(sod - M * 60.);
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
	char buf[32];
	std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d", y, mo, d, H, M, S);
	return std::string(buf);
}

bool ComReader::timeCodeInRange() const {
	if (!yaml_time_filter_enable)
		return true;
	return Time_Code >= yaml_time_code_min && Time_Code <= yaml_time_code_max;
}

void ComReader::noteTimeCodePass(int feeid)
{
	auto update = [](bool& valid, double& minv, double& maxv, double tc) {
		if (!valid) {
			minv = maxv = tc;
			valid = true;
		} else {
			if (tc < minv)
				minv = tc;
			if (tc > maxv)
				maxv = tc;
		}
	};
	update(time_code_pass_valid_, time_code_pass_min_, time_code_pass_max_, Time_Code);
	if (feeid >= 5 && feeid <= 8)
		update(time_code_calo_pass_valid_, time_code_calo_pass_min_, time_code_calo_pass_max_,
		       Time_Code);
	else if (feeid == 2)
		update(time_code_csi_pass_valid_, time_code_csi_pass_min_, time_code_csi_pass_max_, Time_Code);
}

static void writeTimeRangeTxt(const std::string& path, bool valid, double tmin, double tmax)
{
	if (!valid)
		return;
	std::ofstream ofs(path);
	if (!ofs)
		return;
	const std::string bmin = secSinceBaseToBeijingString(tmin);
	const std::string bmax = secSinceBaseToBeijingString(tmax);
	ofs << "time_code_min_pass " << tmin << "\n";
	ofs << "time_code_max_pass " << tmax << "\n";
	ofs << "beijing_min " << bmin << "\n";
	ofs << "beijing_max " << bmax << "\n";
}

void ComReader::logTimeCodePassRange() const
{
	if (!time_code_pass_valid_) {
		logger->info("TimeCode pass range: (no packet passed filter)");
	} else {
		const std::string tmin = secSinceBaseToBeijingString(time_code_pass_min_);
		const std::string tmax = secSinceBaseToBeijingString(time_code_pass_max_);
		logger->info("TimeCode pass range: [{}, {}] -> Beijing {} .. {}", time_code_pass_min_,
		             time_code_pass_max_, tmin, tmax);
	}
	if (!time_code_calo_pass_valid_) {
		logger->info("Calo TimeCode pass range: (no packet passed filter)");
	} else {
		logger->info("Calo TimeCode pass range: [{}, {}] -> Beijing {} .. {}",
		             time_code_calo_pass_min_, time_code_calo_pass_max_,
		             secSinceBaseToBeijingString(time_code_calo_pass_min_),
		             secSinceBaseToBeijingString(time_code_calo_pass_max_));
	}
	if (!time_code_csi_pass_valid_) {
		logger->info("CsI TimeCode pass range: (no packet passed filter)");
	} else {
		logger->info("CsI TimeCode pass range: [{}, {}] -> Beijing {} .. {}",
		             time_code_csi_pass_min_, time_code_csi_pass_max_,
		             secSinceBaseToBeijingString(time_code_csi_pass_min_),
		             secSinceBaseToBeijingString(time_code_csi_pass_max_));
	}

	std::string base = oname;
	const auto dot = base.rfind('.');
	if (dot != std::string::npos)
		base = base.substr(0, dot);

	writeTimeRangeTxt(base + "_time_range.txt", time_code_pass_valid_, time_code_pass_min_,
	                  time_code_pass_max_);
	writeTimeRangeTxt(base + "_calo_time_range.txt", time_code_calo_pass_valid_,
	                  time_code_calo_pass_min_, time_code_calo_pass_max_);
	writeTimeRangeTxt(base + "_csi_time_range.txt", time_code_csi_pass_valid_,
	                  time_code_csi_pass_min_, time_code_csi_pass_max_);
}

bool ComReader::openFile(const std::string& filename){
    file = new std::ifstream(filename,std::ios::binary);
    if (!file->is_open()) {
        std::cerr << "Can not open file: "<<filename << std::endl;
        return false;
    }
    return true;
}

bool ComReader::findHead()
{
	char buffer[1];
	while (file->rdbuf()->sgetn(buffer, sizeof(buffer)) > 0)
	{
		unsigned char byte = static_cast<unsigned char>(buffer[0]);
		if (byte == 0xeb)
		{
			char buffer2[1];
			if (file->rdbuf()->sgetn(buffer2, sizeof(buffer2)) > 0)
			{
				unsigned char b1 = static_cast<unsigned char>(buffer2[0]);
				if (b1 == 0x90)
				{
					char buffer_length[12];// bag count high byte and low byte, bag length high byte and low byte
					file->rdbuf()->sgetn(buffer_length, sizeof(buffer_length));
					// for (int i = 0; i < 12; i++) {
					// 	auto b = static_cast<unsigned char>(buffer_length[i]);
					// 	std::cout << std::hex << static_cast<unsigned int>(b) << " ";
					// }
					// std::cout << std::endl;
					unsigned char length_0 = static_cast<unsigned char>(buffer_length[2]);
					unsigned char length_1 = static_cast<unsigned char>(buffer_length[3]);
					PackageLength = static_cast<unsigned int>(length_0 << 8 | length_1);
					auto char_0 = static_cast<unsigned char>(buffer_length[11]);
					FEEID = static_cast<unsigned int>((char_0 & uc_F0) >> 4); //FEEID 前四位是FEEID 后四位是包类型标识
					PackageID = static_cast<unsigned int>((char_0 & uc_0F) >> 0); 
					// std::cout << "FEEID : " << FEEID << std::endl;
					//std::cout<<PackageLength<<std::endl;
					// 科学数据包长度：44/116=幅度，2060/6668=波形；与 PackageID 无关
					if(FEEID == 1 || FEEID==3 || FEEID==4 || FEEID==9){ // need skip
						file->seekg(PackageLength-6, std::ios::cur); // skip the rest -8+2
					} else if (PackageLength == 116 || PackageLength == 44 || PackageLength==6668 || PackageLength==2060)
					{ 
						std::cout << "PackageID:" << PackageID << std::endl;
						if (PackageID == yaml_package_mode_id) { // Scientific data package
							file->seekg(-12, std::ios::cur);// back to 0x90 next
							nhead++;
							asum = 0xeb90;
							if (!readFEE()) {logger->error("Failed to read FEE data");}
						} else {
							file->seekg(PackageLength-6	, std::ios::cur); // back to 0x90 next
						}
					} else {
						if (0) {
							unsigned char b2 = static_cast<unsigned char>(buffer_length[0]);
							if ((b2 & 0x0F) == 0x0B)
							{ // HK (temperature/status) package: feeFlag=0x0B
								file->seekg(-3, std::ios::cur);
								unsigned int x = (b2 >> 4) & 0x0F;
								m_statusFeID = static_cast<int>(x);
								FEEID = m_statusFeID;
								if (x >= 5 && x <= 8)
								{ // Calo: feID 5,6,7,8 -> FEEID 1,2,3,4
									// std::cout << "Read the temperature" << std::endl;
									nHKCaloHead++;
									if (!readHK()) {logger->error("Failed to read Calo HK data");}
								} else if (x == 2) { // CsITK: single FEE, use slot 1
									nHKCsIHead++;
									if (!readHK()) {logger->error("Failed to read CsITK HK data");}
								}
								// else: HK flag but feID not Calo/CsITK — keep scanning, do not return false
							}
						}
						// else: length not science and not HK 0x0B — keep scanning (return false would stop decode())
					}
				} else {
					file->seekg(-1, std::ios::cur); // Not a scientific data head, seek back
				}
			} else {
				// EOF: no second byte after 0xeb — leave loop via while condition, then return false below
				break;
			}
		}
	}
	return false;
}

bool ComReader::readFEE(){ // Read length, FEEID and other info
	char buffer[12];
	file->rdbuf()->sgetn(buffer,sizeof(buffer));
	std::memcpy(feeHeader_, buffer, sizeof(feeHeader_));
	//Accumulation summation
	dosum(asum,buffer,12);
	auto char_0 = static_cast<unsigned char>(buffer[0]);
	auto char_1 = static_cast<unsigned char>(buffer[1]);
	Package_Count = static_cast<unsigned int>(char_0<<8 | char_1);
	char_0 = static_cast<unsigned char>(buffer[2]);
	char_1 = static_cast<unsigned char>(buffer[3]);
	Package_Length = static_cast<unsigned int>(char_0<<8 | char_1);
	unsigned int sec = 0;
	for (int i = 4; i < 8; i++) {
		char_0 = static_cast<unsigned char>(buffer[i]);
		sec = sec << 8 | char_0;
	}
	char_0 = static_cast<unsigned char>(buffer[8]);
	char_1 = static_cast<unsigned char>(buffer[9]);
	unsigned int ms = char_0 << 8 | char_1;
	Time_Code = sec + ms / 1000.;
	if (!timeCodeInRange()) {
		file->seekg(PackageLength - 6, std::ios::cur);
		nSkippedTime++;
		return true;
	}
	char_0 = static_cast<unsigned char>(buffer[10]);
	Science_Symbol = static_cast<unsigned int>(char_0);
	char_0 = static_cast<unsigned char>(buffer[11]);
	FEEID = static_cast<unsigned int>((char_0 & uc_F0) >> 4); //FEEID 前四位是FEEID 后四位是包类型标识
	noteTimeCodePass(FEEID);
	DetectorID = FEEID;
	PackageID = static_cast<unsigned int>(char_0 & uc_0F); //Package class
	if(FEEID >=5 && FEEID <=8 ){ // Calo
		nCaloHead++;
		if(!readCalo()){
			logger->error("Failed to read Calo data");
		}
	} else if(FEEID == 2){ // CsI
		nCsIHead++;
		if(!readCsI()){
			logger->error("Failed to read CsI data");
		}
	} else {
		// logger->info("Unknown FEEID: {}", FEEID);
	}
	return true;
}

void ComReader::clear(){
	if (Calo_CellID.size() > 0)Calo_CellID.clear();
	if (Calo_CellADC.size() > 0)Calo_CellADC.clear();
	if (Calo_CellPLAT.size() > 0)Calo_CellPLAT.clear();
	if (CsI_CellID.size() > 0)CsI_CellID.clear();
	if (CsI_CellADC.size() > 0)CsI_CellADC.clear();
	if (CsI_CellPLAT.size() > 0)CsI_CellPLAT.clear();
}

void ComReader::clearHK(){
	// hkTree 已改为每包标量行，无 vector 状态可清
}

bool ComReader::readCsI(){
	if(PackageLength==44){
		// std::cout << "============== Amp ==================" << std::endl;
		const int readsize = PackageLength - 6;//-8+2
		char buf[readsize];
		file->rdbuf()->sgetn(buf,sizeof(buf));
		dosum(asum,buf,readsize-2);
		//Trigger ID
		auto char_0 = static_cast<unsigned char>(buf[readsize-6]);
		auto char_1 = static_cast<unsigned char>(buf[readsize-5]);
		CsI_Trigger_Status = static_cast<unsigned int>((char_0 & uc_F0) >> 4);
		// if(CsI_Trigger_Status)std::cout<<"CsI Trigger Status: "<<CsI_Trigger_Status<<std::endl;
		CsI_TriggerID = static_cast<unsigned int>((char_0 & uc_0F)<<8 | char_1);
		//CRC check
		char_0 = static_cast<unsigned char>(buf[readsize-4]);
		char_1 = static_cast<unsigned char>(buf[readsize-3]);
		CsI_CRC = static_cast<unsigned int>(char_0<<8 | char_1);
		//Accumulation summation check
		char_0 = static_cast<unsigned char>(buf[readsize-2]);
		char_1 = static_cast<unsigned char>(buf[readsize-1]);
		CsI_Asum = static_cast<unsigned int>(char_0<<8 | char_1);

		checkPacketCrc(static_cast<uint16_t>(CsI_CRC), buf, readsize - 4, "CsI-amp");
		uint16_t recv_asum = static_cast<uint16_t>(CsI_Asum);
		uint16_t calc_asum = static_cast<uint16_t>(asum & 0x0000FFFF);
		if (recv_asum != calc_asum) {
			logger->error("Accumulation summation check failed: {} != {}", recv_asum, calc_asum);
		}
		readBuffer(buf, 8);
	} else if(PackageLength==2060) {
		// std::cout << "============== Wave ==================" << std::endl;
		const int readsize = PackageLength - 6;//-8+2
		char buf[readsize];
		file->rdbuf()->sgetn(buf,sizeof(buf));
		dosum(asum,buf,readsize-2);
		//Trigger ID
		auto char_0 = static_cast<unsigned char>(buf[readsize-6]);
		auto char_1 = static_cast<unsigned char>(buf[readsize-5]);
		CsI_Trigger_Status = static_cast<unsigned int>((char_0 & uc_F0) >> 4);
		CsI_TriggerID = static_cast<unsigned int>((char_0 & uc_0F)<<8 | char_1);
		//CRC check
		char_0 = static_cast<unsigned char>(buf[readsize-4]);
		char_1 = static_cast<unsigned char>(buf[readsize-3]);
		CsI_CRC = static_cast<unsigned int>(char_0<<8 | char_1);
		//Accumulation summation check
		char_0 = static_cast<unsigned char>(buf[readsize-2]);
		char_1 = static_cast<unsigned char>(buf[readsize-1]);
		CsI_Asum = static_cast<unsigned int>(char_0<<8 | char_1);
		checkPacketCrc(static_cast<uint16_t>(CsI_CRC), buf, readsize - 4, "CsI-wave");
		uint16_t recv_asum = static_cast<uint16_t>(CsI_Asum);
		uint16_t calc_asum = static_cast<uint16_t>(asum & 0x0000FFFF);
		if (recv_asum != calc_asum) {
			logger->error("Accumulation summation check failed: {} != {}", recv_asum, calc_asum);
		}
		if (yaml_csi_first_cross_hist_enable_)
			fillCsIFirstCrossHist(buf);
		if (waveDrawMatch(yaml_draw_wave_events_, yaml_csi_trigger_id, yaml_csi_draw_time_code_,
		                  CsI_TriggerID, Time_Code)) {
			drawCsIWave(buf, CsI_TriggerID);
		}
		readCsIWave(buf);
	}
	csiTree->Fill();
	clear();
	CsI_EventID++;
	// std::cout << "CsI_EventID: TriggerID " << std::dec << CsI_EventID<<": " <<std::hex << CsI_TriggerID << std::endl;
	// std::cout << "================== CURRENT POSITION: " << file->tellg() << std::endl;
	return true;
}

bool ComReader::readCalo(){
	if(PackageLength==116){
		const int readsize = PackageLength - 6;//-8+2e
		char buf[readsize];
		file->rdbuf()->sgetn(buf,sizeof(buf));
		dosum(asum,buf,readsize-2);
		//Trigger ID
		auto char_0 = static_cast<unsigned char>(buf[104]);
		Calo_Trigger_Status = static_cast<unsigned int>((char_0 & uc_F0) >> 4);
		auto char_1 = static_cast<unsigned char>(buf[105]);
		Calo_TriggerID = static_cast<unsigned int>((char_0 & uc_0F)<<8 | char_1);
		//CRC check
		char_0 = static_cast<unsigned char>(buf[106]);
		char_1 = static_cast<unsigned char>(buf[107]);
		Calo_CRC = static_cast<unsigned int>(char_0<<8 | char_1);
		//Accumulation summation check
		char_0 = static_cast<unsigned char>(buf[108]);
		char_1 = static_cast<unsigned char>(buf[109]);
		Calo_Asum = static_cast<unsigned int>(char_0<<8 | char_1);
		checkPacketCrc(static_cast<uint16_t>(Calo_CRC), buf, readsize - 4, "Calo-amp");
		uint16_t recv_asum = static_cast<uint16_t>(Calo_Asum);
		uint16_t calc_asum = static_cast<uint16_t>(asum & 0x0000FFFF);
		if (recv_asum != calc_asum) {
			logger->error("Accumulation summation check failed: {} != {}", recv_asum, calc_asum);
		}
		readBuffer(buf, 26);
		Calo_EventCount++;
	}
	else if(PackageLength==6668){
		const int readsize = PackageLength - 6;//-8+2
		char buf[readsize];
		file->rdbuf()->sgetn(buf,sizeof(buf));
		dosum(asum,buf,readsize-2);
		//Trigger ID
		auto char_0 = static_cast<unsigned char>(buf[6656]);
		auto char_1 = static_cast<unsigned char>(buf[6657]);
		Calo_Trigger_Status = static_cast<unsigned int>((char_0 & uc_F0) >> 4);
		Calo_TriggerID = static_cast<unsigned int>((char_0 & uc_0F)<<8 | char_1);
		//CRC check
		char_0 = static_cast<unsigned char>(buf[6658]);
		char_1 = static_cast<unsigned char>(buf[6659]);
		Calo_CRC = static_cast<unsigned int>(char_0<<8 | char_1);
		//Accumulation summation check
		char_0 = static_cast<unsigned char>(buf[6660]);
		char_1 = static_cast<unsigned char>(buf[6661]);
		Calo_Asum = static_cast<unsigned int>(char_0<<8 | char_1);
		checkPacketCrc(static_cast<uint16_t>(Calo_CRC), buf, readsize - 4, "Calo-wave");
		uint16_t recv_asum = static_cast<uint16_t>(Calo_Asum);
		uint16_t calc_asum = static_cast<uint16_t>(asum & 0x0000FFFF);
		if (recv_asum != calc_asum) {
			logger->error("Accumulation summation check failed: {} != {}", recv_asum, calc_asum);
		}
		if (waveDrawMatch(yaml_draw_wave_events_, yaml_calo_trigger_id, yaml_calo_draw_time_code_,
		                  static_cast<unsigned int>(Calo_TriggerID), Time_Code)) {
			drawCaloWave(buf, Calo_TriggerID);
		}
		readCaloWave(buf); 
		Calo_EventCount++;
	}
	else{
		logger->error("Unexpected Calo package length: {}", PackageLength);
		return false;
	}
	caloTree->Fill();
	clear();
	// std::cout << "Calo_EventID: TriggerID " << std::dec << Calo_EventID;
	// std::cout << ": " <<std::hex << Calo_TriggerID << std::endl;
	if (Calo_EventCount == 4) {
		Calo_EventID++;
		Calo_EventCount = 0;
	}
	asum=0;
	return true;
}

int ComReader::getCellID(const int& chn_i){ //Get cell id based on channel id
	std::pair<int, int> gid_cryid = std::pair<int, int>(0, 0);
	int feeid=FEEID-4;
	if (channelMap_Calo.count(std::pair<int, int>(feeid, chn_i)) == 0)
	{
		if (feeid == 3 || feeid == 4)
		{
			if (chn_i == 12 || chn_i == 25)
			{
				return -1;
			}
		}
	}
	gid_cryid = channelMap_Calo[std::pair<int, int>(feeid, chn_i)];
	if (gid_cryid.second == 0)
	{
		std::cerr << "No channel map for FEEID: " << feeid << " chn: " << chn_i << std::endl;
	}
	int tmp_cellid = gid_cryid.second * 100000 + 10000 * feeid + 1000 * (feeid % 2) + 100 * (gid_cryid.first) + chn_i;
	return tmp_cellid;
}

void ComReader::fillCsIFirstCrossHist(const char* buf)
{
	if (!yaml_csi_first_cross_hist_enable_)
		return;
	int max_first_sp = -1;
	for (size_t chn = 0; chn < 8; ++chn) {
		const int first_sp =
		    csIFirstConsecutiveCrossIndex(buf, chn, kCsIWaveThreshold[chn]);
		if (first_sp >= 0 && csi_first_cross_hist_[chn])
			csi_first_cross_hist_[chn]->Fill(static_cast<double>(first_sp));
		if (first_sp > max_first_sp)
			max_first_sp = first_sp;
	}
	if (max_first_sp > kCsIFirstCrossLateIndex) {
		const long long tc_key8 = timeCodeKey8(Time_Code);
		csi_first_cross_late_.push_back({static_cast<unsigned int>(CsI_TriggerID), tc_key8});
		logger->info("CsI first-cross late: TriggerID={} TimeCode={} max_first_sp={}",
		             CsI_TriggerID, tc_key8, max_first_sp);
	}
}

void ComReader::writeCsIFirstCrossHist()
{
	if (!yaml_csi_first_cross_hist_enable_)
		return;

	std::string base = oname;
	const auto dot = base.rfind('.');
	if (dot != std::string::npos)
		base = base.substr(0, dot);
	if (!csi_first_cross_late_.empty()) {
		const std::string late_path = base + "_csi_first_cross_late.txt";
		std::ofstream ofs(late_path);
		if (ofs) {
			ofs << "# TriggerID TimeCode_key8 (max over 8 chn: first over-thr point + next 9 all over thr, start > "
			    << kCsIFirstCrossLateIndex << ")\n";
			for (const auto& e : csi_first_cross_late_)
				ofs << e.trigger_id << ' ' << e.time_code_key8 << '\n';
			logger->info("CsI first-cross late: {} entries -> {}", csi_first_cross_late_.size(),
			             late_path);
		} else {
			logger->warn("CsI first-cross late: failed to write {}", late_path);
		}
	}
	csi_first_cross_late_.clear();

	long long total_entries = 0;
	for (int i = 0; i < 8; ++i) {
		if (csi_first_cross_hist_[i])
			total_entries += static_cast<long long>(csi_first_cross_hist_[i]->GetEntries());
	}
	if (total_entries == 0) {
		logger->info("CsI first-cross hist: no crossing entries, skip canvas");
		for (int i = 0; i < 8; ++i) {
			delete csi_first_cross_hist_[i];
			csi_first_cross_hist_[i] = nullptr;
		}
		return;
	}

	if (!fout->GetDirectory("CsIFirstCrossHist"))
		fout->mkdir("CsIFirstCrossHist");
	fout->cd("CsIFirstCrossHist");

	auto* canvas = new TCanvas("cCsI_first_cross",
	                           Form("CsI first over-thr + %d following all over thr",
	                                kCsIFirstCrossRunLength - 1),
	                           1600, 800);
	canvas->Divide(4, 2);
	for (int i = 0; i < 8; ++i) {
		canvas->cd(i + 1);
		if (!csi_first_cross_hist_[i])
			continue;
		csi_first_cross_hist_[i]->SetLineColor(kBlue + 1);
		csi_first_cross_hist_[i]->SetLineWidth(2);
		csi_first_cross_hist_[i]->Draw("HIST");
	}
	canvas->Write();
	delete canvas;

	for (int i = 0; i < 8; ++i) {
		if (csi_first_cross_hist_[i]) {
			csi_first_cross_hist_[i]->Write();
			delete csi_first_cross_hist_[i];
			csi_first_cross_hist_[i] = nullptr;
		}
	}
	logger->info("CsI first-cross hist written to CsIFirstCrossHist/ ({} fills)", total_entries);
}

bool ComReader::drawCsIWave(char *b, int triggerID) {
    const TString tag = waveTag(triggerID, Time_Code);
    TCanvas *cWave = new TCanvas(tag, tag, 1600, 400);
    cWave->Divide(4, 2);
    TGraph *grWave[8] = {nullptr}; // 显式初始化为 nullptr

    for (size_t chn_i = 0; chn_i < 8; chn_i++) {
        cWave->cd(chn_i + 1);
        grWave[chn_i] = new TGraph();
        grWave[chn_i]->SetName(Form("%s_Chn%zu", tag.Data(), chn_i));
        grWave[chn_i]->SetTitle(Form("%s Chn%zu;Sampling Point;ADC", tag.Data(), chn_i));

        double y0 = 0.;
        for (size_t sp = 0; sp < 128; sp++) {
            auto d1 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2]);
            auto d2 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2 + 1]);
            unsigned int sp_data = static_cast<unsigned int>(d1 << 8 | d2);
            if (sp == 0)
                y0 = sp_data;
            grWave[chn_i]->SetPoint(sp, sp, sp_data);
        }
        grWave[chn_i]->SetMarkerStyle(20);
        grWave[chn_i]->SetMarkerSize(0.5);
        grWave[chn_i]->GetYaxis()->SetRangeUser(y0 - kWaveYHalfRange, y0 + kWaveYHalfRange);
        grWave[chn_i]->Draw("APL");
        styleWaveAxes(grWave[chn_i]);
        drawCsIWaveThresholdLegend(chn_i);
    }

	if (!cdWavePlotDir(fout, "CsIWavePlot")) {
		std::cerr << "Error: Failed to cd CsIWavePlot" << std::endl;
		cWave->Delete();
		return false;
	}

    cWave->Write();
    delete cWave; // 删除画布对象

    for (auto &gr : grWave) {
        if (gr) {
            delete gr; // 删除 TGraph 对象
        }
    }
    return true;
}

void ComReader::flushCaloWavePlots() {
	if (!caloWaveCanvasInit_)
		return;

	if (!cdWavePlotDir(fout, "CaloWavePlot")) {
		std::cerr << "Error: Failed to cd CaloWavePlot" << std::endl;
		return;
	}
	for (int k = 0; k < 4; ++k) {
		if (caloWaveCanvas_[k]) {
			caloWaveCanvas_[k]->Write();
			delete caloWaveCanvas_[k];
			caloWaveCanvas_[k] = nullptr;
		}
	}
	caloWaveCanvasInit_ = false;
	caloWaveTriggerID_ = -1;
	caloWaveTimeCodeKey_ = -1;
}

bool ComReader::drawCaloWave(char *b, int triggerID) {
	const long long tc_key8 = timeCodeKey8(Time_Code);
	const unsigned int trig = static_cast<unsigned int>(triggerID);
	if (caloWaveCanvasInit_ &&
	    (tc_key8 != caloWaveTimeCodeKey_ || trig != static_cast<unsigned int>(caloWaveTriggerID_)))
		flushCaloWavePlots();

	const int plotTriggerID = caloWaveCanvasInit_ ? caloWaveTriggerID_ : triggerID;
	const double plotTimeCode =
	    caloWaveCanvasInit_ ? static_cast<double>(caloWaveTimeCodeKey_) : Time_Code;
	const TString tag = waveTag(plotTriggerID, plotTimeCode);
	if (!caloWaveCanvasInit_) {
		static const char *tags[] = {"mh", "ml", "bh", "bl"};
		for (int k = 0; k < 4; ++k) {
			caloWaveCanvas_[k] = new TCanvas(Form("%s_%s", tag.Data(), tags[k]),
			                                 Form("%s_%s", tag.Data(), tags[k]), 1600, 800);
			caloWaveCanvas_[k]->Divide(5, 5);
		}
		caloWaveCanvasInit_ = true;
		caloWaveTriggerID_ = triggerID;
		caloWaveTimeCodeKey_ = tc_key8;
	}

	for (size_t chn_i = 0; chn_i < 26; ++chn_i) {
		const int cellid = getCellID(static_cast<int>(chn_i));
		if (cellid < 0)
			continue;

		const int cry = cellid / 100000 - 1;
		if (cry < 0 || cry >= 25)
			continue;

		const int mbid = (cellid % 10000) / 1000;
		const int hg = (cellid % 1000) / 100;
		int canvasIdx = 3;
		if (mbid == 1 && hg == 1)
			canvasIdx = 0;
		else if (mbid == 1 && hg == 0)
			canvasIdx = 1;
		else if (mbid == 0 && hg == 1)
			canvasIdx = 2;

		caloWaveCanvas_[canvasIdx]->cd(cry + 1);
		auto *grWave = new TGraph();
		grWave->SetName(Form("%s_Cry%d", tag.Data(), cry + 1));
		grWave->SetTitle(Form("%s Cry%d;Sampling Point;ADC", tag.Data(), cry + 1));

		double y0 = 0.;
		for (size_t sp = 0; sp < 128; ++sp) {
			const auto d1 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2]);
			const auto d2 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2 + 1]);
			const unsigned int sp_data = static_cast<unsigned int>(d1 << 8 | d2);
			if (sp == 0)
				y0 = sp_data;
			grWave->SetPoint(static_cast<int>(sp), static_cast<double>(sp), sp_data);
		}
		grWave->SetMarkerStyle(20);
		grWave->SetMarkerSize(0.5);
		grWave->GetYaxis()->SetRangeUser(y0 - kWaveYHalfRange, y0 + kWaveYHalfRange);
		grWave->Draw("APL");
		styleWaveAxes(grWave);
		drawWaveThresholdLegend(canvasIdx, cry);
	}
	return true;
}

bool ComReader::readCsIWave(char *b){
	for(size_t chn_i=0;chn_i<8;chn_i++){ //Loop over 8 channels:
		int plat=0;
		int maxi=0;
		for(size_t sp=0;sp<128;sp++){ //Loop over 128 sampling points
			//FIXED : Data for the sampling point
			auto d1 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2]); // bit 15:8
			auto d2 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2 + 1]); // bit 7:0
			unsigned int sp_data = static_cast<unsigned int>(d1 << 8 | d2); //sampling point value
			if(sp<16){plat+=sp_data;} //Pedestal is calculated by the first 16 points
			maxi = sp_data > maxi ? sp_data : maxi; //Maximum value
		}
		plat/=16; // Calculate the pedestal value
		int cellid = (chn_i+1)*100000;
		if(cellid==-1)continue;
		CsI_CellID.emplace_back(cellid);
		CsI_CellADC.emplace_back(static_cast<unsigned int>(maxi));
		CsI_CellPLAT.emplace_back(plat);
	}
	return true;
}
bool ComReader::readCaloWave(char *b){ // TODO I decided to drop it here, little chick xuan you can add the csi wave part. Damn
	// Shi'! chick zhen bro -_-
	for(size_t chn_i=0;chn_i<26;chn_i++){ //Loop over 26 channels
		int plat=0;
		int maxi=0;
		int column=10;
		// std::cout << "=========================" << std::endl;
		for(size_t sp=0;sp<128;sp++){ //Loop over 128 sampling points
			//Data for the sampling point
			auto d1 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2]); // bit 15:8
			auto d2 = static_cast<unsigned char>(b[chn_i * 128 * 2 + sp * 2 + 1]); // bit 7:0
			unsigned int sp_data = static_cast<unsigned int>(d1 << 8 | d2); //sampling point value
			// std::cout << std::hex << sp_data << " ";
			column++;
			if (column==10)	{std::cout << std::endl;column=0;}
			if(sp<16){plat+=sp_data;} //Pedestal is calculated by the first 16 points
			maxi = sp_data > maxi ? sp_data : maxi; //Maximum value
		}
		// std::cout << std::hex << "plat:" << plat << " " << "maxi:" << maxi << std::endl;
		plat/=16; // Calculate the pedestal value
		int cellid = getCellID(chn_i);
		if(cellid==-1)continue;
		Calo_CellID.emplace_back(cellid);
		Calo_CellADC.emplace_back(static_cast<unsigned int>(maxi));
		Calo_CellPLAT.emplace_back(plat);
		//std::cout << "CellADC : CellPLAT : "<< static_cast<unsigned int>(maxi) << " " << plat << std::endl;
	}
	return true;
}
bool ComReader::readBuffer(char *b,const int& NCHN) {
	if(NCHN==26){ //Calo
		for(size_t chn_i = 0; chn_i < NCHN; chn_i++)
		{
			auto d1 = static_cast<unsigned char>(b[chn_i * 4]);
			auto d2 = static_cast<unsigned char>(b[chn_i * 4 + 1]);
			auto d3 = static_cast<unsigned char>(b[chn_i * 4 + 2]);
			auto d4 = static_cast<unsigned char>(b[chn_i * 4 + 3]);
			unsigned int plat = static_cast<unsigned int>(d1 << 8 | d2);
			unsigned int maxi = static_cast<unsigned int>(d3 << 8 | d4);
			int cellid = getCellID(chn_i);
			if(cellid==-1)continue;
			Calo_CellID.emplace_back(cellid);
			Calo_CellADC.emplace_back(maxi);
			Calo_CellPLAT.emplace_back(plat);
		}
	}
	else if(NCHN==8){//CsI
		for(int chn_i = 0; chn_i < NCHN; chn_i++)
		{
			auto d1 = static_cast<unsigned char>(b[chn_i * 4]);
			auto d2 = static_cast<unsigned char>(b[chn_i * 4 + 1]);
			auto d3 = static_cast<unsigned char>(b[chn_i * 4 + 2]);
			auto d4 = static_cast<unsigned char>(b[chn_i * 4 + 3]);
			unsigned int plat = static_cast<unsigned int>(d1 << 8 | d2);
			unsigned int maxi = static_cast<unsigned int>(d3 << 8 | d4);
			int tmp_cellid = (chn_i+1)*100000 ;
			CsI_CellID.emplace_back(tmp_cellid);
			CsI_CellADC.emplace_back(maxi);
			CsI_CellPLAT.emplace_back(plat);
		}
	}
	else{
		std::cerr << "Unknown channel number: " << NCHN << std::endl;
	}
	return true;
}

// Calo: Steinhart-Hart, exact match to MATLAB Calo_status_unpack_2v0.m
// temp1Voltage = (temp1*2500/4096*2100/2000); temp1R = 10000*temp1Voltage./(2500-temp1Voltage);
// temp01 = 2*(-86421.72414)./(-4622.53337+sqrt(4622.53337^2-4*(-86421.72414)*(-6.01188-log(temp1R))))-273.15;
static float tempFromAdc(unsigned int ti) {
    double denom = 2500.0 - 0.64087 * ti;
    if (denom <= 0) throw std::runtime_error("2500 - 0.64087*tempx must be > 0");
    double log_arg = (6408.7 * ti) / denom;
    if (log_arg <= 0) throw std::runtime_error("log argument must be > 0");
    double val = 4622.53337 * 4622.53337
               + 345686.89656 * (-6.001188 - std::log(log_arg));
    if (val < 0) throw std::runtime_error("sqrt argument must be >= 0");
    return -172843.44828 / (-4622.53337 + std::sqrt(val)) - 273.15;
}

bool ComReader::readHK(){
	char rest[60];
	file->rdbuf()->sgetn(rest,sizeof(rest));
	//for(int i=0; i<60; i++)std::cout << " " << static_cast<unsigned char>(rest[i])<< std::endl;

	const bool isCalo = (m_statusFeID >= 5 && m_statusFeID <= 8);
	const bool isCsITK = (m_statusFeID == 2);
	int c0_off, c1_off, c2_off, t_base;
	if (isCalo || isCsITK) {
		// CsITK: [2,3]=cur1P8V, [4,5]=curPam0, [6,7]=curPam1, [8..15]=temp1..4
		c0_off = 2; c1_off = 4; c2_off = 6; t_base = 8;
	} else {
		logger->error("Unknown FEEID: {}", m_statusFeID);
		return false;
	}

	auto read12 = [rest](int idx) -> unsigned int {
		auto a = static_cast<unsigned char>(rest[idx]);
		auto b = static_cast<unsigned char>(rest[idx+1]);
		return static_cast<unsigned int>((a & 0x0F) << 8 | b);
	};
	auto read16 = [rest](int idx) -> unsigned int {
		auto a = static_cast<unsigned char>(rest[idx]);
		auto b = static_cast<unsigned char>(rest[idx+1]);
		return static_cast<unsigned int>((a & 0xFF) << 8 | b);
	};
	//for (int i = 0; i < 10; i++) {std::cout << std::fixed << std::hex << read16(i) << std::endl;}

	// Current: (adc*2500/4096-245.1)/20 / scale; Calo cur1P8V=0.01, Pam0=0.125 or 0.25, Pam1=0.25; CsITK cur1P8V=0.04, Pam=0.25
	float ci0 = (read12(c0_off)*2500.f/4096.f - 245.1f)/20.f;
	float ci1 = (read12(c1_off)*2500.f/4096.f - 245.1f)/20.f;
	float ci2 = (read12(c2_off)*2500.f/4096.f - 245.1f)/20.f;
	if (isCalo) {
		ci0 /= 0.01f;
		ci1 /= (FEEID == 5 || FEEID == 6 ? 0.125f : 0.25f);
		hkC0 = isnan(ci0) ? -1.f : ci0;
		hkC1 = isnan(ci1) ? -1.f : ci1;
		hkC2 = 0.f;
	} else {
		ci0 /= 0.04f;
		ci1 /= 0.25f;
		ci2 /= 0.25f;
		hkC0 = isnan(ci0) ? -1.f : ci0;
		hkC1 = isnan(ci1) ? -1.f : ci1;
		hkC2 = isnan(ci2) ? -1.f : ci2;
	}
	

	for (int i = 1; i <= 4; i++) {
		int idx = t_base + (i-1)*2;
		unsigned int ti = read12(idx);
		//std::cout << std::fixed << std::hex << read16(idx) << " "<< ti << std::endl;
		float tif = tempFromAdc(ti);
		if (i == 1) hkT0 = isnan(tif) ? -1.f : 	tif;
		else if (i == 2) hkT1 = isnan(tif) ? -1.f : tif;
		else if (i == 3) hkT2 = isnan(tif) ? -1.f : tif;
		else hkT3 = isnan(tif) ? -1.f : tif;
	}

	// TPointID++; FeeId = m_statusFeID; // 2 / 5 / 6 / 7 / 8，供 draw_temp 按遥测 FEE 分图
	hkTree->Fill();
	return true;
}

void ComReader::initCaloTree(){ //Initialize Calo Tree
	caloTree = new TTree("caloTree","caloTree");
	caloTree->Branch("PackageCount",&Package_Count);
	caloTree->Branch("PackageLength",&Package_Length);
	caloTree->Branch("TimeCode", &Time_Code);
	caloTree->Branch("ScienceSymbol",&Science_Symbol);
	caloTree->Branch("DetectorID",&DetectorID);
	caloTree->Branch("PackageID",&PackageID);
	caloTree->Branch("TriggerStatus",&Calo_Trigger_Status);
	caloTree->Branch("TriggerID",&Calo_TriggerID);
	caloTree->Branch("CRC",&Calo_CRC);
	caloTree->Branch("Asum",&Calo_Asum);
	caloTree->Branch("EventID",&Calo_EventID);
	caloTree->Branch("CellID",&Calo_CellID);
	caloTree->Branch("CellADC",&Calo_CellADC);
	caloTree->Branch("CellPLAT",&Calo_CellPLAT);
}

void ComReader::initCsITree(){ //Initialize CsITK Tree
	csiTree = new TTree("csiTree","csiTree");
	csiTree->Branch("PackageCount",&Package_Count);
	csiTree->Branch("PackageLength",&Package_Length);
	csiTree->Branch("TimeCode", &Time_Code);
	csiTree->Branch("ScienceSymbol",&Science_Symbol);
	csiTree->Branch("DetectorID",&DetectorID);
	csiTree->Branch("PackageID",&PackageID);
	csiTree->Branch("TriggerStatus",&CsI_Trigger_Status);
	csiTree->Branch("TriggerID",&CsI_TriggerID);
	csiTree->Branch("CRC",&CsI_CRC);
	csiTree->Branch("Asum",&CsI_Asum);
	csiTree->Branch("EventID",&CsI_EventID);
	csiTree->Branch("CellID",&CsI_CellID);
	csiTree->Branch("CellADC",&CsI_CellADC);
	csiTree->Branch("CellPLAT",&CsI_CellPLAT);
}

void ComReader::initHKTree(){ //Initialize HK Tree (Calo and CsITK)
	hkTree = new TTree("hkTree"	,"hkTree");
	hkTree->Branch("FEEID",&FEEID);
	hkTree->Branch("C0",&hkC0);
	hkTree->Branch("C1",&hkC1);
	hkTree->Branch("C2",&hkC2);
	hkTree->Branch("T0",&hkT0);
	hkTree->Branch("T1",&hkT1);
	hkTree->Branch("T2",&hkT2);
	hkTree->Branch("T3",&hkT3);
}

void ComReader::decode(const std::string& filename, const std::string& yamlname){
	logger->info("Decoding file: {}", filename);
	getOutputName(std::string("result_"),filename);
    fout = new TFile(oname.c_str(),"RECREATE");
	logger->info("Initialized output file: {}", oname);
	initCaloTree();
	initCsITree();
	initHKTree();
	logger->info("Initialized output trees");
    openFile(filename);
	long long total_size = file->seekg(0, std::ios::end).tellg();
	logger->info("Total file size: {}", total_size);
	logger->info("If Wave-mode possibly, read yaml: {}", yamlname);
	config = YAML::LoadFile(yamlname);
	yaml_package_mode_id = config["package_mode_id"].as<unsigned int>();
	yaml_draw_wave_events_.clear();
	yaml_calo_trigger_id.clear();
	yaml_csi_trigger_id.clear();
	yaml_calo_draw_time_code_.clear();
	yaml_csi_draw_time_code_.clear();
	auto loadDrawWaveFromArrays = [&](const YAML::Node& triggers_node,
	                                  const YAML::Node& time_codes_node) {
		const auto triggers = triggers_node.as<std::vector<unsigned int>>();
		const auto time_codes = time_codes_node.as<std::vector<long long>>();
		const size_t n = std::min(triggers.size(), time_codes.size());
		if (triggers.size() != time_codes.size())
			logger->warn("draw_wave trigger_id/time_code length mismatch: {} vs {}, using first {}",
			             triggers.size(), time_codes.size(), n);
		for (size_t i = 0; i < n; ++i)
			yaml_draw_wave_events_.push_back({triggers[i], time_codes[i]});
	};
	if (config["draw_wave_trigger_id"] && config["draw_wave_time_code"]) {
		loadDrawWaveFromArrays(config["draw_wave_trigger_id"], config["draw_wave_time_code"]);
	} else if (config["draw_wave_events"] && config["draw_wave_events"].IsMap() &&
	           config["draw_wave_events"]["trigger_id"] &&
	           config["draw_wave_events"]["time_code"]) {
		loadDrawWaveFromArrays(config["draw_wave_events"]["trigger_id"],
		                       config["draw_wave_events"]["time_code"]);
	} else if (config["draw_wave_events"]) {
		for (const auto& item : config["draw_wave_events"]) {
			if (!item.IsSequence() || item.size() < 2)
				continue;
			WaveDrawEvent e;
			e.trigger_id = item[0].as<unsigned int>();
			e.time_code_key8 = item[1].as<long long>();
			yaml_draw_wave_events_.push_back(e);
		}
	}
	if (config["csi_draw_trigger_id"])
		yaml_csi_trigger_id = config["csi_draw_trigger_id"].as<std::vector<unsigned int>>();
	if (config["calo_draw_trigger_id"])
		yaml_calo_trigger_id = config["calo_draw_trigger_id"].as<std::vector<unsigned int>>();
	if (config["csi_draw_time_code"])
		yaml_csi_draw_time_code_ = config["csi_draw_time_code"].as<std::vector<long long>>();
	if (config["calo_draw_time_code"])
		yaml_calo_draw_time_code_ = config["calo_draw_time_code"].as<std::vector<long long>>();
	if (!yaml_draw_wave_events_.empty())
		logger->info("Wave draw by (TriggerID, TimeCode key8) pairs: {} entries",
		             yaml_draw_wave_events_.size());
	else if (!yaml_csi_draw_time_code_.empty() || !yaml_calo_draw_time_code_.empty())
		logger->info("Wave draw fallback: CsI TimeCode key8={}, Calo TimeCode key8={}",
		             yaml_csi_draw_time_code_.size(), yaml_calo_draw_time_code_.size());
	else if (!yaml_csi_trigger_id.empty() || !yaml_calo_trigger_id.empty())
		logger->info("Wave draw fallback: CsI TriggerID={}, Calo TriggerID={}",
		             yaml_csi_trigger_id.size(), yaml_calo_trigger_id.size());
	time_code_pass_valid_ = false;
	time_code_calo_pass_valid_ = false;
	time_code_csi_pass_valid_ = false;
	yaml_time_filter_enable = false;
	yaml_csi_first_cross_hist_enable_ = false;
	csi_first_cross_late_.clear();
	for (int i = 0; i < 8; ++i) {
		delete csi_first_cross_hist_[i];
		csi_first_cross_hist_[i] = nullptr;
	}
	if (config["csi_first_cross_hist"])
		yaml_csi_first_cross_hist_enable_ = config["csi_first_cross_hist"].as<bool>();
	if (yaml_csi_first_cross_hist_enable_) {
		logger->info("CsI first-cross histogram enabled (first over-thr + {} following, 8 chn, 128 bins)",
		             kCsIFirstCrossRunLength - 1);
		for (int i = 0; i < 8; ++i) {
			csi_first_cross_hist_[i] = new TH1F(
			    Form("hCsI_first_cross_chn%d", i),
			    Form("CsI Chn%d first over-thr + %d following (thr=%d);Start index;Entries", i,
			         kCsIFirstCrossRunLength - 1, kCsIWaveThreshold[i]),
			    128, 0, 128);
			csi_first_cross_hist_[i]->SetDirectory(nullptr);
		}
	}
	if (config["time_code_min"] && config["time_code_max"]) {
		yaml_time_code_min = config["time_code_min"].as<double>();
		yaml_time_code_max = config["time_code_max"].as<double>();
		yaml_time_filter_enable = true;
		logger->info("TimeCode filter: [{}, {}] (base BJT 2025-01-01 00:00:00)", yaml_time_code_min,
		             yaml_time_code_max);
	}
	file->seekg(0, std::ios::beg);
	fout->mkdir("CsIWavePlot"); fout->mkdir("CaloWavePlot");
    while(true){
        if(!findHead()){
			break;
		}
    }

	if (yaml_time_filter_enable)
		logger->info("Skipped by time filter: {}", nSkippedTime);
	logTimeCodePassRange();
	flushCaloWavePlots();
	writeCsIFirstCrossHist();
	logger->info("Total head: {}", nhead);
	logger->info("CALO: {} expected Head{}", Calo_EventID, nCaloHead/4.);
	logger->info("CsITK: {} expected Head {}", CsI_EventID, nCsIHead);
	// logger->info("HK: {} rows (packets), nHKHead {}", TPointID, nHKHead);
	logger->info("HK Calo: {} events", nHKCaloHead/4);
	logger->info("HK CsITK: {} events", nHKCsIHead);
    fout->cd();
	//build a directory to store the waveforms,
    caloTree->Write();
	csiTree->Write();
	hkTree->Write();
    fout->Close();
}

void ComReader::getOutputName(const std::string& prefix,const std::string& filename){
    oname = filename.substr(filename.find_last_of("/")+1);
    oname = oname.substr(0,oname.find_last_of("."));
    oname = prefix+oname+".root";
}

ComReader::~ComReader(){
    file->close();
    // delete file;
    
}
