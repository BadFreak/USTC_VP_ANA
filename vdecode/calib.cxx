#include "Caliber.hh"

int main(int argc,char** argv) {
    auto main_logger = spdlog::stdout_color_mt("main");
    main_logger->info("Welcome to VDECODE CALIBER!");
    auto caliber = new Caliber();
    caliber->calib(std::string("calo"));
	caliber->calib(std::string("csi"));
    return 0;
}