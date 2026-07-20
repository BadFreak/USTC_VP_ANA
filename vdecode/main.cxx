// #include "WaveReader.hh"
#include "ComReader.hh"

int main(int argc,char** argv) {
    auto main_logger = spdlog::stdout_color_mt("main");
    main_logger->info("Welcome to VDECODE!");
    if (argc != 4) {
        main_logger->error("Usage: {} <mode> <filename> <yaml>", argv[0]);
        return 1;
    }
    if(std::string(argv[1]) == "combine"){
        main_logger->info("Combining files...");
        ComReader::getInstance().decode(std::string(argv[2]), std::string(argv[3]));
    }
    else{
        main_logger->error("Unknown mode: {}", argv[1]);
    }
    
    return 1;
}
