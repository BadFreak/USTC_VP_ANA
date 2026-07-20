#include "Plotter.hh"

int main(int argc,char** argv) {
    auto main_logger = spdlog::stdout_color_mt("main");
    main_logger->info("Welcome to VDECODE PLOTTER!");
    if (argc != 2) {
        main_logger->error("Usage: {} <FILENAME>", argv[0]);
        return 1;
    }
    auto plotter = new Plotter();
    plotter->draw(1., argv[1]);
    return 1;
}