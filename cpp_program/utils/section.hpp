#ifndef section_h
#define section_h

#include <string>
#include <iostream>

#include "timer.hpp"
#include "settings.hpp"

void print_hLine(std::ostream & os = std::cout){
	os << "----------------------------------------\n";
};

struct Section {
    std::string section_name;
    Timer timer;
    int min_verbose_level;
    Section(){};
    Section(std::string section_name, int min_verbose_level = 0, std::ostream & os = std::cout) : section_name(section_name), min_verbose_level(min_verbose_level) {
        if(settings.VERBOSE_LEVEL < min_verbose_level) return;
        os << "\n----------------------------------------\n";
        os << "Starting: " << section_name << "\n";
        os << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
        timer.start_timer();
    };
    void end_section(std::ostream & os = std::cout){
        if(settings.VERBOSE_LEVEL < min_verbose_level) return;
        timer.stop_timer();
        os << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        os << "Ending: " << section_name << "\n";
        os << "took: " << timer.getTimeInSec() << " sec" << "\n";
        os << "----------------------------------------\n\n";
    };
};

#endif