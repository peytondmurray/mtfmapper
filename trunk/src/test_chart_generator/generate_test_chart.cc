#include "svg_page.h"
#include "svg_page_grid.h"
#include "svg_page_perspective.h"
#include "config.h"

#include <tclap/CmdLine.h>

#include <string>
using std::string;
using std::stringstream;

int main(int argc, char** argv) {

    stringstream ss;
    ss << mtfmapper_VERSION_MAJOR << "." << mtfmapper_VERSION_MINOR;
    
    vector<string> allowed_types;
    allowed_types.push_back("perspective");
    allowed_types.push_back("grid");
    TCLAP::ValuesConstraint<string> type_constraints(allowed_types);
    
    vector<string> allowed_sizes;
    allowed_sizes.push_back("a4");
    allowed_sizes.push_back("A4");
    allowed_sizes.push_back("a3");
    allowed_sizes.push_back("A3");
    allowed_sizes.push_back("a2");
    allowed_sizes.push_back("A2");
    TCLAP::ValuesConstraint<string> size_constraints(allowed_sizes);
    
    TCLAP::CmdLine cmd("Generate test charts for MTF50 measurements", ' ', ss.str());
    TCLAP::ValueArg<std::string> tc_type("t", "type", "Chart type (currently \"grid\" or \"perspective\")", false, "perspective", &type_constraints );
    cmd.add(tc_type);
    TCLAP::ValueArg<std::string> tc_size("s", "size", "Chart size (currently \"a4\", \"a3\" or \"a2\")", false, "a3", &size_constraints );
    cmd.add(tc_size);
    TCLAP::ValueArg<std::string> tc_ofname("o", "output", "Output file name (default chart.svg)", false, "chart.svg", "string" );
    cmd.add(tc_ofname);
    
    cmd.parse(argc, argv);
    
    printf("char type: %s at %s size\n", tc_type.getValue().c_str(), tc_size.getValue().c_str());
    
    if ( tc_type.getValue().compare("perspective") == 0) {
        Svg_page_perspective p(tc_size.getValue(), tc_ofname.getValue());
    
        //p.set_viewing_parameters(20, 5, 45/180.0*M_PI); // must still become paramers
        p.render();
    } else
    if ( tc_type.getValue().compare("grid") == 0) {
        Svg_page_grid p(tc_size.getValue(), tc_ofname.getValue());
        p.render();
    } else {
        printf("Illegal char type %s\n", tc_type.getValue().c_str());
    }
    
    return 0;
}
