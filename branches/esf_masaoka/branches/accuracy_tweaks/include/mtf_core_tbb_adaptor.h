#ifndef MTF_CORE_TBB_ADAPTOR
#define MTF_CORE_TBB_ADAPTOR

#include "include/mtf_core.h"
#include "include/point_helpers.h"

class Mtf_core_tbb_adaptor {
  public:
    Mtf_core_tbb_adaptor(Mtf_core* core) : mtf_core(core) {
    }

    void operator()(const blocked_range<size_t>& r) const {
        for (size_t i=r.begin(); i != r.end(); ++i) {
            Boundarylist::const_iterator it = mtf_core->cl.get_boundaries().find(mtf_core->valid_obj[i]);
            Point cent = centroid(it->second);
            mtf_core->search_borders(cent, mtf_core->valid_obj[i]);
        }
    }
  
    Mtf_core* mtf_core;  
};

#endif
